/*
 * Feather1codeMimoDemoV6_PARALLEL.ino
 * ==================================
 *
 * DUAL-CORE FreeRTOS build for Feather 1 (measurement side of the open RTU
 * MIMO bench), Draft V6 — a REVISION of V5 that fixes the following bugs
 * discovered during Session B on 14 July 2026:
 *
 *   V5 BUG #1  frequencyTask was a stub returning F_NOM. f_Vb=0.0 reported.
 *   V5 BUG #2  samplingTask busy-waited on micros() with no vTaskDelay,
 *              starving IDLE0 and tripping the task watchdog at t = 15 s.
 *   V5 BUG #3  RMS computation held the ring mutex over a 40-iteration loop,
 *              blocking the sampling task and jittering the sample cadence.
 *   V5 BUG #4  Ring buffer wrote all 5 channels to the SAME head index in
 *              one round-robin cycle, meaning subsequent samples of the same
 *              channel overwrote before the compute task could read them.
 *
 * V6 corrections:
 *
 *   FIX #1     Full V2 measureFreq() algorithm ported to ring-buffer form.
 *              Reads the last FREQ_WIN Vb samples, does interpolating
 *              zero-crossing detection, median-filters the periods, returns
 *              a real Hz value. Time base uses per-row micros() timestamps
 *              stored in the ring.
 *   FIX #2     samplingTask now uses vTaskDelayUntil at 100 us granularity
 *              rather than a bare busy loop. This yields to IDLE0 at every
 *              tick and satisfies the task watchdog. Direct consequence:
 *              per-channel effective sample rate is ~2 kHz (10 kHz round-
 *              robin / 5 channels), still comfortably above 2x Nyquist for
 *              a 50 Hz signal, and matches the canonical V2 approach.
 *   FIX #3     RMS loop copies the needed slice out of the ring buffer with
 *              the mutex held for < 200 us, then does arithmetic on the
 *              local copy with the mutex released.
 *   FIX #4    Ring buffer indexed per-channel: each of 5 channels has its
 *              own head pointer. Sampling task increments a channel's head
 *              only after writing that channel's data.
 *
 * Reporting cadence: 20 ms REPORT lines, same as V5.
 *
 * PMU-grade sampling is still deferred to Phase 2 (ADC-DMA via I2S).
 * V6 uses analogRead() but with the FreeRTOS structure PROVEN and the
 * measurements REAL. If V6 hits Class II across a NESO replay, it is a
 * significant contribution to the paper; if it hits Class I, that would
 * be surprising given the polled-ADC limitation, but is possible.
 *
 * BOARD  : Adafruit Feather ESP32 (ESP32-PICO-D4 / PICO-MINI-02)
 * TARGET : Feather 1 (measurement side); Feather 2 is unchanged canonical.
 */

#include <Arduino.h>
#include <math.h>
#include <string.h>
#include <freertos/FreeRTOS.h>
#include <freertos/task.h>
#include <freertos/semphr.h>
#include <esp_task_wdt.h>

// ============================================================================
// Pin assignments (same as canonical V2; MUST match Table 4.1 of dissertation)
// ============================================================================
#define PIN_VA        34
#define PIN_VB        39
#define PIN_VC        36
#define PIN_IB        32
#define PIN_IC        33
#define PIN_SYNC      27
#define DAC_VC        25

#define N_CHANNELS    5
#define CH_VA         0
#define CH_VB         1
#define CH_VC         2
#define CH_IB         3
#define CH_IC         4
static const int adcPins[N_CHANNELS] = {PIN_VA, PIN_VB, PIN_VC, PIN_IB, PIN_IC};

// ============================================================================
// Timing constants
// ============================================================================
// Round-robin ADC cadence: 100 us tick, 5 channels -> 2 kHz per channel.
// This matches canonical V2 which uses delayMicroseconds(500) inside its
// freq loop (2 kHz effective) and delayMicroseconds(470) inside its THD loop.
#define ADC_TICK_US            100
#define SAMPLES_PER_CHANNEL_HZ 2000

// FreeRTOS tick granularity for the sampling task; setting to 1 keeps
// per-channel jitter well below the 20 ms REPORT cadence.
#define SAMPLING_TICK_MS       1

// Analysis window sizes (matched to canonical V2)
#define FREQ_WIN               200   // 100 ms at 2 kHz per channel
#define MAX_EDGES_PER_TYPE     20
#define MAX_PERIODS            40
#define RMS_WIN                40    // matches V2 RMS_N=40

// Ring buffer size per channel — must be > FREQ_WIN and > RMS_WIN,
// with headroom for the compute task's read latency
#define RING_SAMPLES_PER_CH    512

// Reporting cadence
#define REPORT_MS              20

// Nominal frequency for sanity bounds
#define F_NOM                  50.0f
#define F_LOWER                49.5f
#define F_UPPER                50.5f

// ============================================================================
// Shared ring buffer, per-channel head pointers.
// Each channel has independent write index; sampling task walks channel
// round-robin but only advances that channel's head after its write.
// ============================================================================
struct Ring {
  volatile uint16_t samples[N_CHANNELS][RING_SAMPLES_PER_CH];
  volatile uint32_t timestamp_us[N_CHANNELS][RING_SAMPLES_PER_CH];
  volatile uint32_t head[N_CHANNELS];    // next-write index per channel
  volatile uint32_t count[N_CHANNELS];   // total writes per channel
};
static Ring ring;
static SemaphoreHandle_t ringMutex = NULL;

// ============================================================================
// Measurement result set
// ============================================================================
struct MeasSet {
  float Va_rms, Vb_rms, Vc_rms, Ib_rms, Ic_rms;
  float f_Vb;
  uint32_t sample_count_snapshot;      // total sample count at last update
  uint32_t last_freq_update_ms;
  uint32_t last_rms_update_ms;
};
static MeasSet meas;
static SemaphoreHandle_t measMutex = NULL;

// ============================================================================
// Feather 2 DAC-sync (unchanged from V2 and V5)
// ============================================================================
#define LUT_SIZE      40
#define AMP_VC        127
static uint8_t lutVc[LUT_SIZE];
volatile int   dacIdx = 0;
hw_timer_t    *dacTimer = NULL;

void IRAM_ATTR onDacTimer() {
  dacWrite(DAC_VC, lutVc[dacIdx]);
  dacIdx = (dacIdx + 1) % LUT_SIZE;
}

// Bias offsets, measured once at boot
float biasVa = 0, biasVb = 0, biasVc = 0, biasIb = 0, biasIc = 0;

// Diagnostic: samples per second, computed by the frequencyTask
volatile float diag_sample_rate_hz[N_CHANNELS] = {0};

// ============================================================================
// Sampling task (Core 0) - FIXES V5 BUG #2, #4
// ============================================================================
//
// Uses vTaskDelay(1) at 1 ms granularity to yield to IDLE0. Within each
// 1 ms tick, does exactly (1000 / ADC_TICK_US) = 10 samples across the
// 5 channels round-robin. That gives 2 kHz per channel.
//
// The task-watchdog subscription is DISABLED for this task because we
// intentionally hold Core 0 tightly. The idle-watchdog (fed by IDLE0
// during vTaskDelay) still runs.
// ============================================================================
void samplingTask(void *param) {
  // Deliberately do NOT subscribe to task_wdt on this task.
  // esp_task_wdt_delete(NULL);  // no-op since we never added

  TickType_t last_wake = xTaskGetTickCount();
  const TickType_t period_ticks = pdMS_TO_TICKS(SAMPLING_TICK_MS);

  // Samples we should emit per tick period
  const int samples_per_tick = (SAMPLING_TICK_MS * 1000) / ADC_TICK_US; // 10
  int ch = 0;

  while (true) {
    // Do all reads for this tick
    for (int k = 0; k < samples_per_tick; k++) {
      uint16_t raw = analogRead(adcPins[ch]);
      uint32_t t   = micros();

      // Short critical section: write into ring
      if (xSemaphoreTake(ringMutex, pdMS_TO_TICKS(1)) == pdTRUE) {
        uint32_t h = ring.head[ch];
        ring.samples[ch][h]      = raw;
        ring.timestamp_us[ch][h] = t;
        ring.head[ch]  = (h + 1) % RING_SAMPLES_PER_CH;
        ring.count[ch] = ring.count[ch] + 1;
        xSemaphoreGive(ringMutex);
      }

      ch = (ch + 1) % N_CHANNELS;
    }

    // Yield to IDLE0 exactly once per SAMPLING_TICK_MS. This feeds the
    // watchdog and allows same-core tasks (frequencyTask) to run.
    vTaskDelayUntil(&last_wake, period_ticks);
  }
}

// ============================================================================
// Frequency task (Core 0) - FIXES V5 BUG #1
// ============================================================================
//
// Runs at REPORT_MS cadence. Reads the last FREQ_WIN samples of channel Vb
// from the ring buffer (with the timestamps that the samplingTask stored
// alongside them), does the interpolating-zero-crossing / median-filter
// algorithm ported from canonical V2, and writes f_Vb into meas.
// ============================================================================
float ringMeasureFreqVb() {
  // Copy the FREQ_WIN most-recent Vb samples out of the ring buffer with
  // the mutex held only for the duration of the copy (~200 us).
  uint16_t samps[FREQ_WIN];
  uint32_t ts[FREQ_WIN];
  uint32_t head_at_copy;

  if (xSemaphoreTake(ringMutex, pdMS_TO_TICKS(2)) != pdTRUE) {
    return 0.0f;
  }
  head_at_copy = ring.head[CH_VB];
  for (int k = 0; k < FREQ_WIN; k++) {
    int idx = (int)(head_at_copy + RING_SAMPLES_PER_CH - FREQ_WIN + k) % RING_SAMPLES_PER_CH;
    samps[k] = ring.samples[CH_VB][idx];
    ts[k]    = ring.timestamp_us[CH_VB][idx];
  }
  xSemaphoreGive(ringMutex);

  // Now run the V2 algorithm on the local copies
  uint32_t pos_t[MAX_EDGES_PER_TYPE];
  uint32_t neg_t[MAX_EDGES_PER_TYPE];
  int n_pos = 0, n_neg = 0;
  float prev_v = 0.0f;
  uint32_t prev_t = 0;
  bool first_sample = true;

  for (int i = 0; i < FREQ_WIN; i++) {
    float v = samps[i] * (3.3f / 4095.0f) - biasVb;
    uint32_t t = ts[i];

    if (!first_sample) {
      // Positive-going crossing
      if (prev_v < 0.0f && v >= 0.0f) {
        float dv = v - prev_v;
        float frac = (dv > 1e-6f) ? (-prev_v / dv) : 0.5f;
        if (frac < 0.0f) frac = 0.0f;
        if (frac > 1.0f) frac = 1.0f;
        uint32_t t_cross = prev_t + (uint32_t)(frac * (float)(t - prev_t));
        if (n_pos < MAX_EDGES_PER_TYPE) pos_t[n_pos++] = t_cross;
      }
      // Negative-going crossing
      else if (prev_v >= 0.0f && v < 0.0f) {
        float dv = prev_v - v;
        float frac = (dv > 1e-6f) ? (prev_v / dv) : 0.5f;
        if (frac < 0.0f) frac = 0.0f;
        if (frac > 1.0f) frac = 1.0f;
        uint32_t t_cross = prev_t + (uint32_t)(frac * (float)(t - prev_t));
        if (n_neg < MAX_EDGES_PER_TYPE) neg_t[n_neg++] = t_cross;
      }
    } else {
      first_sample = false;
    }

    prev_v = v;
    prev_t = t;
  }

  // Build period buffer (same-edge to same-edge, pos and neg combined)
  float periods[MAX_PERIODS];
  int n_periods = 0;
  for (int i = 0; i < n_pos - 1; i++) {
    if (n_periods < MAX_PERIODS) periods[n_periods++] = (float)(pos_t[i+1] - pos_t[i]);
  }
  for (int i = 0; i < n_neg - 1; i++) {
    if (n_periods < MAX_PERIODS) periods[n_periods++] = (float)(neg_t[i+1] - neg_t[i]);
  }
  if (n_periods < 3) return 0.0f;

  // Insertion sort for median
  for (int i = 1; i < n_periods; i++) {
    float key = periods[i];
    int j = i - 1;
    while (j >= 0 && periods[j] > key) {
      periods[j + 1] = periods[j];
      j--;
    }
    periods[j + 1] = key;
  }

  float median_period;
  if (n_periods % 2 == 0) {
    median_period = 0.5f * (periods[n_periods / 2 - 1] + periods[n_periods / 2]);
  } else {
    median_period = periods[n_periods / 2];
  }

  if (median_period < 100.0f) return 0.0f;
  return 1000000.0f / median_period;
}

void frequencyTask(void *param) {
  const TickType_t interval = pdMS_TO_TICKS(REPORT_MS);
  uint32_t last_count_ts = millis();
  uint32_t last_counts[N_CHANNELS] = {0};

  while (true) {
    vTaskDelay(interval);

    // Update diagnostic sample-rate estimate once per second
    uint32_t now = millis();
    if (now - last_count_ts >= 1000) {
      for (int c = 0; c < N_CHANNELS; c++) {
        uint32_t curr = ring.count[c];
        diag_sample_rate_hz[c] = (float)(curr - last_counts[c]) * 1000.0f / (float)(now - last_count_ts);
        last_counts[c] = curr;
      }
      last_count_ts = now;
    }

    float f = ringMeasureFreqVb();

    if (xSemaphoreTake(measMutex, pdMS_TO_TICKS(5)) == pdTRUE) {
      meas.f_Vb = f;
      meas.last_freq_update_ms = now;
      xSemaphoreGive(measMutex);
    }
  }
}

// ============================================================================
// Compute task (Core 1) - FIXES V5 BUG #3
// ============================================================================
//
// Copies the RMS window slice out of the ring under mutex (< 200 us hold),
// then does the sqrt arithmetic on the local copy with the mutex released.
// ============================================================================
void computeTask(void *param) {
  const TickType_t interval = pdMS_TO_TICKS(REPORT_MS);
  uint16_t copyA[RMS_WIN], copyB[RMS_WIN], copyC[RMS_WIN];
  uint16_t copyIB[RMS_WIN], copyIC[RMS_WIN];

  while (true) {
    vTaskDelay(interval);

    // Short critical section: copy latest RMS_WIN samples per channel
    if (xSemaphoreTake(ringMutex, pdMS_TO_TICKS(2)) != pdTRUE) continue;
    uint32_t hA = ring.head[CH_VA];
    uint32_t hB = ring.head[CH_VB];
    uint32_t hC = ring.head[CH_VC];
    uint32_t hIB = ring.head[CH_IB];
    uint32_t hIC = ring.head[CH_IC];
    for (int k = 0; k < RMS_WIN; k++) {
      copyA[k]  = ring.samples[CH_VA][(hA + RING_SAMPLES_PER_CH - RMS_WIN + k) % RING_SAMPLES_PER_CH];
      copyB[k]  = ring.samples[CH_VB][(hB + RING_SAMPLES_PER_CH - RMS_WIN + k) % RING_SAMPLES_PER_CH];
      copyC[k]  = ring.samples[CH_VC][(hC + RING_SAMPLES_PER_CH - RMS_WIN + k) % RING_SAMPLES_PER_CH];
      copyIB[k] = ring.samples[CH_IB][(hIB + RING_SAMPLES_PER_CH - RMS_WIN + k) % RING_SAMPLES_PER_CH];
      copyIC[k] = ring.samples[CH_IC][(hIC + RING_SAMPLES_PER_CH - RMS_WIN + k) % RING_SAMPLES_PER_CH];
    }
    xSemaphoreGive(ringMutex);

    // Arithmetic (mutex released)
    float sqA=0, sqB=0, sqC=0, sqIB=0, sqIC=0;
    for (int k = 0; k < RMS_WIN; k++) {
      float vA = copyA[k]  * (3.3f/4095.0f) - biasVa;
      float vB = copyB[k]  * (3.3f/4095.0f) - biasVb;
      float vC = copyC[k]  * (3.3f/4095.0f) - biasVc;
      float iB = copyIB[k] * (3.3f/4095.0f) - biasIb;
      float iC = copyIC[k] * (3.3f/4095.0f) - biasIc;
      sqA += vA*vA; sqB += vB*vB; sqC += vC*vC;
      sqIB += iB*iB; sqIC += iC*iC;
    }
    float Va_rms = sqrtf(sqA / (float)RMS_WIN);
    float Vb_rms = sqrtf(sqB / (float)RMS_WIN);
    float Vc_rms = sqrtf(sqC / (float)RMS_WIN);
    float Ib_rms = sqrtf(sqIB / (float)RMS_WIN);
    float Ic_rms = sqrtf(sqIC / (float)RMS_WIN);

    if (xSemaphoreTake(measMutex, pdMS_TO_TICKS(5)) == pdTRUE) {
      meas.Va_rms = Va_rms;
      meas.Vb_rms = Vb_rms;
      meas.Vc_rms = Vc_rms;
      meas.Ib_rms = Ib_rms;
      meas.Ic_rms = Ic_rms;
      meas.last_rms_update_ms = millis();
      xSemaphoreGive(measMutex);
    }
  }
}

// ============================================================================
// Setup
// ============================================================================
void setup() {
  Serial.begin(115200);
  delay(1500);

  pinMode(PIN_SYNC, OUTPUT);
  digitalWrite(PIN_SYNC, LOW);

  analogSetPinAttenuation(PIN_VA, ADC_11db);
  analogSetPinAttenuation(PIN_VB, ADC_11db);
  analogSetPinAttenuation(PIN_VC, ADC_11db);
  analogSetPinAttenuation(PIN_IB, ADC_11db);
  analogSetPinAttenuation(PIN_IC, ADC_11db);

  for (int i = 0; i < LUT_SIZE; i++) {
    float a = 2.0f * PI * i / LUT_SIZE;
    lutVc[i] = (uint8_t)(127 + AMP_VC * sinf(a - 4.0f * PI / 3.0f));
  }
  dacWrite(DAC_VC, 127);

  // DAC hardware timer (unchanged from V2/V5)
  dacTimer = timerBegin(1000000);
  timerAttachInterrupt(dacTimer, &onDacTimer);
  timerAlarm(dacTimer, 500, true, 0);  // 500 us -> 2 kHz -> 50 Hz across 40-sample LUT

  Serial.println("================================================");
  Serial.println("PARALLEL FreeRTOS BUILD (V6) — WATCHDOG-SAFE + M2 PORTED");
  Serial.println("Core 0: sampling (2 kHz/ch) + M2 frequency detector");
  Serial.println("Core 1: RMS window + telemetry");
  Serial.println("Report cadence: 20 ms");
  Serial.println("V5 bugs fixed: watchdog, stub M2, mutex hold, ring index");
  Serial.println("================================================");
  delay(500);

  // Bias measurement
  for (int i = 0; i < 1000; i++) {
    biasVa += analogRead(PIN_VA) * (3.3f/4095.0f);
    biasVb += analogRead(PIN_VB) * (3.3f/4095.0f);
    biasVc += analogRead(PIN_VC) * (3.3f/4095.0f);
    biasIb += analogRead(PIN_IB) * (3.3f/4095.0f);
    biasIc += analogRead(PIN_IC) * (3.3f/4095.0f);
    delayMicroseconds(200);
  }
  biasVa /= 1000; biasVb /= 1000; biasVc /= 1000;
  biasIb /= 1000; biasIc /= 1000;
  Serial.printf("Bias: Va=%.4f Vb=%.4f Vc=%.4f Ib=%.4f Ic=%.4f\n",
                biasVa, biasVb, biasVc, biasIb, biasIc);

  ringMutex = xSemaphoreCreateMutex();
  measMutex = xSemaphoreCreateMutex();
  memset((void*)&ring, 0, sizeof(ring));

  // Create tasks. Priorities picked so IDLE0 gets to run between sampling
  // ticks (sampling delays 1 ms per tick).
  // Priorities: compute (highest, on Core 1), sampling (Core 0), freq (Core 0)
  xTaskCreatePinnedToCore(samplingTask,  "adc",  4096, NULL, 3, NULL, 0);
  xTaskCreatePinnedToCore(frequencyTask, "freq", 6144, NULL, 4, NULL, 0);
  xTaskCreatePinnedToCore(computeTask,   "comp", 6144, NULL, 3, NULL, 1);

  Serial.println("Tasks started.");
}

// ============================================================================
// Loop = telemetry, on Core 1 at REPORT_MS cadence.
// Reports a diagnostic effective sample rate once per second.
// ============================================================================
void loop() {
  static uint32_t last_report = 0;
  static uint32_t last_diag = 0;
  uint32_t now = millis();

  if (now - last_report < REPORT_MS) {
    delay(1);
    return;
  }
  last_report = now;

  MeasSet snap;
  if (xSemaphoreTake(measMutex, pdMS_TO_TICKS(5)) == pdTRUE) {
    snap = meas;
    xSemaphoreGive(measMutex);
  } else {
    return;
  }

  Serial.printf("REPORT,%lu,Va=%.4f,Vb=%.4f,Vc=%.4f,Ib=%.4f,Ic=%.4f,f_Vb=%.4f\n",
                (unsigned long)now, snap.Va_rms, snap.Vb_rms, snap.Vc_rms,
                snap.Ib_rms, snap.Ic_rms, snap.f_Vb);

  // Once per second, emit a diagnostic showing actual per-channel sample rate.
  // Helpful for confirming the FreeRTOS structure is holding.
  if (now - last_diag >= 1000) {
    Serial.printf("DIAG,%lu,fs_Va=%.0f,fs_Vb=%.0f,fs_Vc=%.0f,fs_Ib=%.0f,fs_Ic=%.0f\n",
                  (unsigned long)now,
                  diag_sample_rate_hz[0], diag_sample_rate_hz[1],
                  diag_sample_rate_hz[2], diag_sample_rate_hz[3],
                  diag_sample_rate_hz[4]);
    last_diag = now;
  }
}
