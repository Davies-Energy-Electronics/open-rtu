/* ============================================================================
 * characterise.ino  —  ADC chain characterisation for the open-RTU
 * ----------------------------------------------------------------------------
 * SEPARATE BUILD. Do NOT merge any of this into Feather1codeMimoDemoV2.ino or
 * Feather1codeMimoDemoV6_PARALLEL.ino. Those two are the provenance of every
 * archived capture and must stay byte-identical to the released sources.
 *
 * Purpose: dump raw ADC codes so that the input-referred noise of the
 * conditioning chain (sigma_v) and the IEEE 1241 SINAD / ENOB of the
 * converter can be measured rather than inferred.
 *
 * Acquisition matches the canonical V2 estimator: same pins, same default
 * 12-bit / 11 dB attenuation, same 2 kHz per-channel rate. Two modes:
 *
 *   MODE_SINGLE  - Vb only, no channel switching.
 *   MODE_ROUNDROBIN - all five channels in sequence at 100 us, Vb retained.
 *                     This is what V2 and V6 actually do.
 *
 * Running both and differencing isolates multiplexer settling / charge
 * injection from the rest of the noise. That is a firmware-fixable term if
 * it turns out to dominate, so it is worth separating.
 *
 * Acquisition is silent — nothing is transmitted until the buffer is full,
 * because USB traffic during acquisition perturbs the supply rail and you
 * would be measuring the serial bridge.
 *
 * Output on the serial monitor at 921600 baud:
 *     # header lines beginning with '#'
 *     one decimal ADC code per line
 *
 * Save the whole lot to a .txt and feed it to analyse_noise.py or
 * analyse_sinad.py.
 *
 * Jack Davies, September 2026
 * ========================================================================== */

// ---- pins, identical to the measurement firmware -------------------------
#define PIN_VA   34
#define PIN_VB   39
#define PIN_VC   36
#define PIN_IB   32
#define PIN_IC   33

static const int adcPins[5] = {PIN_VA, PIN_VB, PIN_VC, PIN_IB, PIN_IC};
#define CH_VB    1

// ---- acquisition --------------------------------------------------------
#define CHAR_N            16384      // 8.192 s at 2 kHz. Power of two for the FFT.
#define SAMPLE_PERIOD_US  500        // 2 kHz per channel
#define RR_TICK_US        100        // 5 channels x 100 us = 500 us per channel

// Set ONE of these before flashing.
#define MODE_SINGLE       1
#define MODE_ROUNDROBIN   0

// Free-text label written into the header. CHANGE IT for every run.
#define RUN_LABEL  "A1_grounded_single"

static uint16_t buf[CHAR_N];

// Running accumulators for the other four channels in round-robin mode.
// Cheap, and they show at a glance which inputs are floating: a floating
// ESP32 ADC pin reads noisy and drifts, and if it does it is a candidate
// source for multiplexer charge injection into Vb.
static double chSum[5], chSumSq[5];

// ---------------------------------------------------------------------------
void acquire(uint32_t *t_start, uint32_t *t_end) {
  uint32_t next;
  for (int c = 0; c < 5; c++) { chSum[c] = 0.0; chSumSq[c] = 0.0; }

  // Discard the first few conversions: the first analogRead after a mode
  // change is not representative.
  for (int i = 0; i < 32; i++) (void)analogRead(PIN_VB);

  next = micros();
  *t_start = next;

  for (uint32_t i = 0; i < CHAR_N; i++) {
#if MODE_ROUNDROBIN
    // Walk all five channels exactly as samplingTask does, keep Vb.
    for (int ch = 0; ch < 5; ch++) {
      uint16_t r = analogRead(adcPins[ch]);
      if (ch == CH_VB) buf[i] = r;
      chSum[ch]   += r;
      chSumSq[ch] += (double)r * r;
      next += RR_TICK_US;
      while ((int32_t)(micros() - next) < 0) { /* spin */ }
    }
#else
    buf[i] = analogRead(PIN_VB);
    next += SAMPLE_PERIOD_US;
    while ((int32_t)(micros() - next) < 0) { /* spin */ }
#endif
  }
  *t_end = micros();
}

// ---------------------------------------------------------------------------
void setup() {
  Serial.begin(921600);
  delay(3000);                       // let the USB enumerate

  Serial.println();
  Serial.println("# waiting 600 s for thermal settle - press any key to skip");
  uint32_t t0 = millis();
  while (millis() - t0 < 600000UL) {
    if (Serial.available()) { while (Serial.available()) Serial.read(); break; }
    if ((millis() - t0) % 30000UL < 5) { Serial.print("# warm ");
                                         Serial.print((millis()-t0)/1000);
                                         Serial.println(" s"); delay(10); }
  }

  uint32_t ts, te;
  acquire(&ts, &te);

  // ---- derive the achieved rate from the acquisition window --------------
  double span_s = (double)(te - ts) / 1e6;
  double fs_hz  = (double)(CHAR_N - 1) / span_s;

  // ---- summary, so the bench notebook can be filled in immediately -------
  double sum = 0.0, sumsq = 0.0;
  uint16_t mn = 65535, mx = 0;
  for (uint32_t i = 0; i < CHAR_N; i++) {
    sum += buf[i]; sumsq += (double)buf[i] * buf[i];
    if (buf[i] < mn) mn = buf[i];
    if (buf[i] > mx) mx = buf[i];
  }
  double mean = sum / CHAR_N;
  double var  = sumsq / CHAR_N - mean * mean;
  double sd   = (var > 0.0) ? sqrt(var) : 0.0;

  Serial.println("# ----------------------------------------------------------");
  Serial.print  ("# label,");        Serial.println(RUN_LABEL);
#if MODE_ROUNDROBIN
  Serial.println("# mode,roundrobin");
#else
  Serial.println("# mode,single");
#endif
  Serial.print  ("# n,");            Serial.println(CHAR_N);
  Serial.print  ("# fs_hz,");        Serial.println(fs_hz, 4);
  Serial.print  ("# span_s,");       Serial.println(span_s, 6);
  Serial.print  ("# mean_code,");    Serial.println(mean, 4);
  Serial.print  ("# sd_code,");      Serial.println(sd, 4);
  Serial.print  ("# min_code,");     Serial.println(mn);
  Serial.print  ("# max_code,");     Serial.println(mx);
  Serial.print  ("# mean_volts,");   Serial.println(mean * 3.3 / 4095.0, 6);
  Serial.print  ("# sd_millivolts,");Serial.println(sd * 3.3 / 4095.0 * 1000.0, 4);
  Serial.println("# adc,12bit_default_11dB_attenuation");
  Serial.println("# conv,linear_3.3_over_4095");
#if MODE_ROUNDROBIN
  // Per-channel summary. Any channel whose mean sits near 0 or near 4095,
  // or whose sd is far larger than Vb's, is unterminated — worth recording
  // because a floating neighbour is a plausible charge-injection source.
  static const char *chName[5] = {"Va", "Vb", "Vc", "Ib", "Ic"};
  for (int c = 0; c < 5; c++) {
    double m = chSum[c] / CHAR_N;
    double v = chSumSq[c] / CHAR_N - m * m;
    Serial.print("# ch_"); Serial.print(chName[c]);
    Serial.print(",mean_code,");  Serial.print(m, 2);
    Serial.print(",sd_code,");    Serial.print((v > 0.0) ? sqrt(v) : 0.0, 4);
    Serial.print(",mean_volts,"); Serial.println(m * 3.3 / 4095.0, 4);
  }
#endif
  Serial.println("# ----------------------------------------------------------");

  for (uint32_t i = 0; i < CHAR_N; i++) Serial.println(buf[i]);

  Serial.println("# end");
}

void loop() { }
