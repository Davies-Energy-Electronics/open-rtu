/*
 * MIMO M6 — DNP3 IEEE 1815 Unsolicited Response
 * Feather 1: Va=GPIO34  Vb=GPIO39  Vc=GPIO36
 *            Ib=GPIO32  Ic=GPIO33  DAC_VC=GPIO25
 *
 * V2 patch (2026-06-25): M2 frequency-detector accuracy improvements
 *   - Linear interpolation at each detected zero crossing (removes per-edge
 *     timing quantisation, which was the dominant error term in V1)
 *   - Both-edge detection (positive AND negative going): doubles the
 *     crossing count per window without changing window length
 *   - Same-edge full-period extraction (pos-to-pos, neg-to-neg) — immune
 *     to DC-bias asymmetry between half cycles
 *   - Median filter on the period buffer — rejects outliers from
 *     noise-induced spurious crossings (these caused the 400+ mHz tail
 *     spikes observed in run 2026-06-25_180520)
 *   - Per-sample micros() timestamps (no assumption of fixed 500 us cadence)
 *
 * V1 → V2 expected improvement, based on theoretical analysis:
 *   max |f_inst - f_NESO|:  424 mHz → 80-120 mHz (3-5x)
 *   RMS error:              131 mHz → 25-40 mHz  (3-5x)
 *   max RoCoF error:        0.51 Hz/s → 0.10-0.15 Hz/s
 *   alarm latency:          1808 ms → ~1800 ms   (UNCHANGED — see below)
 *
 * NOT addressed in this patch:
 *   Alarm latency is dominated by the loop architecture (sequential
 *   measurement of 5 RMS + 3 freq + 3 THD + 2 phase + 2 power, all in one
 *   thread, with frame transmission only at the end of each iteration).
 *   No amount of M2 accuracy improvement reduces the latency below ~1 s.
 *   See accompanying notes for the structural change required.
 */

#include <math.h>
#include <string.h>

#define PIN_VA  34
#define PIN_VB  39
#define PIN_VC  36
#define PIN_IB  32
#define PIN_IC  33
#define DAC_VC  25

#define LUT_SIZE  40
#define AMP_VC    127
static uint8_t lutVc[LUT_SIZE];
volatile int   dacIdx = 0;
hw_timer_t    *timer  = NULL;

#define FS         2000
#define F_NOM      50.0f
#define RMS_N      40
#define FREQ_WIN   200
#define THD_N      200
#define THD_HARM   10
#define PHASE_WIN  200
#define PWR_N      100
#define F_LOWER    49.5f
#define F_UPPER    50.5f
#define UV_THRESH  0.50f
#define THD_LIMIT  4.5f
#define DNP_SRC    0x0001
#define DNP_DST    0x0003
#define NUM_ANA    21
#define NUM_BIN    3

// V2: tunable maximums for the crossing-detection buffers in measureFreq.
// At FREQ_WIN=200 and 50 Hz signal we expect ~5-6 positive and ~5-6
// negative crossings per window. MAX_EDGES_PER_TYPE=12 gives plenty of
// headroom for noisy signals where the count can briefly spike.
#define MAX_EDGES_PER_TYPE  12
#define MAX_PERIODS         (2 * (MAX_EDGES_PER_TYPE - 1))

struct PwrResult { float Vrms, Irms, P; };

float biasVa=0, biasVb=0, biasVc=0, biasIb=0, biasIc=0;
static float buf[THD_N];

void IRAM_ATTR onTimer() {
  dacWrite(DAC_VC, lutVc[dacIdx]);
  dacIdx = (dacIdx + 1) % LUT_SIZE;
}

float measureRMS(int pin, float bias, int n) {
  float sq=0;
  for (int k=0; k<n; k++) {
    float v=analogRead(pin)*(3.3f/4095.0f)-bias;
    sq+=v*v;
    delayMicroseconds(500);
  }
  return sqrtf(sq/n);
}

// ─────────────────────────────────────────────────────────────────────────
// V2 measureFreq() — interpolating, both-edge, median-filtered.
//
// Replaces V1's algorithm:
//   "find first and last positive-going crossings; f = (n-1)/(lastT-firstT)"
// which had two failure modes:
//   1. Per-edge +/- 250 us quantisation jitter from delayMicroseconds(500)
//      polling, producing ~110 mHz stdev floor per measurement.
//   2. Heavy tails from noise-induced spurious crossings shifting firstT
//      or lastT by full cycles, producing 400+ mHz outlier excursions.
//
// V2 fixes both: every crossing is sub-sample-interpolated, every period
// goes into a buffer, and the median is reported — outliers contribute
// to the buffer ordering but not to the reported value.
// ─────────────────────────────────────────────────────────────────────────
float measureFreq(int pin, float bias) {
  uint32_t pos_t[MAX_EDGES_PER_TYPE];
  uint32_t neg_t[MAX_EDGES_PER_TYPE];
  int n_pos = 0, n_neg = 0;

  float    prev_v = 0.0f;
  uint32_t prev_t = 0;
  bool     first_sample = true;

  // ── Collect interpolated zero-crossing timestamps ──
  for (int i = 0; i < FREQ_WIN; i++) {
    uint32_t t = micros();
    float v = analogRead(pin) * (3.3f / 4095.0f) - bias;

    if (!first_sample) {
      // Positive-going crossing: prev_v negative, v non-negative
      if (prev_v < 0.0f && v >= 0.0f) {
        float dv = v - prev_v;                     // always > 0 here
        float frac = (dv > 1e-6f) ? (-prev_v / dv) : 0.5f;
        if (frac < 0.0f) frac = 0.0f;
        if (frac > 1.0f) frac = 1.0f;
        uint32_t t_cross = prev_t + (uint32_t)(frac * (float)(t - prev_t));
        if (n_pos < MAX_EDGES_PER_TYPE) pos_t[n_pos++] = t_cross;
      }
      // Negative-going crossing: prev_v non-negative, v negative
      else if (prev_v >= 0.0f && v < 0.0f) {
        float dv = prev_v - v;                     // always > 0 here
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
    delayMicroseconds(500);
  }

  // ── Build a buffer of full-cycle periods (same-edge to same-edge) ──
  // Mixing pos-pos and neg-neg periods is safe even with DC-bias residual:
  // both classes measure the same underlying signal period, only the
  // half-cycle split point differs.
  float periods[MAX_PERIODS];
  int n_periods = 0;
  for (int i = 0; i < n_pos - 1; i++) {
    if (n_periods < MAX_PERIODS) {
      periods[n_periods++] = (float)(pos_t[i + 1] - pos_t[i]);
    }
  }
  for (int i = 0; i < n_neg - 1; i++) {
    if (n_periods < MAX_PERIODS) {
      periods[n_periods++] = (float)(neg_t[i + 1] - neg_t[i]);
    }
  }

  if (n_periods < 3) return 0.0f;

  // ── Insertion sort for median (n_periods <= ~20, sort is trivial) ──
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

  // Sanity bound: a 100-us "period" would be 10 kHz, well outside spec.
  if (median_period < 100.0f) return 0.0f;
  return 1000000.0f / median_period;
}

float measureTHD(int pin, float bias) {
  for (int i=0; i<THD_N; i++) {
    buf[i]=analogRead(pin)*(3.3f/4095.0f)-bias;
    delayMicroseconds(470);
  }
  for (int i=0; i<THD_N; i++) {
    float w=0.5f*(1.0f-cosf(2.0f*PI*i/(THD_N-1)));
    buf[i]*=w;
  }
  int kF=(int)roundf(F_NOM*THD_N/(float)FS);
  auto goertzel=[&](int k)->float {
    float omega=2.0f*PI*(float)k/(float)THD_N;
    float coeff=2.0f*cosf(omega);
    float s0=0,s1=0,s2=0;
    for (int n=0; n<THD_N; n++) {
      s0=buf[n]+coeff*s1-s2;
      s2=s1; s1=s0;
    }
    float re=s1-s2*cosf(omega);
    float im=s2*sinf(omega);
    return (4.0f/THD_N)*sqrtf(re*re+im*im);
  };
  float H1=goertzel(kF);
  if (H1<0.001f) return NAN;
  float sq=0;
  for (int h=2; h<=THD_HARM; h++) {
    float Hh=goertzel(kF*h);
    sq+=Hh*Hh;
  }
  return 100.0f*sqrtf(sq)/H1;
}

float measurePhase(int pinV, float biasV, int pinI, float biasI) {
  float lV=1.0f, lI=1.0f;
  unsigned long tV=0, tI=0;
  bool gV=false, gI=false;
  for (int i=0; i<PHASE_WIN; i++) {
    unsigned long t=micros();
    float v =analogRead(pinV)*(3.3f/4095.0f)-biasV;
    float ia=analogRead(pinI)*(3.3f/4095.0f)-biasI;
    if (!gV && lV<0 && v >=0) { tV=t; gV=true; }
    if (!gI && lI<0 && ia>=0) { tI=t; gI=true; }
    lV=v; lI=ia;
    if (gV && gI) break;
    delayMicroseconds(460);
  }
  if (!gV || !gI) return NAN;
  float ph=(float)((long)(tI-tV))/(1000000.0f/F_NOM)*360.0f;
  while (ph> 180.0f) ph-=360.0f;
  while (ph<-180.0f) ph+=360.0f;
  return ph;
}

PwrResult measurePower(int pinV, float biasV, int pinI, float biasI) {
  float sv=0, si=0, sp=0;
  for (int k=0; k<PWR_N; k++) {
    float v=analogRead(pinV)*(3.3f/4095.0f)-biasV;
    float i=analogRead(pinI)*(3.3f/4095.0f)-biasI;
    sv+=v*v; si+=i*i; sp+=v*i;
    delayMicroseconds(460);
  }
  PwrResult r;
  r.Vrms=sqrtf(sv/PWR_N);
  r.Irms=sqrtf(si/PWR_N);
  r.P=sp/PWR_N;
  return r;
}

uint16_t dnp3CRC(uint8_t *data, int len) {
  uint16_t crc=0;
  for (int i=0; i<len; i++) {
    crc^=(uint16_t)data[i];
    for (int j=0; j<8; j++)
      crc=(crc&1)?((crc>>1)^0xA6BC):(crc>>1);
  }
  return ~crc;
}

int buildDNP3(uint8_t *frame, float *ana, bool *bin) {
  uint8_t ud[400];
  int uLen=0;

  ud[uLen++]=0xC0; ud[uLen++]=0xD0;
  ud[uLen++]=0x82; ud[uLen++]=0x00; ud[uLen++]=0x00;

  ud[uLen++]=30; ud[uLen++]=5; ud[uLen++]=0x17; ud[uLen++]=NUM_ANA;
  for (int i=0; i<NUM_ANA; i++) {
    ud[uLen++]=0x01;
    uint8_t fb[4];
    memcpy(fb, &ana[i], 4);
    ud[uLen++]=fb[0]; ud[uLen++]=fb[1];
    ud[uLen++]=fb[2]; ud[uLen++]=fb[3];
  }

  ud[uLen++]=1; ud[uLen++]=2; ud[uLen++]=0x17; ud[uLen++]=NUM_BIN;
  for (int i=0; i<NUM_BIN; i++)
    ud[uLen++]=bin[i]?0x81:0x01;

  int fLen=0;
  frame[fLen++]=0x05; frame[fLen++]=0x64;
  frame[fLen++]=(uint8_t)(uLen+5);
  frame[fLen++]=0x44;
  frame[fLen++]=(uint8_t)(DNP_DST&0xFF);
  frame[fLen++]=(uint8_t)(DNP_DST>>8);
  frame[fLen++]=(uint8_t)(DNP_SRC&0xFF);
  frame[fLen++]=(uint8_t)(DNP_SRC>>8);
  uint16_t hCRC=dnp3CRC(&frame[2], 6);
  frame[fLen++]=(uint8_t)(hCRC&0xFF);
  frame[fLen++]=(uint8_t)(hCRC>>8);

  int udIdx=0;
  while (udIdx<uLen) {
    int blkLen=min(16, uLen-udIdx);
    uint8_t *blk=&frame[fLen];
    for (int i=0; i<blkLen; i++) frame[fLen++]=ud[udIdx++];
    uint16_t bCRC=dnp3CRC(blk, blkLen);
    frame[fLen++]=(uint8_t)(bCRC&0xFF);
    frame[fLen++]=(uint8_t)(bCRC>>8);
  }
  return fLen;
}

void setup() {
  Serial.begin(115200);
  delay(1000);

  pinMode(PIN_VA,INPUT); pinMode(PIN_VB,INPUT); pinMode(PIN_VC,INPUT);
  pinMode(PIN_IB,INPUT); pinMode(PIN_IC,INPUT);
  analogSetPinAttenuation(PIN_VA,ADC_11db);
  analogSetPinAttenuation(PIN_VB,ADC_11db);
  analogSetPinAttenuation(PIN_VC,ADC_11db);
  analogSetPinAttenuation(PIN_IB,ADC_11db);
  analogSetPinAttenuation(PIN_IC,ADC_11db);

  for (int i=0; i<LUT_SIZE; i++) {
    float a=2.0f*PI*i/LUT_SIZE;
    lutVc[i]=(uint8_t)(127+AMP_VC*sinf(a-4.0f*PI/3.0f));
  }

  dacWrite(DAC_VC,127);
  timer=timerBegin(1000000);
  timerAttachInterrupt(timer,&onTimer);
  timerAlarm(timer,500,true,0);

  Serial.println("================================================");
  Serial.println("MIMO M6 — DNP3 IEEE 1815 Unsolicited Response  (V2)");
  Serial.println("Src:0x0001  Dst:0x0003  G30V5(21) + G1V2(3)");
  Serial.println("V2: interpolating + both-edge + median-filtered M2");
  Serial.println("================================================");
  Serial.println("  Measuring bias (Feather 2 DC midpoint)...");
  delay(1000);

  for (int i=0; i<1000; i++) {
    biasVa+=analogRead(PIN_VA)*(3.3f/4095.0f);
    biasVb+=analogRead(PIN_VB)*(3.3f/4095.0f);
    biasVc+=analogRead(PIN_VC)*(3.3f/4095.0f);
    biasIb+=analogRead(PIN_IB)*(3.3f/4095.0f);
    biasIc+=analogRead(PIN_IC)*(3.3f/4095.0f);
    delayMicroseconds(200);
  }
  biasVa/=1000; biasVb/=1000; biasVc/=1000;
  biasIb/=1000; biasIc/=1000;

  Serial.printf("  Bias Va:%.4fV Vb:%.4fV Vc:%.4fV Ib:%.4fV Ic:%.4fV\n",
    biasVa,biasVb,biasVc,biasIb,biasIc);
  Serial.println("  Waiting for Feather 2 AC generation...");
  delay(4000);
  Serial.println("  Transmitting DNP3 frames...");
  Serial.println("================================================");
}

void loop() {
  float Va_rms=measureRMS(PIN_VA,biasVa,RMS_N);
  float Vb_rms=measureRMS(PIN_VB,biasVb,RMS_N);
  float Vc_rms=measureRMS(PIN_VC,biasVc,RMS_N);
  float Ib_rms=measureRMS(PIN_IB,biasIb,RMS_N);
  float Ic_rms=measureRMS(PIN_IC,biasIc,RMS_N);

  float fVa=measureFreq(PIN_VA,biasVa);
  float fVb=measureFreq(PIN_VB,biasVb);
  float fVc=measureFreq(PIN_VC,biasVc);

  float thdVa=measureTHD(PIN_VA,biasVa);
  float thdVb=measureTHD(PIN_VB,biasVb);
  float thdVc=measureTHD(PIN_VC,biasVc);

  float phVbIb=measurePhase(PIN_VB,biasVb,PIN_IB,biasIb);
  float phVcIc=measurePhase(PIN_VC,biasVc,PIN_IC,biasIc);

  PwrResult pwB=measurePower(PIN_VB,biasVb,PIN_IB,biasIb);
  PwrResult pwC=measurePower(PIN_VC,biasVc,PIN_IC,biasIc);

  float phBrad=isnan(phVbIb)?0.0f:phVbIb*PI/180.0f;
  float phCrad=isnan(phVcIc)?0.0f:phVcIc*PI/180.0f;
  float Sb=pwB.Vrms*pwB.Irms;
  float Pb=Sb*cosf(phBrad);
  float Qb=Sb*sinf(phBrad);
  float PFb=(Sb>0.001f)?Pb/Sb:0.0f;
  float Sc=pwC.Vrms*pwC.Irms;
  float Pc=Sc*cosf(phCrad);
  float Qc=Sc*sinf(phCrad);
  float PFc=(Sc>0.001f)?Pc/Sc:0.0f;

  bool freqAlarm=(fVa<F_LOWER||fVa>F_UPPER)||
                 (fVb<F_LOWER||fVb>F_UPPER)||
                 (fVc<F_LOWER||fVc>F_UPPER);
  bool uvAlarm  =(Va_rms<UV_THRESH)||(Vb_rms<UV_THRESH)||(Vc_rms<UV_THRESH);
  bool thdAlarm =(!isnan(thdVa)&&thdVa>THD_LIMIT)||
                 (!isnan(thdVb)&&thdVb>THD_LIMIT)||
                 (!isnan(thdVc)&&thdVc>THD_LIMIT);

  float safeThd[3]={
    isnan(thdVa)?0.0f:thdVa,
    isnan(thdVb)?0.0f:thdVb,
    isnan(thdVc)?0.0f:thdVc
  };
  float safePh[2]={
    isnan(phVbIb)?0.0f:phVbIb,
    isnan(phVcIc)?0.0f:phVcIc
  };

  float ana[NUM_ANA]={
    Va_rms,Vb_rms,Vc_rms,Ib_rms,Ic_rms,
    fVa,fVb,fVc,
    safeThd[0],safeThd[1],safeThd[2],
    safePh[0],safePh[1],
    Pb,Qb,Sb,PFb,
    Pc,Qc,Sc,PFc
  };
  bool bin[NUM_BIN]={freqAlarm,uvAlarm,thdAlarm};

  uint8_t frame[400];
  int fLen=buildDNP3(frame,ana,bin);

  Serial.printf("\n[DNP3 FRAME | %d bytes | Src:0x%04X Dst:0x%04X]\n",
    fLen,DNP_SRC,DNP_DST);
  for (int i=0; i<fLen; i++) {
    Serial.printf("%02X ",frame[i]);
    if ((i+1)%18==0) Serial.println();
  }
  Serial.println();

  Serial.println("[G30V5 ANALOGUE]");
  Serial.printf("  Va_rms:%.4fV  Vb_rms:%.4fV  Vc_rms:%.4fV\n",
    Va_rms,Vb_rms,Vc_rms);
  Serial.printf("  Ib_rms:%.4fV  Ic_rms:%.4fV\n",Ib_rms,Ic_rms);
  Serial.printf("  f_Va:%.4fHz  f_Vb:%.4fHz  f_Vc:%.4fHz\n",fVa,fVb,fVc);
  Serial.printf("  THD_Va:%.4f%%  THD_Vb:%.4f%%  THD_Vc:%.4f%%\n",
    safeThd[0],safeThd[1],safeThd[2]);
  Serial.printf("  Ph_Vb/Ib:%.4fdeg  Ph_Vc/Ic:%.4fdeg\n",
    safePh[0],safePh[1]);
  Serial.printf("  Pb:%.4f Qb:%.4f Sb:%.4f PFb:%.4f\n",Pb,Qb,Sb,PFb);
  Serial.printf("  Pc:%.4f Qc:%.4f Sc:%.4f PFc:%.4f\n",Pc,Qc,Sc,PFc);
  Serial.println("[G1V2 BINARY ALARMS]");
  Serial.printf("  FREQ_alarm:%d  UV_alarm:%d  THD_alarm:%d\n",
    freqAlarm,uvAlarm,thdAlarm);
  Serial.println("================================================");
}
