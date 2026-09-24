/* ============================================================================
 * characterise_v2loop.ino  -  isolate the V2 acquisition loop's contribution
 *                             to the frequency noise floor
 * ----------------------------------------------------------------------------
 * SEPARATE BUILD. Do NOT merge any of this into Feather1codeMimoDemoV2.ino or
 * Feather1codeMimoDemoV6_PARALLEL.ino. Those two are the provenance of every
 * archived capture and must stay byte-identical to the released sources.
 *
 * WHY THIS SKETCH EXISTS
 * The 13 September characterisation compared a deadline-scheduled,
 * single-channel capture against a June capture taken with V2's own
 * acquisition loop. TWO things differed at once, so the 4.50 mHz attributed to
 * the acquisition loop is an upper bound, not an isolation. This sketch walks
 * from one to the other ONE CHANGE AT A TIME, so each element of V2's loop
 * gets its own number.
 *
 * WHAT ACTUALLY DIFFERS, read from Feather1codeMimoDemoV2.ino:
 *   1. V2 timestamps BEFORE the conversion:
 *          uint32_t t = micros();
 *          float v = analogRead(pin) * (3.3f/4095.0f) - bias;
 *      characterise.ino does not timestamp at all; it spins to a deadline.
 *   2. V2 paces with delayMicroseconds(500) at the END of the iteration, so
 *      the period is FREE-RUNNING: 500 us plus however long micros(), the
 *      conversion, the scaling and the crossing arithmetic took.
 *      characterise.ino spins to a fixed 500 us deadline, which absorbs that
 *      variation instead of passing it on.
 *   3. V2 does float work between samples whose duration DEPENDS ON THE DATA:
 *      a sample that completes a zero crossing costs an interpolation that a
 *      sample in mid-half-cycle does not. That is jitter correlated with
 *      exactly the samples the estimator relies on.
 *
 * NOTE what does NOT differ: the RMS, THD, phase and power measurements do
 * NOT interleave with the frequency samples. void loop() calls measureRMS
 * five times, THEN measureFreq three times, THEN measureTHD, and so on. The
 * 200 samples inside one measureFreq call are back-to-back. So "interleaved
 * acquisition" cannot act inside a frequency window, and this sketch does not
 * test it. Between-window effects are covered by running the real V2 build.
 *
 * MODES - set MODE to 0, 1, 2 or 3 before each flash.
 *   0  SEPT_BASELINE  deadline spin-wait, timestamp after conversion.
 *                     Reproduces the 13 September A1-A3 condition.
 *   1  TS_BEFORE      deadline spin-wait, timestamp BEFORE conversion.
 *                     Isolates the timestamp-ordering term.
 *   2  FREERUN        delayMicroseconds(500), timestamp before.
 *                     Isolates V2's free-running pacing.
 *   3  V2_ARITHMETIC  as mode 2, plus V2's bias subtraction, crossing tests
 *                     and sub-sample interpolation inside the loop.
 *                     Isolates the data-dependent jitter.
 *
 * Every mode logs one line per sample: an absolute microsecond timestamp and
 * the raw converter code. Feed the file to characterise_loop.py, which reports
 * the achieved sampling rate from the mean interval (L41) and the jitter from
 * its spread (L40).
 *
 * Acquisition is silent. Nothing is transmitted until the buffer is full,
 * because USB traffic during acquisition perturbs the shared supply rail and
 * you would be measuring the serial bridge.
 *
 * Output on the serial monitor at 921600 baud:
 *     '#' header lines, then "t_us,code", then one CSV row per sample.
 *
 * Jack Davies, September 2026
 * ========================================================================== */

// ---- pins, identical to the measurement firmware -------------------------
#define PIN_VA   34
#define PIN_VB   39
#define PIN_VC   36
#define PIN_IB   32
#define PIN_IC   33

// ---- acquisition ---------------------------------------------------------
#define CHAR_N            16384      // 8.192 s at 2 kHz, as September
#define SAMPLE_PERIOD_US  500        // the nominal per-channel period
#define WARMUP_MS         600000UL   // 600 s thermal settle, as September

// ======================= SET THESE TWO BEFORE EVERY FLASH ==================
#define MODE        3
#define RUN_LABEL   "G_stim_v2arith"
// ===========================================================================

#define MODE_SEPT_BASELINE   0
#define MODE_TS_BEFORE       1
#define MODE_FREERUN         2
#define MODE_V2_ARITHMETIC   3

static uint16_t codeBuf[CHAR_N];
static uint16_t dtBuf[CHAR_N];       // interval since the previous mark, us

// Bias, measured the way V2 measures it, used only in mode 3.
static float biasVb = 0.0f;

// Crossing bookkeeping for mode 3. The results are discarded: this sketch
// measures the COST of doing the arithmetic, not the frequency. Keeping the
// outputs would invite them to be quoted, and they are not a measurement.
static volatile float  m3_prev_v = 0.0f;
static volatile uint32_t m3_prev_t = 0;
static volatile bool   m3_first = true;
static volatile uint32_t m3_sink = 0;     // stops the optimiser removing the work
static uint32_t m3_nPos = 0, m3_nNeg = 0;

// ---------------------------------------------------------------------------
static const char *modeName() {
  switch (MODE) {
    case MODE_SEPT_BASELINE: return "0_sept_baseline_deadline_ts_after";
    case MODE_TS_BEFORE:     return "1_deadline_ts_before";
    case MODE_FREERUN:       return "2_freerun_ts_before";
    case MODE_V2_ARITHMETIC: return "3_freerun_ts_before_v2_arithmetic";
    default:                 return "unknown";
  }
}

// ---------------------------------------------------------------------------
// V2's bias calibration, reproduced: mean of 200 reads at the DC midpoint.
static float measureBias(int pin) {
  double s = 0.0;
  for (int i = 0; i < 200; i++) {
    s += analogRead(pin) * (3.3 / 4095.0);
    delayMicroseconds(500);
  }
  return (float)(s / 200.0);
}

// ---------------------------------------------------------------------------
// The V2 crossing arithmetic, lifted from measureFreq() so its per-sample
// cost is reproduced exactly. Return value is deliberately thrown away.
static inline void v2CrossingWork(float v, uint32_t t) {
  if (!m3_first) {
    if (m3_prev_v < 0.0f && v >= 0.0f) {
      float dv = v - m3_prev_v;
      float frac = (dv > 1e-6f) ? (-m3_prev_v / dv) : 0.5f;
      if (frac < 0.0f) frac = 0.0f;
      if (frac > 1.0f) frac = 1.0f;
      uint32_t t_cross = m3_prev_t + (uint32_t)(frac * (float)(t - m3_prev_t));
      m3_sink += t_cross;
      m3_nPos++;
    } else if (m3_prev_v >= 0.0f && v < 0.0f) {
      float dv = m3_prev_v - v;
      float frac = (dv > 1e-6f) ? (m3_prev_v / dv) : 0.5f;
      if (frac < 0.0f) frac = 0.0f;
      if (frac > 1.0f) frac = 1.0f;
      uint32_t t_cross = m3_prev_t + (uint32_t)(frac * (float)(t - m3_prev_t));
      m3_sink += t_cross;
      m3_nNeg++;
    }
  } else {
    m3_first = false;
  }
  m3_prev_v = v;
  m3_prev_t = t;
}

// ---------------------------------------------------------------------------
void acquire(uint32_t *t_start, uint32_t *t_end) {
  // Discard the first few conversions: the first analogRead after a mode
  // change is not representative.
  for (int i = 0; i < 32; i++) (void)analogRead(PIN_VB);

  uint32_t next = micros();
  uint32_t prev = next;
  *t_start = next;

  for (uint32_t i = 0; i < CHAR_N; i++) {

#if (MODE == MODE_SEPT_BASELINE)
    // Timestamp AFTER the conversion.
    uint16_t code = (uint16_t)analogRead(PIN_VB);
    uint32_t mark = micros();
#else
    // Timestamp BEFORE the conversion, exactly as V2 measureFreq() does.
    uint32_t mark = micros();
    uint16_t code = (uint16_t)analogRead(PIN_VB);
#endif

    codeBuf[i] = code;
    dtBuf[i]   = (uint16_t)(mark - prev);
    prev       = mark;

#if (MODE == MODE_V2_ARITHMETIC)
    // V2 scales and de-biases every sample, then runs the crossing tests.
    float v = code * (3.3f / 4095.0f) - biasVb;
    v2CrossingWork(v, mark);
#endif

#if (MODE == MODE_FREERUN) || (MODE == MODE_V2_ARITHMETIC)
    // Free-running, as V2: the period is 500 us PLUS whatever the work above
    // cost. This is the pacing whose jitter we are trying to measure.
    delayMicroseconds(SAMPLE_PERIOD_US);
#else
    // Deadline-scheduled, as characterise.ino: absorbs the work's variation.
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
  Serial.print  ("# sketch,characterise_v2loop");
  Serial.println();
  Serial.print  ("# mode_number,");  Serial.println(MODE);
  Serial.print  ("# mode_name,");    Serial.println(modeName());
  Serial.print  ("# label,");        Serial.println(RUN_LABEL);
  Serial.println("# ---- CHECK THE THREE LINES ABOVE BEFORE YOU KEEP THIS CAPTURE ----");

  Serial.println("# waiting 600 s for thermal settle - press any key to skip");
  uint32_t t0 = millis();
  while (millis() - t0 < WARMUP_MS) {
    if (Serial.available()) { while (Serial.available()) Serial.read(); break; }
    if ((millis() - t0) % 30000UL < 5) {
      Serial.print("# warm "); Serial.print((millis() - t0) / 1000);
      Serial.println(" s"); delay(10);
    }
  }

#if (MODE == MODE_V2_ARITHMETIC)
  Serial.println("# measuring bias (Feather 2 DC midpoint), as V2 does");
  biasVb = measureBias(PIN_VB);
  Serial.print("# bias_volts,"); Serial.println(biasVb, 6);
#endif

  uint32_t ts, te;
  acquire(&ts, &te);

  // ---- achieved rate over the whole window -------------------------------
  double span_s = (double)(te - ts) / 1e6;
  double fs_hz  = (double)(CHAR_N - 1) / span_s;

  // ---- interval statistics, computed on the board so the bench notebook
  //      can be filled in before the file is even saved --------------------
  double dsum = 0.0, dsumsq = 0.0;
  uint16_t dmn = 65535, dmx = 0;
  for (uint32_t i = 1; i < CHAR_N; i++) {       // skip [0]: no predecessor
    double d = dtBuf[i];
    dsum += d; dsumsq += d * d;
    if (dtBuf[i] < dmn) dmn = dtBuf[i];
    if (dtBuf[i] > dmx) dmx = dtBuf[i];
  }
  double dn    = (double)(CHAR_N - 1);
  double dmean = dsum / dn;
  double dvar  = dsumsq / dn - dmean * dmean;
  double dsd   = (dvar > 0.0) ? sqrt(dvar) : 0.0;

  // ---- code statistics ---------------------------------------------------
  double sum = 0.0, sumsq = 0.0;
  uint16_t mn = 65535, mx = 0;
  for (uint32_t i = 0; i < CHAR_N; i++) {
    sum += codeBuf[i]; sumsq += (double)codeBuf[i] * codeBuf[i];
    if (codeBuf[i] < mn) mn = codeBuf[i];
    if (codeBuf[i] > mx) mx = codeBuf[i];
  }
  double mean = sum / CHAR_N;
  double var  = sumsq / CHAR_N - mean * mean;
  double sd   = (var > 0.0) ? sqrt(var) : 0.0;

  Serial.println("# ----------------------------------------------------------");
  Serial.print  ("# n,");              Serial.println(CHAR_N);
  Serial.print  ("# fs_hz_window,");   Serial.println(fs_hz, 4);
  Serial.print  ("# span_s,");         Serial.println(span_s, 6);
  Serial.print  ("# dt_mean_us,");     Serial.println(dmean, 4);
  Serial.print  ("# dt_sd_us,");       Serial.println(dsd, 4);
  Serial.print  ("# dt_min_us,");      Serial.println(dmn);
  Serial.print  ("# dt_max_us,");      Serial.println(dmx);
  Serial.print  ("# fs_hz_from_mean,");Serial.println(1.0e6 / dmean, 4);
  Serial.print  ("# mean_code,");      Serial.println(mean, 4);
  Serial.print  ("# sd_code,");        Serial.println(sd, 4);
  Serial.print  ("# min_code,");       Serial.println(mn);
  Serial.print  ("# max_code,");       Serial.println(mx);
  Serial.print  ("# mean_volts,");     Serial.println(mean * 3.3 / 4095.0, 6);
  Serial.print  ("# sd_millivolts,");  Serial.println(sd * 3.3 / 4095.0 * 1000.0, 4);
  Serial.println("# adc,12bit_default_11dB_attenuation");
  Serial.println("# conv,linear_3.3_over_4095");
#if (MODE == MODE_V2_ARITHMETIC)
  Serial.print  ("# crossings_pos,"); Serial.println(m3_nPos);
  Serial.print  ("# crossings_neg,"); Serial.println(m3_nNeg);
  Serial.print  ("# sink,");          Serial.println(m3_sink);   // ignore
#endif
  Serial.println("# ----------------------------------------------------------");

  // ---- the data ----------------------------------------------------------
  Serial.println("t_us,code");
  uint32_t t = 0;
  for (uint32_t i = 0; i < CHAR_N; i++) {
    t += dtBuf[i];
    Serial.print(t); Serial.print(','); Serial.println(codeBuf[i]);
  }

  Serial.println("# end");
}

void loop() { }
