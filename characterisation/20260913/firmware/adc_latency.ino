/* ============================================================================
 * adc_latency.ino  —  how long analogRead() takes, and how much that varies
 * ----------------------------------------------------------------------------
 * SEPARATE BUILD. Flash to FEATHER 1 (COM6).
 *
 * Why this matters. Both measurement builds timestamp the sample BEFORE the
 * conversion:
 *
 *     uint32_t t = micros();                          // V2 measureFreq()
 *     float v = analogRead(pin) * (3.3f/4095.0f) - bias;
 *
 * so the recorded instant precedes the actual sample-and-hold by the
 * conversion latency. If that latency were constant it would be a fixed
 * offset and would cancel in a period. It is not constant, and the varying
 * part lands directly on every interpolated zero crossing as timing noise.
 *
 * A1 measured the converter's input-referred noise at 0.687 mV, which
 * accounts for only 1.75 mHz of the 5.46 mHz realised floor. Simulation says
 * about 6 us of per-sample timestamp jitter closes the remaining gap. This
 * sketch measures that jitter directly, so the attribution stops being a
 * hypothesis.
 *
 * Two modes, matching the two builds:
 *   MODE_ROUNDROBIN 0 - repeated reads of Vb alone (V2's freq loop)
 *   MODE_ROUNDROBIN 1 - five channels in sequence (V6's samplingTask)
 *
 * Output: one latency in microseconds per line. Capture with
 *   python log_capture.py COM6 D1_latency_single.txt --skip-warmup --expect 20000
 *
 * Jack Davies, September 2026
 * ========================================================================== */

#define PIN_VA   34
#define PIN_VB   39
#define PIN_VC   36
#define PIN_IB   32
#define PIN_IC   33
static const int adcPins[5] = {PIN_VA, PIN_VB, PIN_VC, PIN_IB, PIN_IC};
#define CH_VB    1

#define N_LAT            20000
#define MODE_ROUNDROBIN  0
#define RUN_LABEL        "D1_latency_single"

static uint16_t lat[N_LAT];      // microseconds; 16 bits is ample

void setup() {
  Serial.begin(921600);
  delay(3000);

  // Warm the ADC: the first conversions after boot are not representative.
  for (int i = 0; i < 200; i++) (void)analogRead(PIN_VB);

  uint32_t t0 = micros();
  for (uint32_t i = 0; i < N_LAT; i++) {
#if MODE_ROUNDROBIN
    // Time the Vb conversion in its real position in the round-robin, so
    // any multiplexer settling cost is included.
    for (int ch = 0; ch < 5; ch++) {
      if (ch == CH_VB) {
        uint32_t a = micros();
        (void)analogRead(adcPins[ch]);
        uint32_t b = micros();
        lat[i] = (uint16_t)(b - a);
      } else {
        (void)analogRead(adcPins[ch]);
      }
    }
#else
    uint32_t a = micros();
    (void)analogRead(PIN_VB);
    uint32_t b = micros();
    lat[i] = (uint16_t)(b - a);
#endif
  }
  uint32_t t1 = micros();

  double sum = 0.0, sumsq = 0.0;
  uint16_t mn = 65535, mx = 0;
  for (uint32_t i = 0; i < N_LAT; i++) {
    sum += lat[i]; sumsq += (double)lat[i] * lat[i];
    if (lat[i] < mn) mn = lat[i];
    if (lat[i] > mx) mx = lat[i];
  }
  double mean = sum / N_LAT;
  double var  = sumsq / N_LAT - mean * mean;

  Serial.println("# ----------------------------------------------------------");
  Serial.print  ("# label,");        Serial.println(RUN_LABEL);
#if MODE_ROUNDROBIN
  Serial.println("# mode,latency_roundrobin");
#else
  Serial.println("# mode,latency_single");
#endif
  Serial.print  ("# n,");            Serial.println(N_LAT);
  Serial.print  ("# total_s,");      Serial.println((t1 - t0) / 1e6, 6);
  Serial.print  ("# mean_us,");      Serial.println(mean, 4);
  Serial.print  ("# sd_us,");        Serial.println((var > 0.0) ? sqrt(var) : 0.0, 4);
  Serial.print  ("# min_us,");       Serial.println(mn);
  Serial.print  ("# max_us,");       Serial.println(mx);
  Serial.println("# note,timestamp_taken_before_conversion_in_both_builds");
  Serial.println("# ----------------------------------------------------------");

  for (uint32_t i = 0; i < N_LAT; i++) Serial.println(lat[i]);
  Serial.println("# end");
}

void loop() { }
