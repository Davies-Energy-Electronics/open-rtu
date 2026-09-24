/* ============================================================================
 * characterise_sweep.ino  —  long capture for the L45 stepped-DC sweep
 * ----------------------------------------------------------------------------
 * SEPARATE BUILD. Do NOT merge any of this into Feather1codeMimoDemoV2.ino or
 * Feather1codeMimoDemoV6_PARALLEL.ino.
 *
 * This is characterise.ino with ONE change: CHAR_N is doubled to 32768, so the
 * record spans 16.384 s at the same 2 kHz per-channel rate. The staircase in
 * feather2_dc_sweep.ino has a period of exactly 16.000 s, so a 16.384 s record
 * is guaranteed to contain every one of the sixteen levels regardless of when
 * the capture starts relative to the staircase. Nothing else differs: same
 * pins, same 12-bit default 11 dB attenuation, same single-channel Vb read,
 * same silent acquisition.
 *
 * RAM: buf is 65,536 bytes. That links comfortably on the ESP32-PICO-MINI-02,
 * but if the build ever fails with a .bss overflow, set CHAR_N back to 16384
 * and SAMPLE_PERIOD_US to 1000. The record still spans 16.384 s; only the rate
 * changes, and a white-noise floor does not care. Record which you used.
 *
 * Output on the serial monitor at 921600 baud:
 *     '#' header lines, then one decimal ADC code per line.
 * Save the whole lot to a .txt and feed it to analyse_sweep.py.
 *
 * Jack Davies, September 2026
 * ========================================================================== */

// ---- pins, identical to the measurement firmware -------------------------
#define PIN_VA   34
#define PIN_VB   39
#define PIN_VC   36
#define PIN_IB   32
#define PIN_IC   33

// ---- acquisition --------------------------------------------------------
#define CHAR_N            32768      // 16.384 s at 2 kHz; covers one full staircase
#define SAMPLE_PERIOD_US  500        // 2 kHz per channel, as the canonical build

// ======================= SET THIS BEFORE EVERY FLASH =======================
#define RUN_LABEL  "C1_short_run3"
// ===========================================================================

static uint16_t buf[CHAR_N];

// ---------------------------------------------------------------------------
void acquire(uint32_t *t_start, uint32_t *t_end) {
  uint32_t next;

  // Discard the first few conversions: the first analogRead after a mode
  // change is not representative.
  for (int i = 0; i < 32; i++) (void)analogRead(PIN_VB);

  next = micros();
  *t_start = next;

  for (uint32_t i = 0; i < CHAR_N; i++) {
    buf[i] = analogRead(PIN_VB);
    next += SAMPLE_PERIOD_US;
    while ((int32_t)(micros() - next) < 0) { /* spin */ }
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

  double span_s = (double)(te - ts) / 1e6;
  double fs_hz  = (double)(CHAR_N - 1) / span_s;

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
  Serial.println("# mode,sweep_single");
  Serial.print  ("# n,");            Serial.println(CHAR_N);
  Serial.print  ("# fs_hz,");        Serial.println(fs_hz, 4);
  Serial.print  ("# span_s,");       Serial.println(span_s, 6);
  Serial.print  ("# mean_code,");    Serial.println(mean, 4);
  Serial.print  ("# sd_code,");      Serial.println(sd, 4);
  Serial.print  ("# min_code,");     Serial.println(mn);
  Serial.print  ("# max_code,");     Serial.println(mx);
  Serial.println("# NOTE sd_code above is the WHOLE-RECORD spread and is");
  Serial.println("#      dominated by the staircase, not by noise. The");
  Serial.println("#      measurement is the spread WITHIN each step, which");
  Serial.println("#      analyse_sweep.py computes.");
  Serial.println("# adc,12bit_default_11dB_attenuation");
  Serial.println("# conv,linear_3.3_over_4095");
  Serial.println("# ----------------------------------------------------------");

  for (uint32_t i = 0; i < CHAR_N; i++) Serial.println(buf[i]);

  Serial.println("# end");
}

void loop() { }
