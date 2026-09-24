/* ============================================================================
 * feather2_dc_sweep.ino  —  stepped-DC staircase source for L45
 * ----------------------------------------------------------------------------
 * SEPARATE BUILD. Do not merge into feather_2_code_mimoDEMO.ino or
 * feather2_dc_hold.ino. Those are the provenance of archived captures and
 * must stay byte-identical to the released sources.
 *
 * WHY THIS SKETCH EXISTS
 * The quiet-input noise floor (A1, 13 September) was measured at ONE DC level:
 * the midpoint, DAC code 127. That gives one number for sigma_v and the error
 * budget of Section 3.12 assumes it holds everywhere. If the converter's noise
 * is code-dependent — differential non-linearity, reference coupling, or a
 * settling term that varies with the step the multiplexer has to make — then
 * the single midpoint figure is not the figure that applies during a replay,
 * where the signal traverses the whole range every cycle.
 *
 * This sketch drives DAC_VB through sixteen evenly spaced levels, one second
 * each, repeating forever. Feather 1 captures straight through it with
 * characterise_sweep.ino and the levels are recovered offline by their own
 * step changes, so NO trigger wire and NO cross-board synchronisation is
 * needed. The steps are hundreds of ADC codes apart and the noise is under one
 * code, so segmenting the capture is unambiguous.
 *
 * DAC_VA is held at the midpoint throughout, exactly as feather2_dc_hold does,
 * so nothing about the other channel's source impedance changes.
 *
 * DWELL: 25 ms, NOT one second. Stage two of the conditioning chain is an
 * alternating-current coupling whose time constant is 533.7 ms, measured on
 * this bench. A level held for one second decays 85 % of the way back to the
 * bias point before the next arrives, so a one-second staircase measures the
 * high-pass transient and not the converter. At 25 ms the level falls by only
 * 4.6 % across a dwell, which the analyser removes with a linear detrend, and
 * the staircase reaches the converter essentially intact. Validated in
 * simulation against the measured time constant before any bench time was
 * spent on it.
 *
 * A 25 ms dwell also means the signal is FLAT inside each dwell, which is what
 * makes this work where measuring on a sine did not: there is no waveform
 * curvature for a local fit to mis-track and no generator harmonic content to
 * leak into the residual and imitate code-dependent noise.
 *
 * LEVELS: 8, 24, 40, ... 248. The two extremes of the DAC are deliberately
 * avoided: the ESP32 DAC is non-linear in the last few codes at each rail and
 * a noise figure taken there would be measuring the DAC, not the ADC chain.
 *
 * Flash to FEATHER 2 (the generator board). Leave Feather 1 alone.
 * Reflash feather_2_code_mimoDEMO.ino afterwards.
 *
 * Jack Davies, September 2026
 * ========================================================================== */

#define DAC_VA        25
#define DAC_VB        26
#define DC_MID        127       // same midpoint constant the replay LUT uses

#define N_LEVELS      16
#define LEVEL_FIRST   8         // lowest DAC code in the staircase
#define LEVEL_STEP    16        // spacing, so the top level is 8 + 15*16 = 248
#define DWELL_MS      25UL      // 25 ms per level - see note below

void setup() {
  Serial.begin(115200);
  delay(1000);

  dacWrite(DAC_VA, DC_MID);     // the untouched channel stays at midpoint
  dacWrite(DAC_VB, LEVEL_FIRST);

  Serial.println("# feather2_dc_sweep");
  Serial.print  ("# n_levels,");   Serial.println(N_LEVELS);
  Serial.print  ("# level_first,");Serial.println(LEVEL_FIRST);
  Serial.print  ("# level_step,"); Serial.println(LEVEL_STEP);
  Serial.print  ("# dwell_ms,");   Serial.println(DWELL_MS);
  Serial.print  ("# period_s,");   Serial.println(N_LEVELS * DWELL_MS / 1000.0, 3);
  Serial.println("# DAC_VB steps; DAC_VA held at 127. Repeats until reset.");
  Serial.println("# Start the Feather 1 capture at any time - the levels are");
  Serial.println("# recovered offline from the step changes themselves.");
}

void loop() {
  for (int i = 0; i < N_LEVELS; i++) {
    uint8_t code = LEVEL_FIRST + i * LEVEL_STEP;
    dacWrite(DAC_VB, code);
    dacWrite(DAC_VA, DC_MID);   // belt and braces, as in feather2_dc_hold
    Serial.print("# level,"); Serial.print(i);
    Serial.print(",dac,");    Serial.println(code);
    delay(DWELL_MS);
  }
}
