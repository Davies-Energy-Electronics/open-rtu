/* ============================================================================
 * feather2_dc_sweep.ino  —  stepped-DC sweep source for ADC characterisation
 * ----------------------------------------------------------------------------
 * SEPARATE BUILD. Do not merge into feather_2_code_mimoDEMO.ino.
 *
 * Replaces the sine test, which is not usable on this bench: Feather 2's
 * 8-bit DAC driving a 40-point LUT at AMPLITUDE 100 has a quantisation SNR
 * of 20*log10(70.7/0.289) = 47.8 dB, which is WORSE than the ~50 dB ADC
 * chain it would be measuring. Its 2 kHz staircase images at 1950 and
 * 2050 Hz also alias straight back onto 50 Hz in a 2 kHz sampler. Any
 * spectral figure through that path characterises Feather 2.
 *
 * A stepped-DC sweep has none of those problems. Every point is DC, so
 * source distortion cannot enter, and it gives what Section 2.8.6 actually
 * needs: input-referred noise as a function of input level across the whole
 * range, rather than at the bias point alone.
 *
 * No synchronisation is needed. The sweep free-runs; Feather 1 captures
 * with the UNCHANGED characterise.ino; analyse_sweep.py finds the plateaus
 * by detection. The ADC's own mean at each plateau is the x-axis, so the
 * analyser never needs to know which DAC code produced which level.
 *
 * DAC_VA is parked at midpoint throughout so the neighbouring channel is
 * not also moving.
 *
 * Flash to FEATHER 2 (COM5). Leave Feather 1 running characterise.ino with
 * MODE_ROUNDROBIN 0 and RUN_LABEL "C1_sweep".
 *
 * Jack Davies, September 2026
 * ========================================================================== */

#define DAC_VA     25
#define DAC_VB     26

#define DC_MID     127
#define STEP       4      // DAC codes per level -> 64 levels over 0..252
#define DWELL_MS   128    // 256 ADC samples per level at 2 kHz
                          // 64 levels x 128 ms = 8.192 s = one full capture

void setup() {
  Serial.begin(115200);
  delay(1000);
  dacWrite(DAC_VA, DC_MID);        // park the neighbour
  dacWrite(DAC_VB, DC_MID);
  Serial.println("# feather2_dc_sweep: DAC_VB stepped 0..252 in 4s, 128 ms dwell.");
  Serial.println("# DAC_VA parked at 127. Free-running, no sync needed.");
  Serial.println("# Capture on Feather 1 with characterise.ino, MODE_ROUNDROBIN 0.");
}

void loop() {
  for (int code = 0; code <= 252; code += STEP) {
    dacWrite(DAC_VB, code);
    delay(DWELL_MS);
  }
}
