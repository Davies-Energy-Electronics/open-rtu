/* ============================================================================
 * feather2_dc_hold.ino  —  hold both DACs at DC midpoint, indefinitely
 * ----------------------------------------------------------------------------
 * SEPARATE BUILD. Do not merge into feather_2_code_mimoDEMO.ino.
 *
 * Feather 2's normal firmware holds DC midpoint for 5 s at boot and then
 * starts the 50 Hz output. Five seconds is shorter than an 8.192 s
 * characterisation record, so the quiet-input captures need this instead.
 *
 * Holding the DACs at 127 reproduces exactly the condition the measurement
 * firmware calibrates its bias against ("Measuring bias (Feather 2 DC
 * midpoint)"), with the real source impedance and the real bias point still
 * in circuit. That is a better quiet input than shorting anything, because
 * nothing about the signal path changes except that the AC is absent.
 *
 * Flash to FEATHER 2 (the generator board). Leave Feather 1 alone.
 *
 * Jack Davies, September 2026
 * ========================================================================== */

#define DAC_VA  25
#define DAC_VB  26
#define DC_MID  127        // same midpoint constant the replay LUT is built around

void setup() {
  Serial.begin(115200);
  delay(1000);
  dacWrite(DAC_VA, DC_MID);
  dacWrite(DAC_VB, DC_MID);
  Serial.println("# feather2_dc_hold: DAC_VA and DAC_VB held at 127, no AC.");
  Serial.println("# Characterisation source. Reflash the demo firmware afterwards.");
}

void loop() {
  // Rewrite every second. dacWrite latches, so this is belt and braces
  // against anything else touching the pins.
  dacWrite(DAC_VA, DC_MID);
  dacWrite(DAC_VB, DC_MID);
  delay(1000);
}
