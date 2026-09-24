Open RTU - ADC chain characterisation
Bench session of 13 September 2026
=====================================================================

CONDITIONS
  Build            MIMO stripboard, unchanged since June 2026
  Ambient          21 degC
  Power            USB, both boards bus-powered from one laptop
                   (Asus Vivobook Go 14). Shared supply rail: the
                   generator and the measurement board are not
                   independently powered.
  Ports            COM5 = Feather 2 (generator), COM6 = Feather 1
                   (measurement). Board serials not marked; boards
                   identified by port assignment.
  ADC              12-bit, Arduino-ESP32 default 11 dB attenuation,
                   no analogSetPinAttenuation or analogReadResolution
                   call in any build. Linear 3.3/4095 conversion.
  Channel          Vb on GPIO39, the primary M2 metric channel.

FIRMWARE (all separate builds; neither measurement firmware was modified)
  characterise.ino        raw ADC code dump, deadline-scheduled 2 kHz
  adc_latency.ino         analogRead() duration, 20000 repetitions
  feather2_dc_hold.ino    DACs held at midpoint 127 (quiet input)
  feather_2_code_mimoDEMO.ino   production generator, for B1

CAPTURES
  A1_quiet_single.txt        quiet input, single channel
  A2_quiet_single_run2.txt   repeat of A1
  A3_quiet_single_run3.txt   repeat of A1
  D1_latency_single.txt      analogRead() latency
  D2_latency_single_run2.txt repeat of D1
  B1_ac_running.txt          production 50 Hz stimulus present

  NOTE ON NAMING. A2/A3 and D2 were captured intending round-robin and
  sweep modes, but the sketch was not re-edited between uploads, so all
  ran in the previous single-channel mode. They are named here for what
  they ACTUALLY contain, not what was intended. Treated as repeatability
  runs they are useful; they are not the round-robin or sweep tests,
  which remain undone.

RESULTS
  Converter input-referred noise, three independent runs
      0.6871, 0.6332, 0.6555 mV     mean 0.6586 mV, spread 8.2 %
      = 0.817 LSB, 61.4 dB SNR at A = 1.0928 V
      = 10.50 effective bits of a nominal 12

  analogRead() duration
      mean 39.21 us, sd 0.86 us, range 38-44 us
      (0.29 us of that sd is micros() 1 us resolution)
      Sample-instant uncertainty is therefore NOT a significant
      contributor to frequency error.

  Operating conditions, production stimulus present (B1)
      residual after per-window sinusoid fit  12.28 mV = 15.24 LSB
      that is 17.9x the quiet-input figure
      estimator dispersion on these samples   3.047 mHz
      successive-difference floor             3.100 mHz

ERROR BUDGET (frequency, at the FREQ_WIN = 200 estimator basis)
      converter noise                 1.68 mHz    9.4 % of variance
      stimulus, 8-bit DAC             2.61 mHz   22.8 %
      V2 acquisition loop             4.50 mHz   67.8 %
      -----------------------------------------------------------
      measured today, clean sampling  3.10 mHz
      archived June capture (V2 loop) 5.46 mHz

      Terms 2 and 3 are obtained by difference. Term 3 in particular
      compares today's deadline-scheduled single-channel acquisition
      against a June capture taken with V2's analogRead-plus-delay
      loop: two things differ, so 4.50 mHz is an UPPER BOUND on the
      acquisition loop's contribution, not an isolation of it.

DISCREPANCIES FOUND
  1. V2 measureFreq() calls analogRead() then delayMicroseconds(500),
     so its real per-channel period is about 540 us - roughly 1850 Hz,
     not the 2000 Hz stated in the paper and hardcoded in
     crb_analysis.py. Shifts the Cramer-Rao bound by about 8 %.
  2. The reported frequency field is quantised at 2.5 mHz, because the
     median period is an integer number of microseconds and
     df = 1e6/T^2 per us at T = 20 ms. Confirmed directly in
     replay_20260627_102507V2.csv. Not in the paper's budget.
     Firmware-fixable by accumulating the period in float.
  3. The local reproduce_all.py (16836 bytes) differs from the one
     released in the GitHub bundle (7724 bytes). The archived harness
     is not the one in use.

WHAT REMAINS UNDONE
  Round-robin noise comparison (A2 as intended) and the stepped-DC
  sweep (C1). Neither is on the critical path now that the single-
  channel noise is measured at 0.66 mV, but both would tighten the
  budget.
