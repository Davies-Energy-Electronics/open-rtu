# Section 3.12 - direct characterisation of the acquisition chain

Every raw code and timing record quoted in Section 3.12 of the manuscript,
the characterisation firmware sketches that produced them, and the analysis
scripts that produce every value and figure quoted there: Figure 19 (noise
and residual), Figure 20 (error budget) and an additional figure,
Figure_21_noise_vs_code.png (noise against code), which is not in the paper. Figures 19 and 20 were
numbered the other way round in earlier drafts of the manuscript.

No measurement firmware was modified for any characterisation reported here.
The sketches in each firmware/ folder are separate builds, so the provenance
of every archived capture under captures/ at the repository root is unaffected.

## 20260913/

The bench session of 13 September 2026, exactly as it was archived on the day,
including its own README.txt and hashes.txt. The three single-channel
quiet-input records here - A1_quiet_single.txt, A2_quiet_single_run2.txt and
A3_quiet_single_run3.txt - are the records whose mean of 0.659 mV is the
figure the uncertainty budget of Section 2.8 carries and that Figure 19(a) is
drawn from. Its README records the conditions and, honestly, the naming
error that made A2 and A3 repeats of A1 rather than the round-robin and
sweep tests they were intended to be.

The signal-to-noise ratio of 61.4 dB quoted in Sections 2.8 and 3.12 is
referred to the conditioned amplitude measured on this bench, A = 1.0928 V,
that is 0.7727 V RMS against 0.6586 mV of noise.

## 20260923/

The bench sessions of 23 and 24 September 2026, which completed what the
13 September session left undone.

  captures/   A1_*        single-channel quiet input
              A2_quiet_roundrobin_run1..3
                          round-robin quiet input, the multiplexer comparison
                          (run1 warmed, run2/3 cold; a second warmed record
                          was overwritten during capture, see log_capture.py)
              A_b*, B_b*, C_b*, D_b*
                          the four-mode timing walk, three records of 16,384
                          samples in each mode
              F_*, G_*    the fourth mode repeated with the production stimulus
              C0_*        conditioned amplitude
              B1_*, ac.json
                          production stimulus present, single-tone fit residual
              C1_sweep_16step_run*
                          the 1000 ms stepped-DC staircase, which measures the
                          high-pass transient rather than noise against code
              C1_short_run*, sweep_noise.json
                          the 25 ms staircase, which resolves fifteen of the
                          sixteen levels (the top two merge) and yields the
                          upper bound of about 2.4x
              D1_, D2_    conversion latency
              runE_v2.txt, runE.json
                          the unmodified canonical build re-measured on
                          23 September (7.97 mHz over 305 windows)
              runA_*      working runs from the same sessions, retained for
                          completeness; none is quoted in Section 3.12
              A1_quiet_single.txt, A2_quiet_roundrobin.txt, C1_sweep.txt
                          byte-identical copies of the three 13 September
                          quiet-input records, kept so that earlier hashes
                          still verify; they are not additional records and
                          A2_quiet_roundrobin.txt is not a round-robin record

  firmware/   characterise.ino         raw code dump, deadline-scheduled 2 kHz
              characterise_v2loop.ino  the four-mode acquisition-loop walk
              characterise_sweep.ino   32,768-sample capture for the staircase
              feather2_dc_hold.ino     generator held at DC midpoint
              feather2_dc_sweep.ino    sixteen-level staircase, 25 ms dwell

Two files from these sessions are deliberately omitted:
C1_short_STALELABEL_do_not_archive.txt, whose capture banner did not match the
firmware mode actually running, and runB_1.txt, an aborted 262-byte run.

## Scripts added in release v1.2.0 (20260923/analysis/)

  common_312.py                    shared loaders (imports the released ones)
  verify_quiet_noise_312.py        0.659/0.729/0.646 mV, round robin 0.911/0.736 mV,
                                   multiplexer term 0.546/0.354 mV, 1.13/1.24
  verify_stimulus_residual_312.py  12.28 mV residual, 98.7 % harmonic, 1.57 % THD,
                                   1.399 mV broadband; distortion simulation
                                   (+0.27 mHz dispersion, mean < 0.03 mHz)
  verify_budget_312.py             1.70/3.60/3.08/5.06/6.46 mHz, timing walk,
                                   crossing counts, 7.97 mHz re-measurement
  verify_local_poly_312.py         local-polynomial estimator (x4, 550-950 Hz)
  verify_staircase_312.py          1.00 codes at bias, 1.6-2.3 codes elsewhere,
                                   r = +0.64 / +0.03, bound about 2.4
  make_figure_19_noise_residual.py Figure 19
  make_figure_20_error_budget.py   Figure 20
Outputs are written to 20260923/analysis/out_312/ with fixed seeds.

## Verifying

    python analysis/make_hashes.py --check

from the repository root covers this directory along with everything else.
The 13 September bundle additionally carries its own hashes.txt, which still
verifies inside 20260913/.
