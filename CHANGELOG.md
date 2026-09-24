# Changelog

## v1.2.0 — 2026-09 (release cited by the submitted manuscript)

Added
- `analysis/paper_tables.py` — regenerates, from the released captures, every value in
  Tables 2–8 and the per-region statistics of Sections 2.3, 3.1, 3.2 and 3.8 that no
  earlier script produced (exact Student-t TOST p-values and 90 % intervals, 95 %
  bootstrap intervals, the serial-correlation check, all five TVE regions, the
  stationary/dynamic columns of Table 4, Table 5, Table 7 on the Section 3.1 partition,
  Table 8 with its frame accounting, and the Section 3.8 repeatability figures), and
  checks each against the manuscript (`--check`). Writes `paper_tables_output.json`.
- `analysis/paper_figures.py` — redraws Figures 11–15, 17 and 18 with the manuscript's
  Tier 1/Tier 2 nomenclature, the Section 3.1 region partition and "NESO reference"
  labelling; Figure 16 is redrawn when the DMA capture is supplied with `--dma`.
- `reproduce_all.py` steps 12 (tables) and 13 (figures).
- `CHANGELOG.md`.
- `bom.csv` — the bill of materials of the paper's Table 9 as data (manifested).
- `requirements.txt` — the NumPy/SciPy/matplotlib versions that produced the manifested outputs;
  README explains that newer versions reproduce every checked value but not every sidecar byte.
- `hardware/MIMO_LCL_conditioning.net` and `hardware/README.md` — the LTspice netlist of the
  six-channel conditioning model (Supplementary Data S11; op-amp macromodel not redistributed).
- `analysis/kalman_sweep.py` and `kalman_sweep_output.json` — the Kalman Q/R sensitivity sweep of
  Section 3.5 (¼×–4× the first-principles ratio: 9.8–4.4 %, positive throughout).
- `analysis/crb_confirm.py` and `crb_confirm_output.json` — the three confirmations of the
  real-sinusoid Cramér–Rao bound stated in Section 3.5 (closed form, inversion of the 3×3
  Fisher matrix for (A, f, φ), Monte-Carlo least-squares fit with a fixed seed) at 2 kHz and
  1816.5 Hz, 40 dB, and 1816.5 Hz, 61.4 dB. Run separately; not a step of `reproduce_all.py`.
- `tost_metrics.py`: optional `--json` output path, so that runs with non-default options
  need not overwrite the manifested `captures/..._tost.json`. The default is unchanged.
- `characterisation/20260923/analysis/`: `verify_*_312.py`, `common_312.py`,
  `make_figure_19_noise_residual.py`, `make_figure_20_error_budget.py` — every value and figure of
  Section 3.12 recomputed from the released records (outputs in `out_312/`). Figures 19 and 20
  renumbered to first-citation order and regenerated from released data.

Fixed
- `tost_sweep.py`, `tost_metrics.py`: docstrings corrected. The normal approximation to
  the t distribution is mildly anti-conservative, not conservative; `tost_sweep.py`
  reproduces the 5 mHz grid crossovers of Section 3.2, while the exact Student-t p-values
  of Table 2 come from `paper_tables.py`; the bootstrap compares 95 % intervals and does
  not confirm far-tail p-values. Console text and the `--plot` figure now say Tier 1
  (±10 mHz) and Tier 2 (±100 mHz). JSON keys and values are unchanged.
- `bland_altman.py`: docstring refers to the manuscript's Figure 13 (drawn by
  `paper_figures.py`); the script's own plot labels the shaded band "Tier 2 (±100 mHz)"
  instead of an IEC 61400-21 class. The default plot file name is unchanged.
- `reproduce_all.py`: under `--skip-figures`, step 2 now passes the null device to
  `bland_altman.py --figure`, so no `Figure_10_Bland_Altman.png` is written. Without the
  flag the behaviour is unchanged.
- `tve_metrics.py`: the comment claiming that frames sharing a replay index are averaged is
  corrected (every frame is scored, the full-stream basis); the stale "§5.3.4" reference
  now points to Section 3.4; the standard is cited as IEC/IEEE 60255-118-1 (formerly IEEE
  C37.118.1-2011). No computed value changes.
- `crb_analysis.py`: docstring corrected. Its whole-replay ratios are not the Table 4
  ratios, which are stationary-block ratios from `paper_tables.py`; the only
  `crb_output.json` quantity the paper uses is the 40 dB bound, 3.5408 mHz at 1816.5 Hz;
  the `table_17_rows` key is historical and kept. The console prints f_s and SNR to one
  decimal. `crb_output.json` is unchanged.
- `characterisation/.../analyse_noise.py`, `analyse_sinad.py`: comments note that their
  printed CRB uses the complex-exponential coefficient 6, is low by √2 against the
  real-sinusoid bound, and is not used by the paper. Computation and sidecars unchanged.
- `tve_metrics.py`: comments only. The 0.7683 V magnitude reference is documented as a
  retained convention (0.2364 × 3.25 V, within 0.2 % of the whole-run mean amplitude); the
  as-built conditioning gain is 0.2240 (paper Section 2.1.2). No computed value changes.
- `analyse_v6_repeatability.py`: regions and weights now follow the Section 3.1
  partition and the Table 8 frame counts; a two-column trajectory CSV with a
  `frequency_hz` column is now accepted by `--neso`.
- Every analysis script now writes its JSON sidecar with LF line endings, so the
  manifest hashes of regenerated outputs are the same on Windows, macOS and Linux
  (previously only `inertia_output.json` was written this way). Sidecars in
  `captures/` and the root-level outputs were regenerated accordingly.
- `make_hashes.py` now manifests the root-level derived outputs (`*_output.json`),
  which the Data Availability Statement refers to.
- `tost_sweep.py` and `tve_metrics.py` recorded the absolute path of the input CSV in
  their sidecars (a Windows temp path in v1.1.0), so the manifest could never verify
  after a re-run on another machine. They now record the file name only, like the
  other scripts; the SHA-256 of the input is unchanged.
- `inertia_estimator.py`: the provenance notes in the docstring and in
  `inertia_output.json` now match Section 3.11 of the manuscript - the 0.16 Hz/s
  initial RoCoF is a nominal figure, not a published one; the stale 1378 MW note and
  the obsolete reference numbers are gone. The recovered values are unchanged, but the
  file's text changed, so its SHA-256 is now
  `91AFA54072CB60C38444D0223FC3AE2C41960CDB8703BD84F90D5DB055FCF9B2` (quoted in the Data
  Availability Statement and checked by `reproduce_all.py`).
- Figure 18, panel (a): title and marker label now say "external inputs" and "nominal
  0.16 Hz/s" rather than "published inputs", for the same reason.
- Figures 11 and 15-18: the legend and axis notation now matches the text
  (f with subscript Vb, E with subscript k) instead of the literal `f_Vb` and `E_k`. The
  static Figure 16 carries the same legend correction.
- `crb_analysis.py`, `kalman_freq.py`: console headings cite the manuscript's current table
  numbers (Table 4) instead of numbers from earlier drafts.
- Docstrings, comments and console text across `analysis/`, `characterisation/` and `hmi/`
  now cite the manuscript's final section, table, figure and reference numbers; references
  to items that no longer exist in the paper, to internal drafts and to task codes are
  replaced by plain descriptions. `characterisation/20260913/` is left exactly as archived.
  Strings written into JSON outputs are unchanged, so every released output is reproduced
  byte for byte.
- Figures 11 and 15-17 no longer carry an in-image working title (capture file name and
  summary statistics): the manuscript captions carry that information. `paper_figures.py`
  prints the working title to the console, and `--titles` draws it on the figure again. The
  static Figure 16 has had its title removed to match.
- `run_live_verification.py`: archived TVE-equivalent mean/maximum updated to the
  values printed in the manuscript (3.18 % / 8.72 %).
- `figures/`: stale v1.1.0 renders with the old figure numbering and "Class I/II"
  labels replaced by the manuscript's Figures 11–18.
- `characterisation/20260923/figures/`: the error budget is now `Figure_20_error_budget.png` and
  the noise/residual figure `Figure_19_noise_and_residual.png` (first-citation order), both drawn
  by released scripts from released records; the earlier untracked renders are removed.

## v1.1.0

Deposit of everything the Data Availability Statement names (see git history).
