# Open RTU — an open reference instrument for grid-side compliance monitoring

Firmware, bill of materials (`bom.csv`), SCADA HMI, the complete analysis
toolchain and every capture reported in:

> Davies, J., Kang, L., Hettiarachchige Don, A. (2026). *Replay Validation Against NESO
> Records: An Open Reference Instrument for Offshore Wind Frequency Compliance.*
> Submitted to MDPI Sensors. [DOI on acceptance]

This is release **v1.2.0**, the release the submitted manuscript cites. What changed
since v1.1.0 is listed in `CHANGELOG.md`.

Archived at Zenodo, concept DOI [10.5281/zenodo.21575262](https://doi.org/10.5281/zenodo.21575262),
under GPL-3.0.

## One command reproduces the paper

```bash
git clone https://github.com/Davies-Energy-Electronics/open-rtu.git
cd open-rtu
python analysis/reproduce_all.py
```

Expected last line: `ALL CHECKS PASSED`.

**Prerequisites:** Python 3.10+, `numpy`, `scipy`, `matplotlib`. Add `--skip-figures`
to run the checks without writing PNGs.

**Library versions.** The manifested outputs were produced with Python 3.11.15,
NumPy 2.4.4, SciPy 1.17.1 and matplotlib 3.10.9, pinned in `requirements.txt`. With
those versions every derived JSON and CSV file regenerates byte for byte:

```bash
python -m venv .venv && .venv/bin/pip install -r requirements.txt   # Windows: .venv\Scripts\pip
```

With newer versions, every value `reproduce_all.py` checks still reproduces within the
tolerance it prints. However, the bootstrap- and SciPy-derived sidecars can then differ
in trailing digits, so `git status` lists them as modified. That is expected, and it is
why the checks compare values rather than bytes. The one output hash quoted in the
paper, `inertia_output.json`, has been confirmed to reproduce byte for byte on NumPy 2.5.3 /
SciPy 1.18.1 under Windows as well.

`reproduce_all.py` resolves its own inputs, so it runs correctly from the repository
root, from inside `analysis/`, or from anywhere else. It does not care where you
start it.

## Layout

```
open-rtu/
├── analysis/          reproduce_all.py and the analysis toolchain
├── captures/          raw replay captures and their derived sidecars
├── trajectories/      NESO reference trajectories, CSV and .h
├── firmware/          Feather 1 measurement builds, Feather 2 generator
├── hmi/               SCADA_HMI.py — DNP3 ingest, live traces, capture
├── figures/           regenerable paper figures
├── characterisation/  Section 3.12 acquisition-chain characterisation
├── hardware/          LTspice netlist of the conditioning model (Data S11)
├── hashes.txt         SHA-256 manifest of every released file
├── .gitattributes     line-ending conversion OFF — see below
├── CITATION.cff       citation metadata (read by GitHub and Zenodo)
├── CHANGELOG.md
├── README.md
└── LICENSE            GPL-3.0
```

**The canonical capture is locked.** `captures/replay_20260627_102507V2.csv` hashes to
`dadbcfe247dd4e538b71a891b2b26f070d02c42e3dd4622b5d2125aff169560f`, step 0 refuses to
run if that has moved, and pull requests modifying it are rejected on integrity grounds.

**`.gitattributes` is not cosmetic.** Without it git rewrites the JSON sidecars to CRLF
on checkout for Windows users and the derived-output hash quoted in the paper's Data
Availability Statement fails on a clean clone — which is precisely the failure the
statement exists to prevent.

## What the harness checks

It deletes the outputs it is about to verify before it runs, so a stale file from an
earlier run cannot be mistaken for a fresh computation. That guard was added on
8 August 2026 after two scripts were found to write their sidecar only when `--json`
was passed, meaning the wrapper had been validating files on disk rather than new
results.

| Step | Reproduces | Held to |
|---|---|---|
| 0 | Canonical CSV integrity | SHA-256 `dadbcfe2…69560f` |
| 1 | TOST equivalence + bootstrap — §3.2, Table 2 | per-region crossovers |
| 2 | Bland–Altman — §3.3, Figure 13 | aggregate RMSE 208.66 mHz ± 1.0 |
| 3 | TVE-equivalent statistics — §3.4, Table 3 | mean 3.18 % ± 0.02, max 8.72 % ± 0.02 |
| 4 | TOST sensitivity sweep — §3.2, Figure 12 | crossover per region |
| 5 | Replay metrics — §3.1, Figure 11 | max RoCoF error 0.721 Hz/s ± 0.005 |
| 6 | Runtime characterisation — §3.6 | module profile |
| 7 | V3 hybrid — §3.5, Table 4 | std 116.75 mHz ± 1.0 |
| 8 | V4 Kalman — §3.5, Table 4 | std 194.25 mHz ± 1.0 |
| 9 | Swing-equation inversion — §3.11, Figure 18 | 231.4 GVA·s ± 0.5 |
| 10 | Cramér–Rao bound — §3.5, Table 4 ratios | 3.541 mHz ± 0.01 |
| 11 | Derived-output hashes — Data Availability Statement | `inertia_output.json` |
| 12 | Per-region tables — Tables 2–8, §2.3, §3.1, §3.2, §3.8 (`paper_tables.py --check`) | every tabulated value, at the precision printed |
| 13 | Figures 11–15, 17, 18 (`paper_figures.py`) | written to `figures_regenerated/`; not hashed (see below) |

Live-hardware verification (`run_live_verification.py`) is released here and run
separately; it needs the 7 July 2026 live capture rather than the canonical replay.

`analysis/crb_confirm.py` is also run separately. It reproduces the three independent
evaluations of the Cramér–Rao bound stated in Section 3.5 — the closed form of
`crb_analysis.py`, numerical inversion of the 3×3 Fisher information matrix for
(A, f, φ), and a Monte-Carlo single-tone least-squares fit (fixed seed) — at 2 kHz and
1816.5 Hz, 40 dB, and at 1816.5 Hz, 61.4 dB. It writes `crb_confirm_output.json` at the
repository root and takes a few seconds:

```bash
python analysis/crb_confirm.py
```

## Figures and naming

`figures/` holds the paper's Figures 11–18 under the manuscript's numbering. Figures
11–15, 17 and 18 are redrawn by step 13 from the hashed captures into
`figures_regenerated/`, which is not manifested because a rendered PNG is not
byte-reproducible across plotting-library versions; compare the two folders by eye. Figure 16 (the V7.2
DMA replay) is a static file: its raw capture is available from the corresponding
author on request, and `paper_figures.py --dma <capture>` redraws it when supplied.
Figures 19 and 20, and an additional figure not in the paper
(`Figure_21_noise_vs_code.png`), live in `characterisation/20260923/figures/`.

The paper calls its two accuracy thresholds **Tier 1** (±10 mHz) and **Tier 2**
(±100 mHz). Some JSON sidecars written by the v1.0/v1.1 scripts still use the keys
`Class I` and `Class II` for the same two thresholds; the keys are kept so that
existing outputs remain comparable.

Two figures have moved since the first release, and the reasons are in the source
rather than left for a reader to reconstruct:

- **The Cramér–Rao bound** was first computed with the Rife–Boorstyn coefficient 6, the
  complex-exponential form. A real sampled sinusoid takes 12, so the bound is larger by
  exactly √2. Confirmed three ways — in closed form, by inverting the 3×3 Fisher
  information matrix for (A, f, φ), and by Monte-Carlo simulation of a least-squares
  fit (`analysis/crb_confirm.py`).
- **The sampling rate in that bound** was an assumed 2000 Hz. The canonical loop measures
  1816.5 Hz on the characterisation sketch that emulates it (Section 3.12), while the
  build's own runtime profile implies at least 1849.7 Hz. The bound is linear in the rate
  at fixed window length; the paper uses 1816.5 Hz, the conservative end, so the bound is
  3.541 mHz (3.61 mHz at 1849.7 Hz) and every ratio in Table 4 is, if anything, overstated.

## If a check fails

- **Step 0** — the CSV was corrupted in download, or line-ending conversion is on.
  Confirm `.gitattributes` came down with the clone.
- **Step 11 on a Windows clone** — same cause. `git config core.autocrlf false`, delete
  the clone, clone again.
- **A numerical check off by a wide margin** — a script ran against a stale sidecar. The
  harness clears them, so this means one was written outside it. Delete every
  `*_tost.json`, `*_tve.json`, `*_bland_altman.json`, `*_tost_sweep.json`,
  `*_metrics.json`, `kalman_output.json`, `inertia_output.json`, `crb_output.json` and
  `hybrid_output.json`, then rerun.
- **A check off by a hair** — please report it. The tolerances are set at the precision
  the paper quotes, and a platform-dependent difference in the last digit is worth
  knowing about.

## Verifying the manifest

```bash
python analysis/make_hashes.py --check
```

Recomputes every entry in `hashes.txt` and reports any file that has changed, is
missing, or is present but unlisted.

## Where the paper's Supplementary Materials are in this repository

| Supplementary item | Location |
|---|---|
| Software S8 — V6 parallel firmware source | `firmware/Feather1codeMimoDemoV6_PARALLEL/` |
| Software S9 — V7–V7.2 DMA firmware sources | supplied with the paper, not in this repository |
| Data S10 — multi-event trajectory file (December 2023) | `trajectories/neso_trajectory_December2023Peak.csv` |
| Data S11 — LTspice netlist of the conditioning model | `hardware/MIMO_LCL_conditioning.net` |
| Documents S3–S7 — firmware, TOST, Bland–Altman, TVE and CRB documentation | describe the scripts in `analysis/` and the captures in `captures/` |

The NESO trajectories were extracted with `analysis/extract_neso_window.py` from the
August 2019 (`f-2019-8.zip`) and December 2023 (`fnew-2023-12.csv`) monthly files of the
NESO System Frequency dataset (https://www.neso.energy/data-portal/system-frequency-data).
The source files are NESO's and are not redistributed here.

**Known state of the bench board.** On the three-phase board the first-buffer inputs of
the Va and Vc chains are unconnected, so only Vb carries the stimulus; every result in the
paper is computed on Vb (Section 2.5). The `f_Va`, `f_Vc`, `THD_Va`, `THD_Vc` and alarm
columns of the released three-phase captures therefore reflect mains pickup on those two
channels, not the replay.

## Citation

Please cite the paper and the deposit. Issues, questions and corrections are welcome.
Pull requests against the analysis scripts are accepted under the same licence; pull
requests modifying the canonical capture are not.
