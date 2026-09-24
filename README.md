# Open RTU — an open reference instrument for grid-side compliance monitoring

Firmware, conditioning model, bill of materials, SCADA HMI, the complete analysis
toolchain and every capture reported in:

> Davies, J., Hettiarachchige Don, A., Kang, L. (2026). *An Open, Formally-Derived
> Reference Instrument for Grid-Side Compliance Monitoring of Offshore Wind Turbines:
> Replay Validation Against NESO Frequency-Event Records.* MDPI Sensors.
> [DOI on acceptance]

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
├── hashes.txt         SHA-256 manifest of every released file
├── .gitattributes     line-ending conversion OFF — see below
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

Live-hardware verification (`run_live_verification.py`) is released here and run
separately; it needs the 7 July 2026 live capture rather than the canonical replay.

Two figures have moved since the first release, and the reasons are in the source
rather than left for a reader to reconstruct:

- **The Cramér–Rao bound** was first computed with the Rife–Boorstyn coefficient 6, the
  complex-exponential form. A real sampled sinusoid takes 12, so the bound is larger by
  exactly √2. Confirmed twice — in closed form, and by inverting the 3×3 Fisher
  information matrix for (A, f, φ).
- **The sampling rate in that bound** was an assumed 2000 Hz. It is a measured 1816.5 Hz
  for the canonical build. The bound is linear in the rate at fixed window length, so it
  falls to 3.541 mHz and every ratio in Table 4 rises — in the conservative direction.

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

## Citation

Please cite the paper and the deposit. Issues, questions and corrections are welcome.
Pull requests against the analysis scripts are accepted under the same licence; pull
requests modifying the canonical capture are not.
