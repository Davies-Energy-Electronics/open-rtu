"""
analyse_v6_repeatability.py
===========================

Reproduces the Section 3.9 run-to-run repeatability result (Roadmap C.2) from
the three raw V6 parallel-firmware captures of 24 July 2026. This is the
released, reproducible record of the analysis first performed for the paper.

It replaces run_repeatability_check.py for the V6 case: that script wraps the
SCADA-CSV toolchain and cannot parse raw V6 REPORT streams.

WHAT IT DOES
------------
1. Loads the canonical NESO 9 Aug 2019 trajectory (360 x 1 s samples), either
   from a two-column trajectory CSV or by reconstructing it from the archived
   canonical V2 SCADA capture (replay_index -> f_target_hz).
2. Parses each raw V6 capture (lines: REPORT,<millis>,...,f_Vb=<Hz>),
   placing REPORTs on a regular 20 ms grid from each capture's first record.
3. Aligns each capture to the trajectory by normalised cross-correlation of
   (f_Vb - 50) against the 20 ms-expanded (f_target - 50).
4. Scores |f_Vb - f_target| per frame; reports per-region and overall mean/max,
   Class II (<=100 mHz) and Class I (<=10 mHz) coverage.
5. Compares each run's overall mean to the archived Table 6 result
   (78.3 mHz) on BOTH bases: the full 360 s window, and re-weighted to the
   archived Table 6 regional frame counts (the archived capture's settled
   tail was truncated, so the weighted basis is the like-for-like one).
6. Writes repeatability_summary.json with all stats and input SHA-256 hashes.

USAGE
-----
    python analyse_v6_repeatability.py \
        --captures replay_20260724_130842_v6_run1.csv \
                   replay_20260724_131622_v6_run2.csv \
                   replay_20260724_132305_v6_run3.csv \
        --neso replay_20260627_102507V2.csv \
        --outdir repeatability_v6

--neso accepts either the archived canonical V2 SCADA CSV (trajectory is
reconstructed from it) or a trajectory CSV from extract_neso_window.py.

Requires numpy. Author: Jack Davies (analysis specification: paper Section 3.9).
"""

from __future__ import annotations
import argparse, csv, hashlib, json, os, re, sys
import numpy as np

ARCHIVED_OVERALL_MHZ = 78.3
TOLERANCE_MHZ = 8.0
# Table 6 regional frame counts of the archived canonical V6 capture
ARCHIVED_WEIGHTS = {"Pre-event": 2299, "RoCoF descent": 2750, "Nadir": 4000,
                    "Recovery": 4850, "Settled tail": 1997}
REGIONS = [(0, 45, "Pre-event"), (45, 100, "RoCoF descent"), (100, 180, "Nadir"),
           (180, 285, "Recovery"), (285, 360, "Settled tail")]
CADENCE_MS = 20
WINDOW_S = 360


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for c in iter(lambda: fh.read(65536), b""):
            h.update(c)
    return h.hexdigest().upper()


def load_trajectory(path):
    """Return the 360-sample trajectory (Hz). Accepts a SCADA V2 capture
    (replay_index + f_target_hz columns) or a 2-column trajectory CSV."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    cols = rows[0].keys()
    if "replay_index" in cols and "f_target_hz" in cols:
        traj = {}
        for r in rows:
            i = int(float(r["replay_index"]))
            traj.setdefault(i, float(r["f_target_hz"]))
        f = np.array([traj[i] for i in sorted(traj)])
    else:
        fcol = next(c for c in cols if c.lower() in
                    ("f_neso_hz", "f", "f_hz", "frequency"))
        f = np.array([float(r[fcol]) for r in rows])
    if len(f) != WINDOW_S:
        print(f"  warning: trajectory has {len(f)} samples, expected {WINDOW_S}",
              file=sys.stderr)
    return f


def load_capture(path):
    """Raw V6 stream -> f_Vb on a regular 20 ms grid (NaN where missing)."""
    pat = re.compile(r"^REPORT,(\d+),.*f_Vb=([0-9.]+)")
    t, f = [], []
    with open(path) as fh:
        for ln in fh:
            m = pat.match(ln)
            if m:
                t.append(int(m.group(1))); f.append(float(m.group(2)))
    t, f = np.array(t), np.array(f)
    idx = np.round((t - t[0]) / CADENCE_MS).astype(int)
    grid = np.full(idx.max() + 1, np.nan)
    grid[idx] = f
    return grid, len(t)


def align(grid, tgt_20ms):
    """Best integer 20 ms lag by normalised cross-correlation."""
    m = np.where(np.isnan(grid), 50.0, grid) - 50.0
    g = tgt_20ms - 50.0
    best_lag, best_c = 0, -9.0
    for lag in range(0, len(m) - len(g) + 1):
        seg = m[lag:lag + len(g)]
        c = float(np.dot(seg, g) /
                  (np.linalg.norm(seg) * np.linalg.norm(g) + 1e-12))
        if c > best_c:
            best_c, best_lag = c, lag
    return best_lag, best_c


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--captures", nargs="+", required=True)
    ap.add_argument("--neso", required=True)
    ap.add_argument("--outdir", default="repeatability_v6")
    args = ap.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)

    traj = load_trajectory(args.neso)
    tgt = np.repeat(traj, 1000 // CADENCE_MS)  # 20 ms expansion
    print(f"Trajectory: {len(traj)} samples, nadir {traj.min():.3f} Hz "
          f"at index {int(np.argmin(traj))}")

    runs, fulls, wtds = {}, [], []
    W, N = ARCHIVED_WEIGHTS, sum(ARCHIVED_WEIGHTS.values())
    for k, path in enumerate(args.captures, 1):
        grid, nrep = load_capture(path)
        lag, cc = align(grid, tgt)
        err = np.abs(grid[lag:lag + len(tgt)] - tgt) * 1000.0
        valid = ~np.isnan(err)
        regions = {}
        for lo, hi, name in REGIONS:
            e = err[lo * 1000 // CADENCE_MS: hi * 1000 // CADENCE_MS]
            regions[name] = {"n": int(np.sum(~np.isnan(e))),
                             "mean_mHz": round(float(np.nanmean(e)), 1),
                             "max_mHz": round(float(np.nanmax(e)), 1)}
        overall = round(float(np.nanmean(err)), 1)
        weighted = round(sum(regions[n]["mean_mHz"] * W[n] for n in W) / N, 1)
        fulls.append(overall); wtds.append(weighted)
        runs[f"run{k}"] = {
            "file": os.path.basename(path), "sha256": sha256(path),
            "n_report_lines": nrep,
            "replay_start_s": round(lag * CADENCE_MS / 1000.0, 2),
            "align_corr": round(cc, 4),
            "in_replay_frames": int(valid.sum()),
            "overall_mean_mHz": overall,
            "overall_max_mHz": round(float(np.nanmax(err)), 1),
            "classII_pct": round(100 * float(np.mean(err[valid] <= 100)), 1),
            "classI_pct": round(100 * float(np.mean(err[valid] <= 10)), 1),
            "table6_weighted_mean_mHz": weighted,
            "within_tolerance": bool(abs(weighted - ARCHIVED_OVERALL_MHZ)
                                     <= TOLERANCE_MHZ),
            "regions": regions,
        }
        print(f"\n{os.path.basename(path)}: start {lag*CADENCE_MS/1000:.2f}s, "
              f"corr {cc:.4f}, {nrep} REPORT")
        print(f"  overall {overall} mHz | weighted {weighted} mHz | "
              f"within +/-{TOLERANCE_MHZ} of {ARCHIVED_OVERALL_MHZ}: "
              f"{runs[f'run{k}']['within_tolerance']}")
        for name, s in regions.items():
            print(f"    {name:14s} n={s['n']:5d} mean={s['mean_mHz']:6.1f} "
                  f"max={s['max_mHz']:6.1f}")

    verdict = all(r["within_tolerance"] for r in runs.values())
    print(f"\nRun-to-run spread: full-window {max(fulls)-min(fulls):.1f} mHz, "
          f"weighted {max(wtds)-min(wtds):.1f} mHz")
    print(f"VERDICT: {'PASS' if verdict else 'FAIL'} - "
          f"{'all' if verdict else 'not all'} weighted means within "
          f"+/-{TOLERANCE_MHZ} mHz of {ARCHIVED_OVERALL_MHZ} mHz")

    out = {"schema_version": "1.0", "tool": "analyse_v6_repeatability.py",
           "archived_overall_mHz": ARCHIVED_OVERALL_MHZ,
           "tolerance_mHz": TOLERANCE_MHZ,
           "archived_weights": ARCHIVED_WEIGHTS,
           "neso_reference": os.path.basename(args.neso),
           "neso_sha256": sha256(args.neso),
           "verdict": "PASS" if verdict else "FAIL",
           "spread_full_window_mHz": round(max(fulls) - min(fulls), 1),
           "spread_weighted_mHz": round(max(wtds) - min(wtds), 1),
           "runs": runs}
    dst = os.path.join(args.outdir, "repeatability_summary.json")
    with open(dst, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print(f"Wrote {dst}")
    return 0 if verdict else 1


if __name__ == "__main__":
    raise SystemExit(main())
