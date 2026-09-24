#!/usr/bin/env python3
"""
kalman_sweep.py
===============

Noise-covariance sensitivity sweep for the V4 Kalman estimator
(Section 3.5 of the open-RTU paper; companion to kalman_freq.py).

kalman_freq.py derives the V4 filter's process-noise intensity q0 and
measurement-noise variance r0 from first principles (see its docstring).
This script asks how much the V4 dispersion reduction against V2 depends on
those two numbers. It reuses kalman_freq.py's loader, Q/R derivation, filter
and error statistics by import, so the filter under test is byte-for-byte the
one that produces the Table 4 V4 row.

Two sweep definitions are reported:

  (a) 2-D independent sweep
      q = q0 * a,  r = r0 * b,  with a, b each on a log grid over
      [1/4, 4] (a sixteen-fold range per axis; default 9 x 9 points,
      i.e. factors 2^(-2), 2^(-1.5), ..., 2^(+2)).
      Note: because the filter's gain depends (after the initial transient)
      only on the ratio q/r, this grid spans q/r from 1/16x to 16x of the
      first-principles ratio (a 256-fold range in the ratio).

  (b) 1-D ratio sweep
      q/r = (q0/r0) * c,  with c on a log grid over [1/4, 4]
      (a sixteen-fold range in the ratio; default 33 points). Implemented as
      q = q0 * sqrt(c), r = r0 / sqrt(c).

For every point the metric is the one kalman_freq.py reports:
    improvement_pct = (1 - std(V4 - ref) / std(V2 - ref)) * 100
over the whole replay (replay_index >= 0), with the NESO reference
substituted for f_target_hz exactly as kalman_freq.py does.

Usage:
    python kalman_sweep.py captures/replay_20260627_102507V2.csv \
        --neso trajectories/neso_trajectory.csv
    python kalman_sweep.py <csv> --neso <neso> --n2d 9 --n1d 33 \
        --json kalman_sweep_output.json

Deterministic, standard library only. Writes a JSON sidecar with LF line
endings and basenames only (no absolute paths).

Author: Jack Davies
"""

from __future__ import annotations

import argparse
import json
import math
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kalman_freq as kf  # noqa: E402  (reuse the released V4 filter)


SWEEP_LO = 0.25   # x1/4
SWEEP_HI = 4.0    # x4   -> sixteen-fold range


def log_grid(lo: float, hi: float, n: int):
    """n log-spaced factors from lo to hi inclusive (exact endpoints)."""
    if n < 2:
        raise ValueError("grid needs at least 2 points")
    a, b = math.log2(lo), math.log2(hi)
    return [2.0 ** (a + (b - a) * k / (n - 1)) for k in range(n)]


def improvement(rows, q: float, r: float, v2_std: float) -> tuple[float, float]:
    """Return (V4 std in mHz, improvement % vs V2) for one (q, r) pair."""
    filtered = kf.kalman_run(rows, q, r)
    stats, _ = kf.error_stats(filtered, rows)
    return stats["std_mhz"], (1.0 - stats["std_mhz"] / v2_std) * 100.0


def summarise(values):
    return {"min_pct": min(values),
            "max_pct": max(values),
            "median_pct": statistics.median(values),
            "mean_pct": statistics.mean(values),
            "n_points": len(values)}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", help="canonical V2 replay CSV")
    ap.add_argument("--neso", default="neso_trajectory.csv",
                    help="NESO reference CSV (same as kalman_freq.py)")
    ap.add_argument("--n2d", type=int, default=9,
                    help="grid points per axis for the 2-D sweep (default 9)")
    ap.add_argument("--n1d", type=int, default=33,
                    help="grid points for the Q/R-ratio sweep (default 33)")
    ap.add_argument("--json", default="kalman_sweep_output.json")
    args = ap.parse_args(argv)

    # ---- Load exactly as kalman_freq.main() does --------------------------
    rows = kf.load_capture(args.csv)
    if not rows:
        raise SystemExit(f"No replay rows in {args.csv}")
    neso = kf.load_neso(args.neso)
    if neso:
        for d in rows:
            if d["i"] in neso:
                d["f_target"] = neso[d["i"]]
    r0 = kf.derive_R(rows)
    q0 = kf.derive_Q(rows, neso)

    v2_stats, _ = kf.error_stats([d["f_Vb"] for d in rows], rows)
    v2_std = v2_stats["std_mhz"]
    base_std, base_imp = improvement(rows, q0, r0, v2_std)

    # ---- (a) 2-D independent sweep ---------------------------------------
    fac2 = log_grid(SWEEP_LO, SWEEP_HI, args.n2d)
    grid = []
    for a in fac2:
        for b in fac2:
            s, imp = improvement(rows, q0 * a, r0 * b, v2_std)
            grid.append({"q_factor": a, "r_factor": b,
                         "qr_ratio_factor": a / b,
                         "v4_std_mhz": s, "improvement_pct": imp})
    sum2 = summarise([g["improvement_pct"] for g in grid])
    # Subset of the 2-D grid whose q/r ratio stays within the 1/4..4 band
    in_band = [g["improvement_pct"] for g in grid
               if SWEEP_LO - 1e-12 <= g["qr_ratio_factor"] <= SWEEP_HI + 1e-12]
    sum2_band = summarise(in_band)

    # ---- (b) 1-D Q/R-ratio sweep -----------------------------------------
    fac1 = log_grid(SWEEP_LO, SWEEP_HI, args.n1d)
    ratio = []
    for c in fac1:
        s, imp = improvement(rows, q0 * math.sqrt(c), r0 / math.sqrt(c), v2_std)
        ratio.append({"qr_ratio_factor": c, "v4_std_mhz": s,
                      "improvement_pct": imp})
    sum1 = summarise([p["improvement_pct"] for p in ratio])

    # ---- Report -----------------------------------------------------------
    print()
    print("=" * 72)
    print("  V4 KALMAN NOISE-COVARIANCE SENSITIVITY SWEEP (Section 3.5)")
    print("=" * 72)
    print(f"  Input CSV:  {os.path.basename(args.csv)}   frames: {len(rows)}")
    print(f"  q0 = {q0:.6e}   r0 = {r0:.6e} Hz^2   (first-principles)")
    print(f"  V2 std = {v2_std:.2f} mHz   V4 std @ (q0,r0) = {base_std:.2f} mHz"
          f"   improvement = {base_imp:.2f} %")
    print()
    print(f"  (a) 2-D independent sweep, q and r each x1/4..x4, "
          f"{args.n2d}x{args.n2d} log grid")
    print("      rows: q factor, cols: r factor, cells: improvement %")
    print("      " + "".join(f"{b:>7.3g}" for b in fac2))
    for i, a in enumerate(fac2):
        cells = grid[i * len(fac2):(i + 1) * len(fac2)]
        print(f"  {a:>5.3g} " + "".join(f"{c['improvement_pct']:>7.2f}" for c in cells))
    print(f"      min {sum2['min_pct']:.2f}  max {sum2['max_pct']:.2f}  "
          f"median {sum2['median_pct']:.2f}  (n={sum2['n_points']})")
    print(f"      subset with q/r within x1/4..x4: min {sum2_band['min_pct']:.2f}  "
          f"max {sum2_band['max_pct']:.2f}  median {sum2_band['median_pct']:.2f}"
          f"  (n={sum2_band['n_points']})")
    print()
    print(f"  (b) Q/R-ratio sweep, ratio x1/4..x4, {args.n1d}-point log grid")
    for p in ratio[::max(1, (len(ratio) - 1) // 8)]:
        print(f"      q/r x{p['qr_ratio_factor']:<7.4g} -> "
              f"{p['improvement_pct']:6.2f} %")
    print(f"      min {sum1['min_pct']:.2f}  max {sum1['max_pct']:.2f}  "
          f"median {sum1['median_pct']:.2f}  (n={sum1['n_points']})")
    print("=" * 72)

    out = {
        "schema_version": "1.0",
        "tool": "kalman_sweep.py",
        "source_file": os.path.basename(args.csv),
        "reference_file": os.path.basename(args.neso) if neso else None,
        "filter": "kalman_freq.kalman_run (imported, unmodified)",
        "metric": "improvement_pct = (1 - std(V4-ref)/std(V2-ref)) * 100, "
                  "whole replay (replay_index >= 0)",
        "n_frames": len(rows),
        "baseline": {"q0": q0, "r0_hz2": r0, "v2_std_mhz": v2_std,
                     "v4_std_mhz": base_std, "improvement_pct": base_imp},
        "sweep_2d_independent": {
            "definition": "q = q0*a, r = r0*b; a, b each log-spaced over "
                          "[1/4, 4] (q/r therefore spans 1/16..16 x q0/r0)",
            "factors": fac2,
            "summary": sum2,
            "summary_qr_ratio_within_quarter_to_4x": sum2_band,
            "points": grid,
        },
        "sweep_qr_ratio": {
            "definition": "q/r = (q0/r0)*c, c log-spaced over [1/4, 4]; "
                          "q = q0*sqrt(c), r = r0/sqrt(c)",
            "factors": fac1,
            "summary": sum1,
            "points": ratio,
        },
    }
    with open(args.json, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")
    print(f"  -> {os.path.basename(args.json)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
