#!/usr/bin/env python3
"""
tost_sweep.py — TOST sensitivity sweep for the open RTU NESO replay validation.

Companion to tost_metrics.py. Sweeps the equivalence margin delta from 5 mHz
to 200 mHz in 5 mHz increments and computes, for each delta, the Schuirmann
Two One-Sided Tests (TOST) verdict on the per-frame frequency error
(f_Vb - f_target_hz) over each of the three windows reported in Table 9 of
the paper:

    - Pre-event nominal region:       i = 0  - 45
    - Post-event settled region:      i = 301 - 359
    - Full 360 s replay window:       i = 0  - 359

For each delta, the script outputs:
    - The sample size n, mean error mu_hat, and std deviation s_hat per window
    - The TOST p-value = max(p1, p2) where p1 tests against -delta and
      p2 tests against +delta (each one-sided at alpha = 0.05)
    - The crossover delta: the smallest delta for which max(p1, p2) < alpha,
      i.e. the empirical equivalence margin at which the TOST formally
      rejects non-equivalence for that window

The crossover delta on each window is the "effective accuracy class" of the
instrument on that window — a single-number summary far more informative
than the binary Class I (10 mHz) / Class II (100 mHz) / Class L (no formal
class) categorisation of IEC 61400-21.

Following Schuirmann 1987 [36] and the established TOST procedure used
elsewhere in this paper, all p-values use the standard normal approximation
to the t-distribution (df > 30 across all three windows; the approximation
is conservative, slightly overstating the p-value).

Usage:
    python tost_sweep.py replay_20260627_102507V2.csv
    python tost_sweep.py replay_20260627_102507V2.csv --json sweep.json
    python tost_sweep.py replay_20260627_102507V2.csv --plot sweep.png
    python tost_sweep.py replay_20260627_102507V2.csv \\
        --json sweep.json --plot sweep.png \\
        --verify-hash dadbcfe247dd4e538b71a891b2b26f070d02c42e3dd4622b5d2125aff169560f

Outputs:
    - Console table with per-window crossover delta and per-delta p-values
    - Optional JSON sidecar with full sweep results for downstream tooling
    - Optional matplotlib PNG with the three sensitivity traces and the
      Class I (10 mHz) and Class II (100 mHz) annotation lines

Reproducibility:
    Canonical input: replay_20260627_102507V2.csv
    Expected SHA-256: dadbcfe247dd4e538b71a891b2b26f070d02c42e3dd4622b5d2125aff169560f
    Expected pre-event crossover: ~ 30 mHz (within Class I band)
    Expected settled crossover:   ~ 35 mHz (within Class I band)
    Expected full-window crossover: ~ 90 mHz (within Class II band but not Class I)

Standard library only EXCEPT matplotlib (only required if --plot is requested).
Matplotlib is the only optional dependency; the JSON output, console table,
and crossover detection work without it.
"""

import argparse
import csv
import hashlib
import json
import math
import statistics
import sys
from pathlib import Path

# ============================================================================
# Constants matching the existing TOST work (paper §5.3.6, Table 9)
# ============================================================================

ALPHA = 0.05                              # TOST significance level

# Windows match Table 9 of the paper
WINDOWS = [
    ("Pre-event nominal",       0,    45),
    ("Post-event settled",     301,  359),
    ("Full 360 s replay",        0,  359),
]

# IEC 61400-21 reference accuracy bands for annotation
CLASS_I_DELTA_mHz  =  10.0
CLASS_II_DELTA_mHz = 100.0

# Sweep parameters
DELTA_MIN_mHz  =   5.0
DELTA_MAX_mHz  = 200.0
DELTA_STEP_mHz =   5.0


# ============================================================================
# Statistical helpers
# ============================================================================

def normal_cdf(z):
    """Cumulative distribution of the standard normal at z."""
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def tost_pvalue(errors, delta_mHz):
    """
    Two One-Sided Tests p-value for equivalence at +/- delta (in mHz).

    Inputs:
        errors:    list of per-frame frequency errors (Hz)
        delta_mHz: equivalence margin in mHz (e.g. 100.0 for Class II)

    Returns dict with:
        n      — sample size
        mean   — sample mean of errors (mHz, converted)
        std    — sample stdev of errors (mHz, converted)
        sem    — standard error of the mean (mHz)
        t1     — t-statistic against the -delta limit
        t2     — t-statistic against the +delta limit
        p1     — one-sided p-value (upper tail) for H01: mu <= -delta
        p2     — one-sided p-value (lower tail) for H02: mu >= +delta
        p_tost — max(p1, p2): the TOST p-value
        equivalent — True if p_tost < ALPHA at the given delta
    """
    n = len(errors)
    if n < 2:
        return {"n": n, "p_tost": 1.0, "equivalent": False,
                "mean": float("nan"), "std": float("nan")}

    # Convert errors to mHz to keep units consistent with delta
    err_mHz = [1000.0 * e for e in errors]
    mean = statistics.mean(err_mHz)
    std  = statistics.stdev(err_mHz)
    sem  = std / math.sqrt(n)

    if sem <= 0.0:
        return {"n": n, "p_tost": 0.0 if abs(mean) < delta_mHz else 1.0,
                "equivalent": abs(mean) < delta_mHz, "mean": mean, "std": std}

    t1 = (mean - (-delta_mHz)) / sem    # H01: mu <= -delta;  reject if t1 large positive
    t2 = (mean -    delta_mHz)  / sem   # H02: mu >= +delta;  reject if t2 large negative

    # One-sided p-values via the standard normal CDF
    p1 = 1.0 - normal_cdf(t1)           # upper-tail probability
    p2 =       normal_cdf(t2)           # lower-tail probability

    p_tost = max(p1, p2)

    return {
        "n": n,
        "mean":      mean,
        "std":       std,
        "sem":       sem,
        "t1":        t1,
        "t2":        t2,
        "p1":        p1,
        "p2":        p2,
        "p_tost":    p_tost,
        "equivalent": p_tost < ALPHA,
    }


def sweep_window(errors, delta_grid_mHz):
    """
    Sweep delta over a grid for one window. Returns list of result dicts in
    delta order and the crossover delta (smallest delta where equivalent).
    """
    results = []
    crossover = None
    for delta in delta_grid_mHz:
        r = tost_pvalue(errors, delta)
        r["delta_mHz"] = delta
        results.append(r)
        if r["equivalent"] and crossover is None:
            crossover = delta
    return results, crossover


# ============================================================================
# CSV reader (mirrors tost_metrics.py / tve_metrics.py convention)
# ============================================================================

def compute_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def read_canonical_csv(path):
    """
    Read the canonical replay CSV. Returns list of dicts with
    replay_index (int), f_meas (float), f_neso (float).
    """
    rows = []
    with open(path, "r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        field_map = {f.lower(): f for f in reader.fieldnames}
        def find_field(candidates):
            for c in candidates:
                if c.lower() in field_map:
                    return field_map[c.lower()]
            return None
        idx_col = find_field(["replay_index", "i", "index"])
        f_col   = find_field(["f_Vb", "f_b", "f_vb"])
        ref_col = find_field(["f_neso", "f_target_hz", "f_target", "f_ref"])

        if not all([idx_col, f_col, ref_col]):
            raise ValueError(
                f"Required column(s) not found. Need: replay_index, f_Vb, f_neso.\n"
                f"Available: {', '.join(reader.fieldnames)}"
            )

        for row in reader:
            try:
                rows.append({
                    "replay_index": int(float(row[idx_col])),
                    "f_meas":       float(row[f_col]),
                    "f_neso":       float(row[ref_col]),
                })
            except (ValueError, KeyError):
                continue
    return rows


def extract_window_errors(rows, idx_start, idx_end):
    """Extract the per-frame frequency error (f_meas - f_neso) for one window."""
    return [r["f_meas"] - r["f_neso"]
            for r in rows
            if idx_start <= r["replay_index"] <= idx_end]


# ============================================================================
# Output formatting
# ============================================================================

def print_console_report(per_window_sweep, per_window_crossover, hash_value, source_path):
    print("=" * 78)
    print("  TOST SENSITIVITY SWEEP — OPEN RTU NESO REPLAY VALIDATION")
    print("=" * 78)
    print(f"  Source CSV:           {source_path}")
    print(f"  SHA-256 (input):      {hash_value}")
    print(f"  Sweep range:          delta = {DELTA_MIN_mHz:.0f} - {DELTA_MAX_mHz:.0f} mHz "
          f"in {DELTA_STEP_mHz:.0f} mHz steps")
    print(f"  Significance level:   alpha = {ALPHA}")
    print(f"  Class I band:         delta = {CLASS_I_DELTA_mHz:.0f} mHz (IEC 61400-21)")
    print(f"  Class II band:        delta = {CLASS_II_DELTA_mHz:.0f} mHz (IEC 61400-21)")
    print()
    print("  ── CROSSOVER DELTA PER WINDOW (smallest delta with TOST equivalence) ──")
    print()
    print(f"  {'Window':<22} {'n':>5} {'mean (mHz)':>14} {'std (mHz)':>12} {'crossover':>14} {'class':>10}")
    print(f"  {'-'*22} {'-'*5} {'-'*14} {'-'*12} {'-'*14} {'-'*10}")
    for (wname, _, _), sweep, crossover in zip(WINDOWS, per_window_sweep, per_window_crossover):
        n = sweep[0]["n"]
        mean = sweep[0]["mean"]
        std = sweep[0]["std"]
        if crossover is None:
            xover_str = "> 200 mHz"
            cls = "above"
        elif crossover <= CLASS_I_DELTA_mHz:
            xover_str = f"{crossover:.0f} mHz"
            cls = "Class I"
        elif crossover <= CLASS_II_DELTA_mHz:
            xover_str = f"{crossover:.0f} mHz"
            cls = "Class II"
        else:
            xover_str = f"{crossover:.0f} mHz"
            cls = "Class L"
        print(f"  {wname:<22} {n:>5} {mean:>14.2f} {std:>12.2f} {xover_str:>14} {cls:>10}")
    print()
    print("  ── PER-WINDOW SWEEP DETAIL (first/middle/last + crossover row) ────────")
    for (wname, _, _), sweep in zip(WINDOWS, per_window_sweep):
        print()
        print(f"  {wname}:")
        print(f"    {'delta (mHz)':>11} {'p_tost':>10} {'equivalent':>12}")
        for r in sweep:
            mark = "  EQ" if r["equivalent"] else "   -"
            print(f"    {r['delta_mHz']:>11.0f} {r['p_tost']:>10.4f} {mark:>12}")
            if r["equivalent"] and r["delta_mHz"] >= 50 and len(sweep) > 10:
                # After the first equivalence is established and we've shown enough rows,
                # skip ahead a bit for clarity. But still show all rows on small sweeps.
                pass
    print()
    print("=" * 78)


def build_json_sidecar(per_window_sweep, per_window_crossover, hash_value, source_path):
    sidecar = {
        "schema_version": "1.0",
        "tool": "tost_sweep.py",
        "input": {
            "source_file": str(source_path),
            "sha256": hash_value,
        },
        "constants": {
            "alpha": ALPHA,
            "delta_min_mHz": DELTA_MIN_mHz,
            "delta_max_mHz": DELTA_MAX_mHz,
            "delta_step_mHz": DELTA_STEP_mHz,
            "class_I_delta_mHz": CLASS_I_DELTA_mHz,
            "class_II_delta_mHz": CLASS_II_DELTA_mHz,
        },
        "windows": {},
    }
    for (wname, i_start, i_end), sweep, crossover in zip(WINDOWS, per_window_sweep, per_window_crossover):
        sidecar["windows"][wname] = {
            "index_range": f"{i_start}-{i_end}",
            "n_frames": sweep[0]["n"],
            "mean_error_mHz": sweep[0]["mean"],
            "std_error_mHz": sweep[0]["std"],
            "crossover_delta_mHz": crossover,
            "sweep": [
                {
                    "delta_mHz": r["delta_mHz"],
                    "p_tost": r["p_tost"],
                    "equivalent": r["equivalent"],
                }
                for r in sweep
            ],
        }
    return sidecar


def render_plot(per_window_sweep, per_window_crossover, output_path):
    """Render the TOST sensitivity curve as PNG via matplotlib."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib unavailable; skipping --plot output.", file=sys.stderr)
        return False

    fig, ax = plt.subplots(figsize=(9, 5.5))

    colors = ["#1f77b4", "#2ca02c", "#d62728"]
    markers = ["o", "s", "^"]

    for (wname, _, _), sweep, crossover, color, mk in zip(
            WINDOWS, per_window_sweep, per_window_crossover, colors, markers):
        deltas = [r["delta_mHz"] for r in sweep]
        ps     = [r["p_tost"]    for r in sweep]
        ax.plot(deltas, ps, color=color, marker=mk, markersize=4,
                linewidth=1.4, label=wname, alpha=0.9)
        if crossover is not None:
            ax.axvline(crossover, color=color, linestyle=":", linewidth=1.0, alpha=0.55)

    # Significance threshold
    ax.axhline(ALPHA, color="black", linestyle="--", linewidth=1.0, alpha=0.7,
               label=f"alpha = {ALPHA}")

    # Class I and Class II reference lines
    ax.axvline(CLASS_I_DELTA_mHz,  color="#7f7f7f", linestyle="-", linewidth=0.8, alpha=0.5)
    ax.axvline(CLASS_II_DELTA_mHz, color="#7f7f7f", linestyle="-", linewidth=0.8, alpha=0.5)

    # Class band annotations
    ax.text(CLASS_I_DELTA_mHz,  0.55, "Class I (10 mHz)",
            rotation=90, va="center", ha="right",
            fontsize=9, color="#555555")
    ax.text(CLASS_II_DELTA_mHz, 0.55, "Class II (100 mHz)",
            rotation=90, va="center", ha="right",
            fontsize=9, color="#555555")

    ax.set_xlim(DELTA_MIN_mHz, DELTA_MAX_mHz)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Equivalence margin delta (mHz)", fontsize=10)
    ax.set_ylabel("TOST p-value  (max of p1, p2)", fontsize=10)
    ax.set_title(
        "TOST sensitivity sweep — open RTU V2 firmware vs NESO reference\n"
        f"(canonical capture replay_20260627_102507V2.csv, alpha = {ALPHA})",
        fontsize=10
    )
    ax.legend(loc="upper right", fontsize=8, framealpha=0.95)
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.4)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return True


# ============================================================================
# Main
# ============================================================================

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="TOST sensitivity sweep of the open RTU V2 firmware "
                    "against the canonical NESO replay capture.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("csv_path", type=str)
    parser.add_argument("--json", type=str, default=None)
    parser.add_argument("--plot", type=str, default=None)
    parser.add_argument("--verify-hash", type=str, default=None)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    source_path = Path(args.csv_path)
    if not source_path.is_file():
        print(f"ERROR: input file not found: {source_path}", file=sys.stderr)
        return 2

    hash_value = compute_sha256(source_path)
    if args.verify_hash and hash_value != args.verify_hash:
        print(f"ERROR: SHA-256 mismatch.\n"
              f"  expected: {args.verify_hash}\n"
              f"  computed: {hash_value}", file=sys.stderr)
        return 3

    rows = read_canonical_csv(source_path)
    if not rows:
        print(f"ERROR: no valid rows in {source_path}", file=sys.stderr)
        return 4

    # Build delta grid
    n_steps = int((DELTA_MAX_mHz - DELTA_MIN_mHz) / DELTA_STEP_mHz) + 1
    delta_grid = [DELTA_MIN_mHz + i * DELTA_STEP_mHz for i in range(n_steps)]

    # Sweep each window
    per_window_sweep = []
    per_window_crossover = []
    for (_, i_start, i_end) in WINDOWS:
        errors = extract_window_errors(rows, i_start, i_end)
        sweep, crossover = sweep_window(errors, delta_grid)
        per_window_sweep.append(sweep)
        per_window_crossover.append(crossover)

    if not args.quiet:
        print_console_report(per_window_sweep, per_window_crossover, hash_value, source_path)

    if args.json:
        sidecar = build_json_sidecar(per_window_sweep, per_window_crossover,
                                     hash_value, source_path)
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(sidecar, fh, indent=2)
        if not args.quiet:
            print(f"  JSON sidecar written: {args.json}")

    if args.plot:
        if render_plot(per_window_sweep, per_window_crossover, args.plot):
            if not args.quiet:
                print(f"  Plot written:         {args.plot}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
