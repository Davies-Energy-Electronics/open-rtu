"""
runtime_metrics.py
==================

Processes the runtime-characterisation CSV produced by the timing-instrumented
variant of the canonical V2 firmware (Feather1codeMimoDemo_TIMED.ino) and
reports per-module execution-time statistics matching the row schema of
Table 13 in Section 5.4.4 of the open-RTU MDPI Sensors paper.

For each of the eight measured stages (M1 RMS, M2 Frequency, M5 THD,
M4-phase, M4-pwr, alarm evaluation, M6 DNP3 build, Serial output) the script
reports:

    n          — frames measured after warm-up
    mean_us    — sample-mean execution time in microseconds
    std_us     — sample standard deviation
    p99_us     — 99th percentile (the rare-worst-case figure)
    max_us     — sample maximum
    pct_cycle  — mean execution time as a percentage of the measured cycle

The total-cycle row (sum across all eight stages) is computed as a µs/ms
figure and compared against three reference budgets:

    - PMU-class target:        20 ms (IEEE C37.118.2 50 Hz frame rate)
    - IEC 61400-21 Class II:  100 ms (intermediate target)
    - Current V2 architecture: ~ measured (the actual cycle time)

A headline "engineering quantification" line for §6.3 Gate 2 is produced
that breaks the total cycle into acquisition-paced (delayMicroseconds)
versus computation, since the §6.3 Gate 2 argument is that continuous DMA
into a ring buffer removes the acquisition pacing as a cycle-time
constraint.

Usage:
    python runtime_metrics.py timing_capture.csv
    python runtime_metrics.py timing_capture.csv --json timing_capture_runtime.json

Input format (one row per main-loop iteration, CSV; lines beginning with '#'
or 'TIMING_CSV_HEADER:' are treated as comments and skipped):
    iter,M1_RMS_us,M2_freq_us,M5_THD_us,M4_phase_us,M4_pwr_us,alarm_us,M6_DNP3_us,serial_us

Standard library only. No external dependencies.

Reproducibility:
    Expected input  : timing_capture.csv (1000 rows after firmware warm-up)
    Expected output : timing_capture_runtime.json
    PMU target      : 20 000 us (PMU 50 Hz frame convention)

Author: Jack Davies
"""

from __future__ import annotations
import argparse, csv, json, math, os, statistics, sys


# ============================================================================
# Constants
# ============================================================================

PMU_TARGET_US      =  20_000       # PMU-class 50 Hz frame rate
CLASS_II_TARGET_US = 100_000       # IEC 61400-21 Class II intermediate

MODULE_LABELS = [
    ("M1_RMS_us",    "M1  — RMS (5 channels)"),
    ("M2_freq_us",   "M2  — Frequency (3 channels)"),
    ("M5_THD_us",    "M5  — THD / Harmonics (3 channels)"),
    ("M4_phase_us",  "M4  — Phase angle (2 channels)"),
    ("M4_pwr_us",    "M4  — Power + P/Q/S/PF (2 channels)"),
    ("alarm_us",     "Alarm evaluation"),
    ("M6_DNP3_us",   "M6  — DNP3 frame build"),
    ("serial_us",    "Serial output (verbose debug)"),
]
COLUMN_KEYS = [k for k, _ in MODULE_LABELS]


# ============================================================================
# Statistical helpers
# ============================================================================

def percentile(sorted_values, p):
    if not sorted_values:
        return float("nan")
    if len(sorted_values) == 1:
        return float(sorted_values[0])
    rank = (p / 100.0) * (len(sorted_values) - 1)
    lo = int(math.floor(rank))
    hi = int(math.ceil(rank))
    if lo == hi:
        return float(sorted_values[lo])
    frac = rank - lo
    return float(sorted_values[lo]) * (1 - frac) + float(sorted_values[hi]) * frac


def summarise(values):
    n = len(values)
    if n == 0:
        return {"n": 0, "mean_us": float("nan"), "std_us": float("nan"),
                "p99_us": float("nan"), "max_us": float("nan")}
    mean = statistics.mean(values)
    std = statistics.stdev(values) if n > 1 else 0.0
    sorted_v = sorted(values)
    return {
        "n": n,
        "mean_us": mean,
        "std_us":  std,
        "p99_us":  percentile(sorted_v, 99.0),
        "max_us":  float(sorted_v[-1]),
    }


# ============================================================================
# CSV loader — handles TIMING-prefixed lines and pure CSV alike
# ============================================================================

def load_timing_csv(path):
    columns = {key: [] for key in COLUMN_KEYS}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            # Accept both 'TIMING,iter,...' (from raw serial capture) and
            # 'iter,...' (from cleaned CSV)
            if line.startswith("TIMING,"):
                line = line[len("TIMING,"):]
            # Also skip the header line if present
            if "iter," in line or "M1_RMS" in line:
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 9:
                continue
            try:
                module_values = [int(parts[1 + i]) for i in range(8)]
            except (ValueError, IndexError):
                continue
            for key, val in zip(COLUMN_KEYS, module_values):
                columns[key].append(val)
    return columns


# ============================================================================
# Acquisition-pacing estimate
# ============================================================================

def build_acq_estimate(module_stats):
    """Per-module acquisition-pacing budget in microseconds.

    For the fixed-window stages, acquisition time is the analytical product of
    (channels x samples-per-window x per-sample delayMicroseconds); these are
    the delayMicroseconds() calls that pace the ADC reads and dominate each
    measure*() function.

    The M4-phase stage is the exception: its sampling loop early-exits on the
    detected voltage zero-crossing rather than always running the full
    PHASE_WIN=200-sample window, so the analytical worst case (2 x 200 x 460us
    = 184 ms) does not represent the acquisition time actually incurred. Using
    that worst case makes the summed acquisition estimate exceed the measured
    cycle time and yields an impossible negative "residual compute" figure.
    We therefore use the MEASURED per-frame mean for M4-phase, which correctly
    reflects the early-exit behaviour (measured ~25 ms, not 184 ms). The
    firmware constants remain the source for every fixed-window stage, so the
    estimate stays analytical wherever the analytical model is valid.
    """
    return {
        "M1_RMS_us":   5 * 40  * 500,    # 5 channels x RMS_N x 500us
        "M2_freq_us":  3 * 200 * 500,    # 3 channels x FREQ_WIN x 500us
        "M5_THD_us":   3 * 200 * 470,    # 3 channels x THD_N x 470us
        # M4-phase early-exits on zero-crossing: use the measured mean, not the
        # analytical 2 x 200 x 460us worst case (see docstring).
        "M4_phase_us": module_stats["M4_phase_us"]["mean_us"],
        "M4_pwr_us":   2 * 100 * 460,    # 2 channels x PWR_N x 460us
        "alarm_us":    0,
        "M6_DNP3_us":  0,
        "serial_us":   0,
    }


# ============================================================================
# Reporting
# ============================================================================

def build_table_13_rows(module_stats, cycle_stats):
    rows = []
    cycle_mean = cycle_stats["mean_us"]
    for key, label in MODULE_LABELS:
        s = module_stats[key]
        rows.append({
            "module":      label,
            "n":           s["n"],
            "mean_us":     round(s["mean_us"], 1),
            "mean_ms":     round(s["mean_us"] / 1000.0, 2),
            "p99_us":      round(s["p99_us"],  1),
            "p99_ms":      round(s["p99_us"] / 1000.0, 2),
            "max_us":      round(s["max_us"],  1),
            "pct_cycle":   (round(100.0 * s["mean_us"] / cycle_mean, 2)
                            if cycle_mean > 0 else 0.0),
        })
    rows.append({
        "module":      "TOTAL (per-cycle, summed mean)",
        "n":           cycle_stats["n"],
        "mean_us":     round(cycle_stats["mean_us"], 1),
        "mean_ms":     round(cycle_stats["mean_us"] / 1000.0, 2),
        "p99_us":      round(cycle_stats["p99_us"], 1),
        "p99_ms":      round(cycle_stats["p99_us"] / 1000.0, 2),
        "max_us":      round(cycle_stats["max_us"], 1),
        "pct_cycle":   100.00,
    })
    return rows


def print_console_report(rows, cycle_stats, module_stats, source_path):
    cycle_mean_ms = cycle_stats["mean_us"] / 1000.0
    cycle_p99_ms  = cycle_stats["p99_us"]  / 1000.0

    print("=" * 88)
    print("  RUNTIME / COMPUTATIONAL-COST CHARACTERISATION")
    print("  Canonical V2 firmware, Feather 1, ESP32 @ 240 MHz")
    print("=" * 88)
    print(f"  Source CSV:           {os.path.basename(source_path)}")
    print(f"  Iterations measured:  {rows[0]['n']}")
    print(f"  PMU-class target:     {PMU_TARGET_US} us (20 ms, IEEE C37.118.2)")
    print(f"  Class II target:      {CLASS_II_TARGET_US} us (100 ms, IEC 61400-21)")
    print(f"  Measured cycle mean:  {cycle_stats['mean_us']:.0f} us "
          f"({cycle_mean_ms:.1f} ms)")
    print(f"  Measured cycle p99:   {cycle_stats['p99_us']:.0f} us "
          f"({cycle_p99_ms:.1f} ms)")
    print()
    print(f"  {'Stage':<38} {'mean(ms)':>10} {'p99(ms)':>10} {'max(ms)':>10} {'% cycle':>10}")
    print(f"  {'-' * 38} {'-' * 10} {'-' * 10} {'-' * 10} {'-' * 10}")
    for r in rows[:-1]:
        print(f"  {r['module']:<38} {r['mean_ms']:>10.2f} "
              f"{r['p99_ms']:>10.2f} {r['max_us']/1000.0:>10.2f} {r['pct_cycle']:>9.2f}%")
    print(f"  {'-' * 38} {'-' * 10} {'-' * 10} {'-' * 10} {'-' * 10}")
    tot = rows[-1]
    print(f"  {tot['module']:<38} {tot['mean_ms']:>10.2f} "
          f"{tot['p99_ms']:>10.2f} {tot['max_us']/1000.0:>10.2f} {tot['pct_cycle']:>9.2f}%")
    print()

    # Headline for §6.3 Gate 2: how far above PMU-class target
    target_ratio = cycle_stats["mean_us"] / PMU_TARGET_US
    classII_ratio = cycle_stats["mean_us"] / CLASS_II_TARGET_US
    print(f"  Cycle vs PMU target (20 ms):       {target_ratio:.1f}x over target")
    print(f"  Cycle vs Class II target (100 ms): {classII_ratio:.1f}x over target")
    print()

    # Acquisition vs computation decomposition (analytical estimate)
    # Each measure*() function is dominated by delayMicroseconds() — the
    # acquisition pacing. We estimate acquisition-paced µs per call from
    # the firmware constants and report computation as the residual.
    acq_estimate_us = build_acq_estimate(module_stats)
    total_acquisition_us = sum(acq_estimate_us.values())
    total_compute_us = cycle_stats["mean_us"] - total_acquisition_us

    print(f"  Acquisition-pacing budget (analytical sum of delayMicroseconds):")
    print(f"    Estimated total:    {total_acquisition_us / 1000.0:.1f} ms")
    print(f"    Fraction of cycle:  {100.0 * total_acquisition_us / cycle_stats['mean_us']:.1f}%")
    print(f"  Residual (compute + Serial + scheduling overhead):")
    print(f"    {total_compute_us / 1000.0:.1f} ms "
          f"({100.0 * total_compute_us / cycle_stats['mean_us']:.1f}% of cycle)")
    print()

    # Engineering quantification for §6.3 Gate 2
    print("  Headline (engineering quantification for §6.3 Gate 2):")
    print(f"    Current cycle:                       {cycle_mean_ms:.0f} ms")
    print(f"    PMU-class target (C37.118.2 50 Hz):  20 ms")
    print(f"    Acquisition-pacing dominates ({100.0 * total_acquisition_us / cycle_stats['mean_us']:.0f}%);")
    print(f"    continuous DMA into a ring buffer would remove this constraint,")
    print(f"    leaving only ~{total_compute_us / 1000.0:.0f} ms of computation on the main thread.")
    print("=" * 88)


def build_json_sidecar(rows, cycle_stats, module_stats, source_path):
    # Recompute acquisition decomposition for JSON output
    acq_estimate_us = build_acq_estimate(module_stats)
    total_acquisition_us = sum(acq_estimate_us.values())
    total_compute_us = cycle_stats["mean_us"] - total_acquisition_us

    return {
        "schema_version": "1.0",
        "tool": "runtime_metrics.py",
        "source_file": os.path.basename(source_path),
        "n_iterations": rows[0]["n"],
        "targets": {
            "pmu_target_us": PMU_TARGET_US,
            "class_ii_target_us": CLASS_II_TARGET_US,
        },
        "cycle_stats": {
            "n":         cycle_stats["n"],
            "mean_us":   cycle_stats["mean_us"],
            "mean_ms":   cycle_stats["mean_us"] / 1000.0,
            "p99_us":    cycle_stats["p99_us"],
            "p99_ms":    cycle_stats["p99_us"] / 1000.0,
            "max_us":    cycle_stats["max_us"],
            "max_ms":    cycle_stats["max_us"] / 1000.0,
        },
        "table_13_rows": rows,
        "per_module_full_stats": {key: module_stats[key] for key in COLUMN_KEYS},
        "headline": {
            "cycle_mean_ms": cycle_stats["mean_us"] / 1000.0,
            "pmu_target_ratio": cycle_stats["mean_us"] / PMU_TARGET_US,
            "classII_target_ratio": cycle_stats["mean_us"] / CLASS_II_TARGET_US,
            "acquisition_estimate_us": total_acquisition_us,
            "acquisition_pct_of_cycle":
                100.0 * total_acquisition_us / cycle_stats["mean_us"]
                if cycle_stats["mean_us"] > 0 else None,
            "residual_compute_us": total_compute_us,
            "residual_pct_of_cycle":
                100.0 * total_compute_us / cycle_stats["mean_us"]
                if cycle_stats["mean_us"] > 0 else None,
        },
    }


# ============================================================================
# Main
# ============================================================================

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("csv_path", type=str)
    parser.add_argument("--json", type=str, default=None,
                        help="Output JSON sidecar path (default: <csv>_runtime.json)")
    args = parser.parse_args(argv)

    if not os.path.isfile(args.csv_path):
        print(f"ERROR: input file not found: {args.csv_path}", file=sys.stderr)
        return 2

    columns = load_timing_csv(args.csv_path)
    if not columns["M1_RMS_us"]:
        print(f"ERROR: no usable rows in {args.csv_path}", file=sys.stderr)
        print("       Expected lines beginning with 'TIMING,' or '<iter>,<M1>,...'",
              file=sys.stderr)
        return 3

    module_stats = {key: summarise(columns[key]) for key in COLUMN_KEYS}

    n_iter = len(columns["M1_RMS_us"])
    per_iter_totals = []
    for i in range(n_iter):
        per_iter_totals.append(sum(columns[k][i] for k in COLUMN_KEYS))
    cycle_stats = summarise(per_iter_totals)

    rows = build_table_13_rows(module_stats, cycle_stats)
    print_console_report(rows, cycle_stats, module_stats, args.csv_path)

    out_path = (args.json if args.json
                else os.path.splitext(args.csv_path)[0] + "_runtime.json")
    sidecar = build_json_sidecar(rows, cycle_stats, module_stats, args.csv_path)
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(sidecar, fh, indent=2)
    print(f"  JSON sidecar written: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
