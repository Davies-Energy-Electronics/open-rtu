#!/usr/bin/env python3
"""
tve_metrics.py — Total Vector Error (TVE) characterisation of the open RTU
canonical NESO replay capture.

Computes a TVE-equivalent figure of merit per frame and aggregates statistics
over the full 360 s replay window and over the five regions identified in
§5.3.4 (pre-event nominal, RoCoF descent, nadir, late recovery, post-event
settled).

Following IEEE C37.118.1-2011 §5.2, the per-frame TVE-equivalent is:

    TVE(n) = sqrt[ (V_meas(n) - V_ref)^2
                 + ( V_meas(n) * 2*pi * Delta_f(n) * T_report )^2 ]
             / V_ref

where V_ref = K_design * V_nom = 0.2364 * 3.25 = 0.7682 V (constant nominal
RMS), Delta_f(n) = f_meas(n) - f_NESO(n), and T_report = 0.020 s (the SCADA
50 Hz reporting cadence, matching the PMU C37.118 frame rate convention).

The formulation deliberately uses the per-frame Delta_f scaled by T_report
rather than the cumulative integrated phase drift, so the metric reflects
the per-reporting-window phase error contribution (the PMU-relevant
quantity) rather than the unbounded long-baseline drift that would arise
without GPS-locked timestamping.

PMU compliance threshold (IEEE C37.118.1): TVE <= 1.0%.

Usage:
    python tve_metrics.py replay_20260627_102507.csv
    python tve_metrics.py replay_20260627_102507.csv --json tve_output.json
    python tve_metrics.py replay_20260627_102507.csv --verify-hash dadbcfe247...

Output:
    - Console table with aggregate TVE and per-region TVE statistics
    - Optional JSON sidecar with full per-region breakdown

Reproducibility:
    Canonical input: replay_20260627_102507.csv
    Expected SHA-256: dadbcfe247dd4e538b71a891b2b26f070d02c42e3dd4622b5d2125aff169560f
    Expected aggregate TVE max: ~7.9% (at replay index 176, nadir region)
    Expected aggregate TVE mean: ~2.0%
    Expected pre-event TVE max: ~1.0% (below PMU threshold)
    Expected settled TVE max: ~1.1% (at PMU threshold boundary)

Standard-library only. No external dependencies.
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
# Constants (paper §2.3, §3.1, §5.3)
# ============================================================================

K_DESIGN = 0.2364                 # Conditioning circuit divider ratio (paper §2.3)
V_NOM = 3.25                      # Nominal pre-divider RMS (paper §1, §2.2)
V_REF = K_DESIGN * V_NOM          # Reference RMS at ADC input: 0.7682 V
T_REPORT = 0.020                  # SCADA reporting cadence in seconds (50 Hz)
PMU_TVE_THRESHOLD = 0.01          # IEEE C37.118.1 TVE compliance threshold: 1.0%

# Five regions of the 360 s replay window per §5.3.4
REGIONS = [
    ("Pre-event nominal",       0,    45),
    ("RoCoF descent",          46,   119),
    ("Nadir + early recovery", 120,  200),
    ("Late recovery",          201,  300),
    ("Post-event settled",     301,  359),
]

# NESO 9 August 2019 reference frequency trajectory embedded as fallback.
# In normal operation this is read from the CSV's `f_neso` column if present.
# The embedded sequence below is a sparse fallback for testing only.
NESO_FALLBACK_SEQUENCE_PRESENT = False  # set True only if embedded data provided

# ============================================================================
# Core computation
# ============================================================================

def compute_per_frame_tve(v_meas, f_meas, f_neso):
    """
    Per-frame TVE-equivalent.

    Inputs:
        v_meas: measured RMS magnitude (V)
        f_meas: measured frequency (Hz)
        f_neso: NESO reference frequency (Hz)

    Returns:
        Tuple (tve, mag_contrib, phase_contrib) all as dimensionless fractions:
            tve          — full TVE per IEEE C37.118.1
            mag_contrib  — magnitude-only error |V_meas − V_ref| / V_ref
            phase_contrib — phase-only error |V_meas · 2π · Δf · T_report| / V_ref
        Note that tve = sqrt(mag_contrib² + phase_contrib²) by construction.
    """
    delta_v = v_meas - V_REF
    delta_f = f_meas - f_neso
    phase_term = v_meas * 2.0 * math.pi * delta_f * T_REPORT
    mag_contrib = abs(delta_v) / V_REF
    phase_contrib = abs(phase_term) / V_REF
    tve = math.sqrt(delta_v * delta_v + phase_term * phase_term) / V_REF
    return tve, mag_contrib, phase_contrib


def aggregate_tve(tve_values):
    """
    Compute aggregate statistics on a TVE sequence.

    Returns dict with max, mean, median, p95, std, count.
    """
    if not tve_values:
        return {"max": 0.0, "mean": 0.0, "median": 0.0,
                "p95": 0.0, "std": 0.0, "count": 0}
    sorted_v = sorted(tve_values)
    p95_index = max(0, int(0.95 * len(sorted_v)) - 1)
    return {
        "max":    max(tve_values),
        "mean":   statistics.mean(tve_values),
        "median": statistics.median(tve_values),
        "p95":    sorted_v[p95_index],
        "std":    statistics.stdev(tve_values) if len(tve_values) >= 2 else 0.0,
        "count":  len(tve_values),
    }


def classify_compliance(tve_max):
    """
    Classify the TVE-equivalent against IEEE C37.118.1 PMU compliance.

    Returns string verdict: 'PMU-compliant' (<= 1.0%) or 'Above PMU threshold' (> 1.0%).
    """
    return "PMU-compliant" if tve_max <= PMU_TVE_THRESHOLD else "Above PMU threshold"


# ============================================================================
# CSV reader and SHA-256 integrity check
# ============================================================================

def compute_sha256(path):
    """Compute SHA-256 of the input file for reproducibility logging."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def read_canonical_csv(path):
    """
    Read the canonical replay CSV. Returns list of dicts with keys:
        replay_index (int), v_meas (float), f_meas (float), f_neso (float)

    Expected CSV columns (case-insensitive matching for robustness):
        - replay_index OR i  (the NESO-trajectory second index, 0..359)
        - Vb_rms OR V_b_rms OR v_b_rms  (measured RMS magnitude on the canonical Vb channel)
        - f_Vb OR f_b OR f_vb  (measured frequency on Vb)
        - f_neso OR f_target OR f_ref  (NESO reference frequency at this index)
    """
    rows = []
    with open(path, "r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        # Normalise header to lower case for matching
        field_map = {f.lower(): f for f in reader.fieldnames}

        def find_field(candidates):
            for c in candidates:
                if c.lower() in field_map:
                    return field_map[c.lower()]
            return None

        idx_col = find_field(["replay_index", "i", "index"])
        v_col   = find_field(["Vb_rms", "V_b_rms", "v_b_rms", "vb_rms"])
        f_col   = find_field(["f_Vb", "f_b", "f_vb"])
        ref_col = find_field(["f_neso", "f_target_hz", "f_target", "f_ref"])

        missing = []
        if not idx_col: missing.append("replay_index")
        if not v_col:   missing.append("Vb_rms")
        if not f_col:   missing.append("f_Vb")
        if not ref_col: missing.append("f_neso")
        if missing:
            raise ValueError(
                f"Required column(s) not found in CSV: {', '.join(missing)}.\n"
                f"Available columns: {', '.join(reader.fieldnames)}"
            )

        for row in reader:
            try:
                rows.append({
                    "replay_index": int(float(row[idx_col])),
                    "v_meas":       float(row[v_col]),
                    "f_meas":       float(row[f_col]),
                    "f_neso":       float(row[ref_col]),
                })
            except (ValueError, KeyError):
                continue  # Skip malformed rows silently
    return rows


# ============================================================================
# Region analysis
# ============================================================================

def analyse_regions(rows):
    """
    Compute per-region TVE statistics for the five §5.3.4 regions.

    Inputs: list of dicts from read_canonical_csv()
    Returns: dict mapping region name to its stats dict, with separate
             magnitude-contrib and phase-contrib breakdowns.
    """
    # Aggregate by replay_index: if multiple frames share an index, use mean
    by_index = {}
    for r in rows:
        idx = r["replay_index"]
        by_index.setdefault(idx, []).append(r)

    region_stats = {}
    for region_name, idx_start, idx_end in REGIONS:
        tve_values = []
        mag_values = []
        phase_values = []
        for idx in range(idx_start, idx_end + 1):
            if idx not in by_index:
                continue
            for r in by_index[idx]:
                tve, mag, ph = compute_per_frame_tve(
                    r["v_meas"], r["f_meas"], r["f_neso"]
                )
                tve_values.append(tve)
                mag_values.append(mag)
                phase_values.append(ph)
        stats = aggregate_tve(tve_values)
        stats["mag_stats"] = aggregate_tve(mag_values)
        stats["phase_stats"] = aggregate_tve(phase_values)
        stats["compliance"] = classify_compliance(stats["max"])
        stats["index_range"] = f"{idx_start}–{idx_end}"
        region_stats[region_name] = stats

    return region_stats


# ============================================================================
# Output formatting
# ============================================================================

def format_pct(v):
    """Format a TVE value as percentage with 2 decimal places."""
    return f"{100*v:.2f}%"


def print_console_report(rows, region_stats, hash_value, source_path):
    """Print a human-readable summary report to stdout."""
    all_tve = []
    all_mag = []
    all_phase = []
    for r in rows:
        tve, mag, ph = compute_per_frame_tve(r["v_meas"], r["f_meas"], r["f_neso"])
        all_tve.append(tve)
        all_mag.append(mag)
        all_phase.append(ph)
    overall = aggregate_tve(all_tve)
    overall_mag = aggregate_tve(all_mag)
    overall_phase = aggregate_tve(all_phase)

    print("=" * 78)
    print("  TVE-EQUIVALENT CHARACTERISATION — OPEN RTU NESO REPLAY VALIDATION")
    print("=" * 78)
    print(f"  Source CSV:           {source_path}")
    print(f"  SHA-256 (input):      {hash_value}")
    print(f"  Frames processed:     {len(rows)}")
    print(f"  V_ref (nominal):      {V_REF:.4f} V  (K_design = {K_DESIGN}, V_nom = {V_NOM} V)")
    print(f"  T_report (PMU 50Hz):  {T_REPORT*1000:.0f} ms")
    print(f"  PMU TVE threshold:    {100*PMU_TVE_THRESHOLD:.1f}% (IEEE C37.118.1)")
    print()
    print("  ── PER-REGION RESULTS (§5.3.4 five-region partition) ─────────────")
    print()
    print(f"  {'Region':<26} {'count':>6} {'TVE max':>8} {'TVE μ':>8} "
          f"{'|mag| μ':>8} {'|φ| μ':>8} {'verdict':>22}")
    print(f"  {'-'*26} {'-'*6} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*22}")
    for region_name, _, _ in REGIONS:
        s = region_stats[region_name]
        print(f"  {region_name:<26} {s['count']:>6} "
              f"{format_pct(s['max']):>8} {format_pct(s['mean']):>8} "
              f"{format_pct(s['mag_stats']['mean']):>8} {format_pct(s['phase_stats']['mean']):>8} "
              f"{s['compliance']:>22}")
    print()
    print("  ── FULL 360 s AGGREGATE ───────────────────────────────────────────")
    print()
    print(f"  TVE max:                {format_pct(overall['max'])}")
    print(f"  TVE mean:               {format_pct(overall['mean'])}")
    print(f"  TVE median:             {format_pct(overall['median'])}")
    print(f"  TVE 95th percentile:    {format_pct(overall['p95'])}")
    print(f"  TVE σ:                  {format_pct(overall['std'])}")
    print()
    print(f"  Magnitude contrib mean: {format_pct(overall_mag['mean'])}  "
          f"(max {format_pct(overall_mag['max'])})")
    print(f"  Phase contrib mean:     {format_pct(overall_phase['mean'])}  "
          f"(max {format_pct(overall_phase['max'])})")
    print()
    print(f"  Verdict: {classify_compliance(overall['max'])}")
    print("=" * 78)


def build_json_sidecar(rows, region_stats, hash_value, source_path):
    """Build the JSON sidecar dict for downstream tooling and the paper."""
    all_tve = []
    all_mag = []
    all_phase = []
    for r in rows:
        tve, mag, ph = compute_per_frame_tve(r["v_meas"], r["f_meas"], r["f_neso"])
        all_tve.append(tve)
        all_mag.append(mag)
        all_phase.append(ph)
    overall = aggregate_tve(all_tve)
    overall_mag = aggregate_tve(all_mag)
    overall_phase = aggregate_tve(all_phase)
    overall["compliance"] = classify_compliance(overall["max"])

    return {
        "schema_version": "1.1",
        "tool": "tve_metrics.py",
        "input": {
            "source_file": str(source_path),
            "sha256": hash_value,
            "frames_processed": len(rows),
        },
        "constants": {
            "K_design": K_DESIGN,
            "V_nom_V": V_NOM,
            "V_ref_V": V_REF,
            "T_report_s": T_REPORT,
            "pmu_tve_threshold_pct": 100 * PMU_TVE_THRESHOLD,
        },
        "overall": {
            "tve_max_pct":    100 * overall["max"],
            "tve_mean_pct":   100 * overall["mean"],
            "tve_median_pct": 100 * overall["median"],
            "tve_p95_pct":    100 * overall["p95"],
            "tve_std_pct":    100 * overall["std"],
            "magnitude_contrib_mean_pct": 100 * overall_mag["mean"],
            "magnitude_contrib_max_pct":  100 * overall_mag["max"],
            "phase_contrib_mean_pct":     100 * overall_phase["mean"],
            "phase_contrib_max_pct":      100 * overall_phase["max"],
            "compliance":     overall["compliance"],
        },
        "regions": {
            region_name: {
                "index_range":                region_stats[region_name]["index_range"],
                "frame_count":                region_stats[region_name]["count"],
                "tve_max_pct":                100 * region_stats[region_name]["max"],
                "tve_mean_pct":               100 * region_stats[region_name]["mean"],
                "tve_p95_pct":                100 * region_stats[region_name]["p95"],
                "tve_std_pct":                100 * region_stats[region_name]["std"],
                "magnitude_contrib_mean_pct": 100 * region_stats[region_name]["mag_stats"]["mean"],
                "magnitude_contrib_max_pct":  100 * region_stats[region_name]["mag_stats"]["max"],
                "phase_contrib_mean_pct":     100 * region_stats[region_name]["phase_stats"]["mean"],
                "phase_contrib_max_pct":      100 * region_stats[region_name]["phase_stats"]["max"],
                "compliance":                 region_stats[region_name]["compliance"],
            }
            for region_name, _, _ in REGIONS
        },
    }


# ============================================================================
# Command-line entry point
# ============================================================================

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Compute TVE-equivalent characterisation of the open RTU "
                    "canonical NESO replay capture (per IEEE C37.118.1 §5.2).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("csv_path", type=str,
                        help="Path to the canonical replay CSV "
                             "(e.g. replay_20260627_102507.csv)")
    parser.add_argument("--json", type=str, default=None,
                        help="Path to write JSON sidecar output "
                             "(default: stdout console table only)")
    parser.add_argument("--verify-hash", type=str, default=None,
                        help="Expected SHA-256 of the input CSV; abort if mismatch.")
    parser.add_argument("--quiet", action="store_true",
                        help="Suppress console table output (JSON only).")
    args = parser.parse_args(argv)

    source_path = Path(args.csv_path)
    if not source_path.is_file():
        print(f"ERROR: input file not found: {source_path}", file=sys.stderr)
        return 2

    hash_value = compute_sha256(source_path)
    if args.verify_hash and hash_value != args.verify_hash:
        print(f"ERROR: SHA-256 mismatch. Expected {args.verify_hash}, "
              f"got {hash_value}", file=sys.stderr)
        return 3

    rows = read_canonical_csv(source_path)
    if not rows:
        print(f"ERROR: no valid rows found in {source_path}", file=sys.stderr)
        return 4

    region_stats = analyse_regions(rows)

    if not args.quiet:
        print_console_report(rows, region_stats, hash_value, source_path)

    if args.json:
        sidecar = build_json_sidecar(rows, region_stats, hash_value, source_path)
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(sidecar, fh, indent=2)
        if not args.quiet:
            print(f"\n  JSON sidecar written: {args.json}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
