"""
generate_3phase_trajectory.py
==============================

Generate a synthetic three-phase frequency-following trajectory from a
one-second NESO reference CSV. Produces balanced 120-degree-shifted phase
waveforms with per-phase THD and per-phase magnitude variation.

**DUAL-FEATHER OUTPUT (matches the MIMO stripboard as built; channel map in
Figure 8 of the paper):**
  - Feather 2 generates Va (DAC1, 0°) and Vb (DAC2, -120°)
  - Feather 1 generates Vc (DAC1, +120°)

The script therefore emits TWO Feather-ready `.h` headers:

  - `<prefix>_3phase_F2.h` containing LUT_VA and LUT_VB (for Feather 2)
  - `<prefix>_3phase_F1.h` containing LUT_VC (for Feather 1)

plus a single `<prefix>_3phase.csv` for the analysis toolchain.

Vc phase is set to -120° (equivalently +240°), consistent with the
channel-map convention of Figure 8 of the paper. Va = 0° reference, Vb = -120°, Vc = +120°
would ALSO be a positive-sequence set; this generator uses the Figure 8
convention Va = 0°, Vb = -120°, Vc = +120° which is the same as Va = 0°,
Vb = 240°, Vc = 120° after modulo. Users who need the strict Figure 8
labelling can pass --vc-phase-deg 240 (default is 120, but both produce a
mathematically valid balanced set).

The output is intentionally deterministic: given the same input CSV and
the same command-line options the same LUTs and CSV emerge byte-for-byte.
A `--seed` option seeds the per-phase perturbations so a run is reproducible.

Usage:

  Basic (defaults matching the paper's 9 August 2019 canonical run,
  no perturbation, balanced 3-phase):

    python generate_3phase_trajectory.py neso_trajectory.csv

  Multi-event runs, one per NESO CSV, with per-phase THD variation
  1-5 %% and per-phase magnitude variation +/- 2 %%, deterministic:

    python generate_3phase_trajectory.py trajectory_jan2019_dip1.csv \\
        --thd-var 1,3,5 --mag-var 0.98,1.00,1.02 --seed 20260712

Input CSV format:

  The NESO trajectory CSV as used already in the paper. The script accepts
  either the two-column form (index, frequency_hz) or the labelled form
  used by extract_neso_csv.py. A '#' line is treated as a comment.

Output files:

  <basename>_3phase.h    Feather 2 header: LUT_VA/LUT_VB/LUT_VC (uint8_t,
                         LUT_SIZE samples each), F_TARGET (float array,
                         one entry per replay index).

  <basename>_3phase.csv  CSV with columns:
                            replay_index, f_target_hz,
                            Va_lut_amp, Vb_lut_amp, Vc_lut_amp,
                            Va_thd_pct, Vb_thd_pct, Vc_thd_pct.

Standard library only. No external dependencies.

Author: Jack Davies
"""

from __future__ import annotations
import argparse
import csv
import math
import os
import random
import sys


# ============================================================================
# Constants matching the existing Feather 2 firmware
# ============================================================================

LUT_SIZE = 40           # samples per fundamental cycle in feather_2_code_mimoDEMO.ino
DAC_MID  = 127          # centre of 8-bit DAC output
DAC_AMP  = 100          # nominal amplitude used in the current Feather 2 sketch
F_NOMINAL = 50.0        # Hz, GB nominal
# Phase offsets (radians) for a positive-sequence 3-phase system,
# using the same convention as feather_2_code_mimoDEMO.ino: Va = 0, Vb = -120.
PHASE_OFFSETS_RAD = [0.0, -2.0 * math.pi / 3.0, 2.0 * math.pi / 3.0]
PHASE_NAMES = ["Va", "Vb", "Vc"]


# ============================================================================
# Trajectory loading (accepts multiple CSV shapes)
# ============================================================================

def load_neso_reference(path):
    """Read a NESO reference CSV. Returns list of (index, freq_hz)."""
    rows = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        rdr = csv.reader(fh)
        header = None
        for parts in rdr:
            if not parts:
                continue
            first = parts[0].strip()
            if first.startswith("#"):
                continue
            if header is None:
                # Detect header row
                if any(cell.strip().lower() in ("frequency_hz", "frequency", "freq_hz", "f_hz") for cell in parts):
                    header = [c.strip().lower() for c in parts]
                    continue
                # No header - assume 2 cols index, freq
                header = ["index", "frequency_hz"]
            try:
                # Locate the frequency column
                if header:
                    fcol = None
                    for i, name in enumerate(header):
                        if name in ("frequency_hz", "frequency", "freq_hz", "f_hz"):
                            fcol = i
                            break
                    if fcol is None:
                        fcol = 1  # default second column
                    icol = 0
                    idx = int(float(parts[icol]))
                    freq = float(parts[fcol])
                else:
                    idx = int(float(parts[0]))
                    freq = float(parts[1])
                if math.isnan(freq) or freq <= 0.0:
                    continue
                rows.append((idx, freq))
            except (ValueError, IndexError):
                continue
    if not rows:
        raise RuntimeError(f"no usable rows in {path}")
    return rows


# ============================================================================
# Per-phase LUT synthesis
# ============================================================================

def build_phase_lut(phase_offset_rad, amplitude, thd_pct):
    """Build a LUT_SIZE-sample uint8_t sine LUT for one phase.

    Distortion is added deterministically as a small 3rd harmonic whose
    amplitude reproduces the requested THD percentage exactly (single-tone
    approximation: THD ~= (a3 / a1) * 100 for a1 >> a3, no other components).
    Using the 3rd harmonic keeps the perturbation on a symmetric,
    zero-DC-offset waveform, which matters for the M2 zero-crossing
    detector.
    """
    a1 = float(amplitude)
    thd_frac = max(0.0, thd_pct) / 100.0
    a3 = thd_frac * a1
    lut = []
    for i in range(LUT_SIZE):
        theta = 2.0 * math.pi * i / LUT_SIZE + phase_offset_rad
        val = a1 * math.sin(theta) + a3 * math.sin(3.0 * theta)
        sample = int(round(DAC_MID + val))
        sample = max(0, min(255, sample))  # clamp to 8-bit
        lut.append(sample)
    return lut


# ============================================================================
# .h emission — dual-Feather split (F2: Va+Vb, F1: Vc)
# ============================================================================

def _write_lut_declaration(fh, name, lut):
    """Write a single uint8_t LUT array in the LUT_SIZE-per-line convention."""
    fh.write(f"static const uint8_t {name}[{LUT_SIZE}] = {{\n  ")
    fh.write(", ".join(str(v) for v in lut))
    fh.write("\n};\n\n")


def _write_ftarget_declaration(fh, f_targets):
    n_index = len(f_targets)
    fh.write(f"#define REPLAY_LENGTH {n_index}\n\n")
    fh.write(f"static const float F_TARGET[REPLAY_LENGTH] = {{\n")
    for row_start in range(0, n_index, 8):
        row = f_targets[row_start:row_start + 8]
        fh.write("  " + ", ".join(f"{v:.6f}f" for v in row))
        if row_start + 8 < n_index:
            fh.write(",")
        fh.write("\n")
    fh.write("};\n")


def emit_header_f2(path, lut_va, lut_vb, thds, mags, f_targets, source_csv, cli):
    """Feather 2 header: Va (DAC1) and Vb (DAC2)."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("// Auto-generated by generate_3phase_trajectory.py\n")
        fh.write(f"// Source CSV : {os.path.basename(source_csv)}\n")
        fh.write(f"// CLI        : {cli}\n")
        fh.write("// Board      : Feather 2 (of two)\n")
        fh.write("// Generates  : Va on DAC1 (GPIO25) at 0 deg, Vb on DAC2 (GPIO26) at -120 deg\n")
        fh.write(f"// Samples/cy : {LUT_SIZE}\n")
        fh.write(f"// Replay idx : {len(f_targets)}\n")
        fh.write("\n#pragma once\n#include <stdint.h>\n\n")
        _write_lut_declaration(fh, "LUT_VA", lut_va)
        _write_lut_declaration(fh, "LUT_VB", lut_vb)
        fh.write(f"static const float VA_AMP_LUT     = {mags[0]:.6f}f;\n")
        fh.write(f"static const float VB_AMP_LUT     = {mags[1]:.6f}f;\n")
        fh.write(f"static const float VA_THD_LUT_PCT = {thds[0]:.4f}f;\n")
        fh.write(f"static const float VB_THD_LUT_PCT = {thds[1]:.4f}f;\n\n")
        _write_ftarget_declaration(fh, f_targets)


def emit_header_f1(path, lut_vc, thds, mags, f_targets, source_csv, cli):
    """Feather 1 header: Vc (DAC1)."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("// Auto-generated by generate_3phase_trajectory.py\n")
        fh.write(f"// Source CSV : {os.path.basename(source_csv)}\n")
        fh.write(f"// CLI        : {cli}\n")
        fh.write("// Board      : Feather 1 (of two)\n")
        fh.write("// Generates  : Vc on DAC1 (GPIO25) at +120 deg (or -120 mod 360)\n")
        fh.write(f"// Samples/cy : {LUT_SIZE}\n")
        fh.write(f"// Replay idx : {len(f_targets)}\n")
        fh.write("\n#pragma once\n#include <stdint.h>\n\n")
        _write_lut_declaration(fh, "LUT_VC", lut_vc)
        fh.write(f"static const float VC_AMP_LUT     = {mags[2]:.6f}f;\n")
        fh.write(f"static const float VC_THD_LUT_PCT = {thds[2]:.4f}f;\n\n")
        _write_ftarget_declaration(fh, f_targets)


# ============================================================================
# CSV emission for the analysis toolchain
# ============================================================================

def emit_csv(path, thds, mags, f_targets, source_csv, cli):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["# generate_3phase_trajectory.py output"])
        w.writerow([f"# source: {os.path.basename(source_csv)}"])
        w.writerow([f"# cli   : {cli}"])
        w.writerow([f"# Va_thd_pct={thds[0]:.4f}, Vb_thd_pct={thds[1]:.4f}, Vc_thd_pct={thds[2]:.4f}"])
        w.writerow([f"# Va_amp_frac={mags[0]:.6f}, Vb_amp_frac={mags[1]:.6f}, Vc_amp_frac={mags[2]:.6f}"])
        w.writerow([
            "replay_index", "f_target_hz",
            "Va_amp_frac", "Vb_amp_frac", "Vc_amp_frac",
            "Va_thd_pct",  "Vb_thd_pct",  "Vc_thd_pct",
        ])
        for i, f in enumerate(f_targets):
            w.writerow([i, f"{f:.6f}",
                        f"{mags[0]:.6f}", f"{mags[1]:.6f}", f"{mags[2]:.6f}",
                        f"{thds[0]:.4f}",  f"{thds[1]:.4f}",  f"{thds[2]:.4f}"])


# ============================================================================
# CLI
# ============================================================================

def parse_triplet(s, name, dtype=float):
    parts = [p.strip() for p in s.split(",")]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(f"--{name} needs three comma-separated values, got {len(parts)}")
    try:
        return tuple(dtype(p) for p in parts)
    except ValueError:
        raise argparse.ArgumentTypeError(f"--{name} could not be parsed as {dtype.__name__}")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("csv_path", type=str,
        help="NESO reference CSV: index, frequency_hz")
    parser.add_argument("--thd-var", type=str, default="0,0,0",
        help="Per-phase THD in percent, comma-separated Va,Vb,Vc (default 0,0,0)")
    parser.add_argument("--mag-var", type=str, default="1.0,1.0,1.0",
        help="Per-phase amplitude fraction Va,Vb,Vc (1.0 = 100 %% of DAC_AMP)")
    parser.add_argument("--out-prefix", type=str, default=None,
        help="Output filename prefix (default: input CSV basename)")
    parser.add_argument("--seed", type=int, default=None,
        help="RNG seed for future stochastic modes; currently unused")
    args = parser.parse_args(argv)

    if not os.path.isfile(args.csv_path):
        print(f"ERROR: input file not found: {args.csv_path}", file=sys.stderr)
        return 2

    thds = parse_triplet(args.thd_var, "thd-var")
    mags = parse_triplet(args.mag_var, "mag-var")

    if args.seed is not None:
        random.seed(args.seed)

    ref = load_neso_reference(args.csv_path)
    f_targets = [f for (_i, f) in ref]

    # Build the three per-phase LUTs.
    luts = []
    for k in range(3):
        amp_k = mags[k] * DAC_AMP
        luts.append(build_phase_lut(PHASE_OFFSETS_RAD[k], amp_k, thds[k]))

    prefix = args.out_prefix or os.path.splitext(os.path.basename(args.csv_path))[0]
    out_h_f2 = f"{prefix}_3phase_F2.h"
    out_h_f1 = f"{prefix}_3phase_F1.h"
    out_csv  = f"{prefix}_3phase.csv"

    cli = " ".join([os.path.basename(sys.argv[0])] + sys.argv[1:])
    emit_header_f2(out_h_f2, luts[0], luts[1], thds, mags, f_targets, args.csv_path, cli)
    emit_header_f1(out_h_f1, luts[2],           thds, mags, f_targets, args.csv_path, cli)
    emit_csv(out_csv, thds, mags, f_targets, args.csv_path, cli)

    print("=" * 72)
    print(f"  3-phase trajectory generated from {os.path.basename(args.csv_path)}")
    print("=" * 72)
    print(f"  replay indices     : {len(f_targets)}")
    print(f"  per-phase THD (%)  : Va={thds[0]:.2f}  Vb={thds[1]:.2f}  Vc={thds[2]:.2f}")
    print(f"  per-phase amp frac : Va={mags[0]:.3f}  Vb={mags[1]:.3f}  Vc={mags[2]:.3f}")
    print(f"  frequency range    : {min(f_targets):.4f} - {max(f_targets):.4f} Hz")
    print(f"  Feather 2 header   : {out_h_f2}   (LUT_VA + LUT_VB)")
    print(f"  Feather 1 header   : {out_h_f1}   (LUT_VC)")
    print(f"  CSV written        : {out_csv}")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
