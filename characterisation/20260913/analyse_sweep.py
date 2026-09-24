#!/usr/bin/env python3
"""
analyse_sweep.py — input-referred noise as a function of input level.

Reads a characterise.ino capture taken while Feather 2 runs
feather2_dc_sweep.ino, finds the DC plateaus by detection, and reports the
ADC's mean and noise at each level.

Answers three things the single-point grounded capture cannot:
  1. Is the noise code-dependent, or flat across the range? Code-dependent
     noise is differential non-linearity; flat noise is additive.
  2. Where does the chain actually saturate? The ESP32 ADC at 11 dB
     attenuation is usable over roughly 150-2450 mV, not 0-3.3 V, and the
     signal swings +/-1.06 V about the bias.
  3. What is sigma_v at the operating point specifically, as opposed to
     averaged over a range the instrument never uses?

No synchronisation is needed: the ADC's own mean at each plateau is the
x-axis, so the analyser never needs to know which DAC code produced it.

Usage:
    python analyse_sweep.py C1_sweep.txt
    python analyse_sweep.py C1_sweep.txt --json sweep.json --csv sweep.csv

Standard library plus numpy.
"""
from __future__ import annotations
import argparse, json, math, sys
import numpy as np

LSB_V     = 3.3 / 4095.0
SIGMA_Q_V = LSB_V / math.sqrt(12)
FREQ_WIN  = 200
AMPLITUDE = 1.060          # V, conditioned amplitude at the ADC (Sec 2.2)


def read_capture(path):
    codes, meta = [], {}
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                b = s.lstrip("#").strip()
                if "," in b:
                    k, _, v = b.partition(",")
                    meta[k.strip()] = v.strip()
                continue
            try:
                codes.append(int(s))
            except ValueError:
                pass
    if not codes:
        sys.exit(f"{path}: no ADC codes found")
    return np.asarray(codes, dtype=float), meta


def median_filter(x, w=9):
    if w < 3:
        return x
    pad = w // 2
    xp = np.pad(x, pad, mode="edge")
    return np.median(np.lib.stride_tricks.sliding_window_view(xp, w), axis=-1)


def find_plateaus(codes, min_len, jump_thresh, trim=0.2):
    """Segment on large level changes in a median-filtered copy, then trim
    the settling transient from each end."""
    sm = median_filter(codes, 9)
    d = np.abs(np.diff(sm))
    edges = np.flatnonzero(d > jump_thresh)

    bounds = [0] + [int(e) + 1 for e in edges] + [len(codes)]
    segs = []
    for a, b in zip(bounds[:-1], bounds[1:]):
        if b - a < min_len:
            continue
        k = int((b - a) * trim)
        seg = codes[a + k: b - k]
        if len(seg) >= 20:
            segs.append((a + k, b - k, seg))
    return segs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("capture")
    ap.add_argument("--dwell-ms", type=float, default=128.0)
    ap.add_argument("--jump", type=float, default=10.0,
                    help="ADC codes of change that marks a new level (default 10)")
    ap.add_argument("--json", default=None)
    ap.add_argument("--csv", default=None)
    a = ap.parse_args(argv)

    codes, meta = read_capture(a.capture)
    fs = float(meta.get("fs_hz", 2000.0))
    expected = int(fs * a.dwell_ms / 1000.0)
    segs = find_plateaus(codes, min_len=max(30, int(expected * 0.5)), jump_thresh=a.jump)

    if len(segs) < 5:
        sys.exit(f"\n  Only {len(segs)} plateaus found. Either the sweep was not "
                 f"running on Feather 2,\n  or --jump needs lowering. "
                 f"Capture range was {codes.min():.0f}..{codes.max():.0f} codes.\n")

    rows = []
    for i0, i1, seg in segs:
        m = seg.mean()
        s = seg.std(ddof=1)
        rows.append({"mean_code": float(m), "sd_code": float(s),
                     "mean_volts": float(m * LSB_V),
                     "sd_millivolts": float(s * LSB_V * 1e3),
                     "n": int(len(seg))})
    rows.sort(key=lambda r: r["mean_code"])

    print(f"\n{'='*70}")
    print(f"  STEPPED-DC SWEEP — {meta.get('label', a.capture)}")
    print(f"{'='*70}")
    print(f"  samples {len(codes)}   fs {fs:.1f} Hz   plateaus found {len(rows)}")
    print(f"  capture spans codes {codes.min():.0f}..{codes.max():.0f}"
          f"  ({codes.min()*LSB_V:.3f}..{codes.max()*LSB_V:.3f} V)")
    print(f"  {'-'*66}")
    print(f"  {'mean code':>10} {'volts':>8} {'sd code':>9} {'sd mV':>8} "
          f"{'x ideal':>8}  {'n':>5}")
    for r in rows:
        print(f"  {r['mean_code']:>10.1f} {r['mean_volts']:>8.4f} "
              f"{r['sd_code']:>9.3f} {r['sd_millivolts']:>8.3f} "
              f"{r['sd_millivolts']/1e3/SIGMA_Q_V:>8.2f}  {r['n']:>5}")

    sds = np.array([r["sd_millivolts"] for r in rows])
    mns = np.array([r["mean_code"] for r in rows])

    # The instrument's own operating window: bias +/- the conditioned amplitude.
    bias = float(np.median(mns))
    amp_codes = AMPLITUDE / LSB_V
    inwin = (mns > bias - amp_codes) & (mns < bias + amp_codes)
    # Saturation: plateaus where the ADC has stopped following the input.
    sat_lo = mns < 190
    sat_hi = mns > 3040

    print(f"  {'-'*66}")
    print(f"  median sd, all levels           {np.median(sds):.3f} mV"
          f"  ({np.median(sds)/1e3/SIGMA_Q_V:.2f}x ideal)")
    if inwin.sum():
        print(f"  median sd, operating window     {np.median(sds[inwin]):.3f} mV"
              f"  ({inwin.sum()} levels)   <-- the figure for Sec 2.8.6")
    print(f"  spread across levels            {sds.min():.3f} to {sds.max():.3f} mV")
    ratio = sds.max() / max(sds.min(), 1e-9)
    print(f"  max/min                         {ratio:.2f}x")
    if sat_lo.sum() or sat_hi.sum():
        print(f"  levels outside the ADC's usable window: "
              f"{int(sat_lo.sum())} low, {int(sat_hi.sum())} high")

    print(f"{'='*70}")
    if ratio < 1.5:
        print("  >> NOISE IS FLAT across the range. It is additive — front-end")
        print("     or reference noise — not differential non-linearity. Report")
        print("     a single sigma_v and use it for the bound.")
    else:
        print("  >> NOISE IS CODE-DEPENDENT. Part of it is converter DNL, which")
        print("     is exactly the mechanism Section 4.1 already blames for the")
        print("     harmonic result — so one cause now accounts for both, with")
        print("     a measurement behind it rather than an inference.")
    print()

    if a.csv:
        with open(a.csv, "w") as fh:
            fh.write("mean_code,mean_volts,sd_code,sd_millivolts,n\n")
            for r in rows:
                fh.write(f"{r['mean_code']:.3f},{r['mean_volts']:.6f},"
                         f"{r['sd_code']:.4f},{r['sd_millivolts']:.4f},{r['n']}\n")
        print(f"  -> {a.csv}")
    if a.json:
        out = {"label": meta.get("label", a.capture), "fs_hz": fs,
               "n_plateaus": len(rows),
               "median_sd_millivolts_all": float(np.median(sds)),
               "median_sd_millivolts_operating": (
                   float(np.median(sds[inwin])) if inwin.sum() else None),
               "sd_spread_ratio": float(ratio),
               "levels": rows}
        with open(a.json, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"  -> {a.json}")
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
