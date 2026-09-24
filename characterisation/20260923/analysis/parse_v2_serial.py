#!/usr/bin/env python3
"""
parse_v2_serial.py - dispersion of f_Vb from an unmodified V2 serial capture.

Run E of the L40 experiment uses Feather1codeMimoDemoV2.ino EXACTLY as
released, so it cannot log raw converter codes - it prints DNP3 frames and a
human-readable block. What it does print, once per loop() iteration, is

    f_Va:49.9987Hz  f_Vb:50.0013Hz  f_Vc:50.0002Hz

This script pulls f_Vb out of that stream and reports its dispersion, which is
the statistic the 5.46 mHz archived floor is quoted in. It is the only run in
the ladder that measures the released build itself, so it is the one that says
whether anything has drifted since June.

Two dispersions are reported and they answer different questions:
  * standard deviation about the mean - comparable with the archived 5.46 mHz
  * successive-difference floor, sd(diff)/sqrt(2) - insensitive to slow drift
    over the capture, which a standard deviation would count as noise
The 13 September README quotes 3.047 and 3.100 mHz for these two on the same
capture, so the gap between them is normally small. A large gap here means the
board drifted during the run and the capture should be retaken.

Usage:
    python parse_v2_serial.py runE_v2.txt
    python parse_v2_serial.py runE_v2.txt --json runE.json
    python parse_v2_serial.py --self-test
"""
from __future__ import annotations
import argparse, json, math, re, statistics, sys

# V2 prints "f_Vb:50.0013Hz"; V6 prints "...,f_Vb=50.0013". Accept both,
# and make the trailing "Hz" optional, so the same tool reads either build.
PAT = re.compile(r"f_Vb\s*[:=]\s*([0-9]+\.?[0-9]*)\s*(?:Hz)?", re.I)
ARCHIVED_FLOOR_MHZ = 5.46          # Section 3.12, archived June capture
SEPT_CLEAN_MHZ     = 3.047         # Section 3.12, estimator on clean samples


def extract(path):
    vals = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = PAT.search(line)
            if m:
                try:
                    v = float(m.group(1))
                except ValueError:
                    continue
                if v > 0.0:            # V2 returns 0.0 when it finds <3 periods
                    vals.append(v)
    return vals


def analyse(vals):
    n = len(vals)
    if n < 3:
        raise SystemExit(f"ERROR: only {n} usable f_Vb values found.")
    mean = statistics.mean(vals)
    sd_mhz = statistics.stdev(vals) * 1000.0
    d = [b - a for a, b in zip(vals, vals[1:])]
    succ_mhz = (statistics.stdev(d) / math.sqrt(2.0)) * 1000.0 if len(d) >= 2 else float("nan")

    # quantisation check: V2 forms the period as an integer number of
    # microseconds, so f moves in steps of about 1e6/T^2 per us
    step_mhz = 1.0e6 / (1.0e6 / mean) ** 2 * 1000.0 if mean > 0 else float("nan")

    out = {
        "n_windows": n,
        "mean_hz": mean,
        "sd_mhz": sd_mhz,
        "successive_difference_mhz": succ_mhz,
        "min_hz": min(vals),
        "max_hz": max(vals),
        "expected_quantisation_step_mhz": step_mhz,
        "reference": {"archived_floor_mhz": ARCHIVED_FLOOR_MHZ,
                      "september_clean_mhz": SEPT_CLEAN_MHZ},
    }
    for key, val in (("sd", sd_mhz), ("successive", succ_mhz)):
        if val == val and val < ARCHIVED_FLOOR_MHZ:
            out[f"residual_vs_archived_{key}_mhz"] = math.sqrt(
                max(ARCHIVED_FLOOR_MHZ ** 2 - val ** 2, 0.0))
    return out


def report(a):
    W = 74
    print("=" * W)
    print("  RUN E - UNMODIFIED V2 BUILD, f_Vb DISPERSION")
    print("=" * W)
    print(f"  Windows parsed:            {a['n_windows']}")
    print(f"  Mean frequency:            {a['mean_hz']:.5f} Hz")
    print(f"  Range:                     {a['min_hz']:.4f} to {a['max_hz']:.4f} Hz")
    print()
    print(f"  Standard deviation:        {a['sd_mhz']:.3f} mHz")
    print(f"  Successive-difference:     {a['successive_difference_mhz']:.3f} mHz")
    print(f"  Expected quantisation step:{a['expected_quantisation_step_mhz']:.3f} mHz")
    print()
    print(f"  Archived June floor:       {a['reference']['archived_floor_mhz']:.2f} mHz")
    print(f"  September clean sampling:  {a['reference']['september_clean_mhz']:.3f} mHz")
    print()
    d = a["sd_mhz"] - a["reference"]["archived_floor_mhz"]
    pct = d / a["reference"]["archived_floor_mhz"] * 100.0
    print(f"  This run vs the archived floor: {d:+.3f} mHz ({pct:+.1f} %)")
    if abs(pct) <= 10.0:
        print("  -> within 10 %: nothing material has drifted since June.")
        print("     Bench state is EXCLUDED as the cause of the residual.")
    elif a["sd_mhz"] < a["reference"]["archived_floor_mhz"]:
        print("  -> materially BELOW the archived floor. The difference is the")
        print("     bench-state term, and it is now measured rather than assumed.")
    else:
        print("  -> materially ABOVE the archived floor. Something is worse than")
        print("     in June; find out what before using any of these runs.")
    print("=" * W)


def self_test():
    import random, tempfile, os
    random.seed(11)
    lines = []
    for _ in range(400):
        f = 50.0 + random.gauss(0, 5.46e-3)
        lines.append(f"  f_Va:49.9990Hz  f_Vb:{f:.4f}Hz  f_Vc:50.0001Hz")
    fd, p = tempfile.mkstemp(suffix=".txt"); os.close(fd)
    open(p, "w").write("\n".join(lines))
    vals = extract(p); a = analyse(vals); os.unlink(p)
    ok = True
    print("SELF-TEST")
    print("-" * 60)
    c1 = len(vals) == 400
    print(f"  [{'PASS' if c1 else 'FAIL'}] parsed 400 of 400 windows -> {len(vals)}")
    c2 = 4.0 < a["sd_mhz"] < 7.0
    print(f"  [{'PASS' if c2 else 'FAIL'}] recovers ~5.46 mHz -> {a['sd_mhz']:.3f} mHz")
    ok = c1 and c2
    print("-" * 60)
    print("ALL PASS" if ok else "FAILURES PRESENT")
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("log", nargs="?")
    ap.add_argument("--json", type=str, default=None)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)
    if args.self_test:
        return self_test()
    if not args.log:
        ap.error("give a serial log, or --self-test")
    a = analyse(extract(args.log))
    report(a)
    if args.json:
        a["input"] = args.log
        with open(args.json, "w", encoding="utf-8", newline="\n") as fh:
            json.dump(a, fh, indent=2)
        print(f"  JSON written: {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
