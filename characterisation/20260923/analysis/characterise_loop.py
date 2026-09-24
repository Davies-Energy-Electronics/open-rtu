#!/usr/bin/env python3
"""
characterise_loop.py - close L40 and L41 from one bench capture.

Takes the per-sample log produced by the loop-timing characterisation sketch
and answers both open questions from the same data:

  L41  What is the estimator's TRUE per-channel sampling rate?
       Measured as the mean inter-sample interval, with its uncertainty,
       rather than derived from a conversion time plus a coded delay.

  L40  How much of the 4.53 mHz unattributed residual is sample-instant
       jitter?  Measured as the standard deviation of that same interval,
       propagated through the estimator's own geometry.

WHY ONE SCRIPT ANSWERS BOTH
The two questions are the first and second moments of the same measurement.
Section 3.12 derived the sampling rate (39.21 us conversion + 500 us coded
delay = 539.2 us, hence "at most 1855 Hz") and never measured it; and it
named loop jitter as a candidate for the residual without measuring that
either. The inter-sample interval gives both: its mean is the rate, its
standard deviation is the jitter.

THE PROPAGATION, DERIVED
A zero-crossing estimator with sub-sample linear interpolation converts an
amplitude error into a time error through the slope of the waveform at the
crossing:

    S       = 2 * pi * f0 * A                    crossing slope, V/s
    sigma_t = sigma_v / S                        crossing-instant error
    sigma_T = sqrt(2) * sigma_t                  period, two crossings
    sigma_f = f0^2 * sigma_T                     since f = 1/T

A timing error on the sample instant adds to sigma_t directly and in
quadrature, because it displaces the crossing along the same axis:

    sigma_t_total = sqrt( (sigma_v/S)^2 + sigma_j^2 )

This chain reproduces the paper's own quantisation figures exactly: at the
ideal 12-bit step of 0.233 mV it gives 0.700 us, 0.989 us and 2.47 mHz
against the 0.698 us, 0.987 us and 2.47 mHz stated in Section 2.3. That
agreement is the reason the chain is trusted to propagate the measured
jitter.

INPUT
A CSV with at least these columns, one row per frequency-channel sample:
    t_us    integer microsecond timestamp of the sample
    code    raw converter code (optional; used for the noise cross-check)
Extra columns are ignored. The sketch that produces it is described in the
L40/L41 protocol document delivered alongside this script.

Usage:
    python characterise_loop.py loop_log.csv
    python characterise_loop.py loop_log.csv --json loop_out.json
    python characterise_loop.py loop_log.csv --amplitude 1.060 --f0 50.0
    python characterise_loop.py --self-test

Standard library only.
"""

from __future__ import annotations
import argparse, csv, io, json, math, statistics, sys

# --- constants from the manuscript ------------------------------------------
F0_HZ_DEFAULT      = 50.0      # nominal fundamental
AMPLITUDE_V_DEFAULT = 1.060    # conditioned amplitude at the converter input (S2.3)
IDEAL_LSB_MV       = 0.233     # ideal 12-bit quantisation step, Section 3.12
QUIET_NOISE_MV     = 0.659     # measured quiet-input noise, Section 3.12
BROADBAND_MV       = 1.399     # measured broadband noise with stimulus, Section 3.12
ESTIMATOR_CHAR_MHZ = 3.05      # estimator over the characterisation samples
ARCHIVED_FLOOR_MHZ = 5.46      # floor realised in the archived replay capture
STATED_FS_HZ       = 2000.0    # the nominal figure stated in Sections 2.1.2 and 2.3
DERIVED_FS_HZ      = 1854.6    # 1e6 / (39.21 + 500) us, the Section 3.12 derivation


def slope_v_per_s(f0: float, amplitude_v: float) -> float:
    return 2.0 * math.pi * f0 * amplitude_v


def mhz_from_sigma_t(sigma_t_s: float, f0: float) -> float:
    """Frequency dispersion in mHz from a crossing-instant error in seconds."""
    return f0 * f0 * math.sqrt(2.0) * sigma_t_s * 1000.0


def mhz_from_sigma_v(sigma_v_v: float, f0: float, amplitude_v: float) -> float:
    return mhz_from_sigma_t(sigma_v_v / slope_v_per_s(f0, amplitude_v), f0)


def read_log(path: str):
    """Read the per-sample log.

    The sketch writes '#' metadata lines BEFORE the "t_us,code" header, so the
    '#' lines are stripped first - csv.DictReader would otherwise take the
    first of them as the column header and find no timestamp column. Fixed
    23 Sep 2026 against a real capture.
    """
    rows = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        body = "".join(ln for ln in fh if not ln.lstrip().startswith("#"))
    with io.StringIO(body, newline="") as fh:
        reader = csv.DictReader(fh)
        fmap = {f.lower().strip(): f for f in reader.fieldnames or []}

        def col(*names):
            for n in names:
                if n in fmap:
                    return fmap[n]
            return None

        tcol = col("t_us", "micros", "timestamp_us", "t")
        ccol = col("code", "adc", "raw", "counts")
        if tcol is None:
            raise SystemExit(
                "ERROR: no timestamp column found. Expected one of "
                "t_us / micros / timestamp_us / t."
            )
        for r in reader:
            try:
                t = int(float(r[tcol]))
            except (TypeError, ValueError):
                continue
            c = None
            if ccol is not None:
                try:
                    c = float(r[ccol])
                except (TypeError, ValueError):
                    c = None
            rows.append((t, c))
    return rows


def intervals(rows):
    """Inter-sample intervals in microseconds, wrap-safe for a 32-bit micros()."""
    out = []
    WRAP = 1 << 32
    for (a, _), (b, _) in zip(rows, rows[1:]):
        d = b - a
        if d < 0:
            d += WRAP
        out.append(d)
    return out


def robust_stats(xs):
    """Median, MAD-based sigma, and the ordinary mean and sd.

    The MAD figure is reported alongside the standard deviation because a
    scheduling loop produces occasional long intervals - a housekeeping pass,
    a Wi-Fi task - and those outliers inflate a standard deviation without
    representing the jitter that acts on every sample. Both are given and
    neither is chosen for you.
    """
    n = len(xs)
    med = statistics.median(xs)
    mad = statistics.median([abs(x - med) for x in xs])
    sigma_mad = 1.4826 * mad
    mean = statistics.mean(xs)
    sd = statistics.stdev(xs) if n >= 2 else 0.0
    return {"n": n, "median": med, "mad_sigma": sigma_mad, "mean": mean, "sd": sd,
            "min": min(xs), "max": max(xs)}


def decompose_alternation(iv):
    """Split the interval series into a deterministic period-2 component and
    the random residual.

    Measured on this bench 23 Sep 2026: a mode-0 capture can show a perfect
    two-level alternation of +/-4 us, 8191 intervals of each, arising from the
    converter's conversion time alternating between two values. It is NOT
    jitter. It is deterministic, it displaces alternate samples by a fixed
    amount, and it largely cancels over many crossings - whereas a standard
    deviation counts it at full weight and reports 14 mHz of phantom noise
    against a residual we are trying to resolve at 4.5 mHz.

    So both are reported: the alternation amplitude, and the RANDOM jitter
    that remains once it is removed. The second is the one that propagates.
    """
    ev = iv[0::2]
    od = iv[1::2]
    if len(ev) < 2 or len(od) < 2:
        return 0.0, statistics.stdev(iv) if len(iv) > 1 else 0.0
    alt = (statistics.mean(ev) - statistics.mean(od)) / 2.0
    resid = [x - (alt if i % 2 == 0 else -alt) for i, x in enumerate(iv)]
    return abs(alt), statistics.stdev(resid)


def analyse(rows, f0, amplitude_v):
    if len(rows) < 3:
        raise SystemExit("ERROR: need at least three samples.")
    iv = intervals(rows)
    st = robust_stats(iv)
    alt_us, rand_us = decompose_alternation(iv)

    fs_mean = 1.0e6 / st["mean"]
    fs_med = 1.0e6 / st["median"]

    # jitter propagated through the estimator geometry
    sj_sd = st["sd"] * 1e-6
    sj_mad = st["mad_sigma"] * 1e-6
    mhz_sd = mhz_from_sigma_t(sj_sd, f0)
    mhz_mad = mhz_from_sigma_t(sj_mad, f0)

    # what the residual demands
    residual = math.sqrt(ARCHIVED_FLOOR_MHZ ** 2 - ESTIMATOR_CHAR_MHZ ** 2)
    sigma_t_needed = residual / 1000.0 / (f0 * f0) / math.sqrt(2.0)

    # optional noise cross-check from the codes
    noise = None
    codes = [c for _, c in rows if c is not None]
    if len(codes) >= 32:
        noise = {"n": len(codes), "sd_codes": statistics.stdev(codes),
                 "mean_codes": statistics.mean(codes)}

    return {
        "alternation": {
            "amplitude_us": alt_us,
            "random_jitter_us": rand_us,
            "mhz_from_alternation": mhz_from_sigma_t(alt_us * 1e-6, f0),
            "mhz_from_random_jitter": mhz_from_sigma_t(rand_us * 1e-6, f0),
        },
        "interval_us": st,
        "sampling_rate_hz": {"from_mean": fs_mean, "from_median": fs_med,
                             "stated_in_paper": STATED_FS_HZ,
                             "derived_in_3_12": DERIVED_FS_HZ},
        "jitter": {"sigma_sd_us": st["sd"], "sigma_mad_us": st["mad_sigma"],
                   "mhz_from_sd": mhz_sd, "mhz_from_mad": mhz_mad},
        "residual": {"archived_floor_mhz": ARCHIVED_FLOOR_MHZ,
                     "estimator_char_mhz": ESTIMATOR_CHAR_MHZ,
                     "unattributed_mhz": residual,
                     "sigma_t_that_would_explain_it_us": sigma_t_needed * 1e6},
        "noise_cross_check": noise,
    }


def report(a, f0, amplitude_v):
    S = slope_v_per_s(f0, amplitude_v)
    st = a["interval_us"]; fs = a["sampling_rate_hz"]; j = a["jitter"]; r = a["residual"]
    W = 78
    print("=" * W)
    print("  ACQUISITION-LOOP CHARACTERISATION  (L40 jitter, L41 sampling rate)")
    print("=" * W)
    print(f"  Samples logged:        {st['n'] + 1}")
    print(f"  Crossing slope:        S = 2*pi*{f0:.1f}*{amplitude_v:.3f} = {S:.1f} V/s")
    print()
    print("  -- L41  TRUE SAMPLING RATE, MEASURED -------------------------------")
    print(f"  Mean interval:         {st['mean']:.3f} us   ->  {fs['from_mean']:.1f} Hz")
    print(f"  Median interval:       {st['median']:.3f} us   ->  {fs['from_median']:.1f} Hz")
    print(f"  Range:                 {st['min']} to {st['max']} us")
    print(f"  Stated in the paper:   {fs['stated_in_paper']:.0f} Hz")
    print(f"  Derived in S3.12:      {fs['derived_in_3_12']:.1f} Hz")
    d_stated = (fs['from_mean'] - fs['stated_in_paper']) / fs['stated_in_paper'] * 100
    d_derived = (fs['from_mean'] - fs['derived_in_3_12']) / fs['derived_in_3_12'] * 100
    print(f"  Measured vs stated:    {d_stated:+.2f} %")
    print(f"  Measured vs derived:   {d_derived:+.2f} %")
    print()
    print("  -- L40  SAMPLE-INSTANT JITTER, MEASURED ----------------------------")
    al = a["alternation"]
    print(f"  sigma (standard dev):  {j['sigma_sd_us']:.4f} us  ->  {j['mhz_from_sd']:.3f} mHz")
    print(f"  sigma (MAD-based):     {j['sigma_mad_us']:.4f} us  ->  {j['mhz_from_mad']:.3f} mHz")
    print()
    print("     decomposed:")
    print(f"       deterministic period-2 alternation  +/-{al['amplitude_us']:.4f} us"
          f"  ({al['mhz_from_alternation']:.3f} mHz if counted as noise)")
    print(f"       RANDOM jitter after removing it      {al['random_jitter_us']:.4f} us"
          f"  ->  {al['mhz_from_random_jitter']:.4f} mHz   <- the one that propagates")
    if al["amplitude_us"] > 5 * max(al["random_jitter_us"], 1e-9):
        print("       The interval series is dominated by a DETERMINISTIC two-level")
        print("       alternation, not by jitter. Use the random figure.")
    print()
    print(f"  Unattributed residual: {r['unattributed_mhz']:.3f} mHz")
    print(f"    = sqrt({r['archived_floor_mhz']:.2f}^2 - {r['estimator_char_mhz']:.2f}^2), both measured")
    print(f"  Jitter that would explain it in full: {r['sigma_t_that_would_explain_it_us']:.3f} us")
    print()
    for label, mhz in (("standard deviation", j["mhz_from_sd"]),
                       ("MAD-based sigma", j["mhz_from_mad"]),
                       ("random jitter", a["alternation"]["mhz_from_random_jitter"])):
        frac = mhz / r["unattributed_mhz"]
        if mhz >= r["unattributed_mhz"]:
            verdict = "accounts for the WHOLE residual"
            left = 0.0
        else:
            left = math.sqrt(r["unattributed_mhz"] ** 2 - mhz ** 2)
            verdict = f"leaves {left:.2f} mHz still unattributed"
        print(f"  Using the {label:19s}: {frac * 100:5.1f} % of the residual in amplitude, {verdict}")
    nc = a["noise_cross_check"]
    if nc:
        print()
        print("  -- noise cross-check from the logged codes -------------------------")
        print(f"  {nc['n']} codes, mean {nc['mean_codes']:.2f}, sd {nc['sd_codes']:.3f} codes")
    print("=" * W)


def self_test():
    """Verify the propagation against the manuscript's own published figures."""
    f0, A = F0_HZ_DEFAULT, AMPLITUDE_V_DEFAULT
    S = slope_v_per_s(f0, A)
    checks = []
    st = IDEAL_LSB_MV * 1e-3 / S
    checks.append(("crossing-instant error at the ideal LSB", st * 1e6, 0.698, "us", 0.01))
    checks.append(("period error at the ideal LSB", math.sqrt(2) * st * 1e6, 0.987, "us", 0.01))
    checks.append(("frequency error at the ideal LSB",
                   mhz_from_sigma_v(IDEAL_LSB_MV * 1e-3, f0, A), 2.47, "mHz", 0.02))
    checks.append(("crossing slope", S, 333.0, "V/s", 0.5))
    checks.append(("unattributed residual",
                   math.sqrt(ARCHIVED_FLOOR_MHZ ** 2 - ESTIMATOR_CHAR_MHZ ** 2),
                   4.53, "mHz", 0.01))
    ok = True
    print("SELF-TEST against the published figures in Sections 2.3 and 3.12")
    print("-" * 72)
    for name, got, want, unit, tol in checks:
        good = abs(got - want) <= tol
        ok &= good
        print(f"  [{'PASS' if good else 'FAIL'}] {name:46s} {got:8.3f} vs {want:7.3f} {unit}")
    print("-" * 72)
    print("ALL PASS" if ok else "FAILURES PRESENT")
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("log", nargs="?", help="per-sample CSV from the timing sketch")
    ap.add_argument("--f0", type=float, default=F0_HZ_DEFAULT)
    ap.add_argument("--amplitude", type=float, default=AMPLITUDE_V_DEFAULT,
                    help=f"conditioned amplitude in V (default {AMPLITUDE_V_DEFAULT})")
    ap.add_argument("--json", type=str, default=None)
    ap.add_argument("--self-test", action="store_true",
                    help="check the propagation against the published figures and exit")
    args = ap.parse_args(argv)

    if args.self_test:
        return self_test()
    if not args.log:
        ap.error("give a log file, or --self-test")

    rows = read_log(args.log)
    a = analyse(rows, args.f0, args.amplitude)
    report(a, args.f0, args.amplitude)
    if args.json:
        a["input"] = {"file": args.log, "f0_hz": args.f0, "amplitude_v": args.amplitude}
        with open(args.json, "w", encoding="utf-8", newline="\n") as fh:
            json.dump(a, fh, indent=2)
        print(f"  JSON written: {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
