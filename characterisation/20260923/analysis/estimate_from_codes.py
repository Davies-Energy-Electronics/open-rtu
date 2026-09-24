#!/usr/bin/env python3
"""
estimate_from_codes.py - run the zero-crossing estimator offline over a raw
code log, so the SOURCE's own frequency dispersion can be separated from the
instrument's.

WHY THIS EXISTS
Run E measured the unmodified V2 build against the Feather 2 stimulus and
returned 7.97 mHz where the archived June capture recorded 5.46 mHz. The
ladder (characterise_loop.py over modes 0-3) accounts for the archived figure
to 3.9 %, so the extra 6.00 mHz is something Run E has and the ladder did not:
the ladder ran on a quiet input, Run E ran on a synthesised waveform.

A synthesised 50 Hz reference is not a perfect reference. Feather 2 paces its
waveform with the same microsecond-granular primitives, so its OUTPUT
frequency is quantised and dithers, and that dither adds in quadrature to
anything measuring it. Run E cannot tell the two apart because V2 reports only
its own estimate.

This script can. Given a mode-0 capture (deadline-scheduled, timestamps
uniform to 0.010 us) taken against the SAME stimulus, the only terms left are
the converter noise and the source itself. Subtract the known noise term in
quadrature and what remains is the source.

    sigma_measured^2 = sigma_noise^2 + sigma_source^2   (mode 0, timing ~ 0)

METHOD
Linear-interpolated zero crossings on the de-biased codes; periods formed
between consecutive same-direction crossings; pooled within a window of N
samples exactly as measureFreq() does, capped at MAX_EDGES per direction.
The pooling statistic is selectable because V2 medians and the offline
September estimator meaned, and the two differ by about 25 % on the same data.

INPUT
The same "t_us,code" CSV the characterisation sketch writes. '#' lines are
stripped. Usage:

    python estimate_from_codes.py F_stim_baseline.txt
    python estimate_from_codes.py F_stim_baseline.txt --json F.json
    python estimate_from_codes.py F_stim_baseline.txt --pooling mean
    python estimate_from_codes.py --self-test

Standard library only.

Jack Davies, September 2026
"""

from __future__ import annotations
import argparse, io, json, math, random, statistics, sys

F0_HZ_DEFAULT    = 50.0
WINDOW_N_DEFAULT = 200      # FREQ_WIN in Feather1codeMimoDemoV2.ino
MAX_EDGES        = 12       # MAX_EDGES_PER_TYPE in the same file
VREF             = 3.3
FULL_SCALE       = 4095.0

# reference figures, for the comparison block only
ARCHIVED_FLOOR_MHZ  = 5.46
SEPT_CLEAN_MHZ      = 3.047
RUN_E_MHZ           = 7.973
LADDER_TIMING_MHZ   = 4.271      # mode 3, mean of three records
PREDICTED_CANON_MHZ = 5.247      # sqrt(SEPT_CLEAN^2 + LADDER_TIMING^2)


# ---------------------------------------------------------------------------
def read_log(path):
    """Read t_us,code, tolerating the sketch's '#' header and footer lines."""
    keep = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            keep.append(s)
    if not keep:
        raise SystemExit("ERROR: no data rows in %s" % path)

    import csv as _csv
    buf = io.StringIO("\n".join(keep) + "\n")
    rdr = _csv.reader(buf)
    rows = []
    for r in rdr:
        if len(r) < 2:
            continue
        try:
            rows.append((int(float(r[0])), float(r[1])))
        except ValueError:
            continue          # the header line "t_us,code" lands here
    if len(rows) < 1000:
        raise SystemExit("ERROR: only %d usable rows in %s" % (len(rows), path))
    return rows


def crossings(ts_us, vs, rising=True):
    """Linear-interpolated crossing instants, in microseconds."""
    out = []
    for i in range(1, len(vs)):
        a, b = vs[i - 1], vs[i]
        hit = (a < 0.0 <= b) if rising else (a >= 0.0 > b)
        if not hit:
            continue
        dv = b - a
        frac = (-a / dv) if abs(dv) > 1e-12 else 0.5
        frac = min(1.0, max(0.0, frac))
        out.append(ts_us[i - 1] + frac * (ts_us[i] - ts_us[i - 1]))
    return out


def pool(values, how):
    if not values:
        return None
    return statistics.median(values) if how == "median" else statistics.mean(values)


def estimate(rows, window_n, pooling, f0):
    ts = [r[0] for r in rows]
    codes = [r[1] for r in rows]
    bias = statistics.mean(codes)
    vs = [(c - bias) * (VREF / FULL_SCALE) for c in codes]

    freqs, edge_counts, raw_counts = [], [], []
    for start in range(0, len(rows) - window_n + 1, window_n):
        wt = ts[start:start + window_n]
        wv = vs[start:start + window_n]
        periods, raw = [], 0
        for rising in (True, False):
            allx = crossings(wt, wv, rising)
            raw += len(allx)
            xs = allx[: MAX_EDGES]
            periods += [(xs[i] - xs[i - 1]) * 1e-6 for i in range(1, len(xs))]
        edge_counts.append(len(periods))
        raw_counts.append(raw)          # uncapped, for the is-this-a-signal test
        T = pool(periods, pooling)
        if T and T > 0:
            freqs.append(1.0 / T)

    if len(freqs) < 10:
        raise SystemExit("ERROR: only %d usable windows - is this a sinusoid?"
                         % len(freqs))

    mean_hz = statistics.mean(freqs)
    sd_mhz = statistics.stdev(freqs) * 1000.0
    diffs = [freqs[i] - freqs[i - 1] for i in range(1, len(freqs))]
    succ_mhz = math.sqrt(statistics.mean([d * d for d in diffs]) / 2.0) * 1000.0

    # crossings expected from a clean fundamental in one window
    span_s = (ts[window_n - 1] - ts[0]) * 1e-6
    expected_periods = 2.0 * (f0 * span_s - 1.0)
    expected_raw = 2.0 * f0 * span_s
    noisy = statistics.mean(raw_counts) > 3.0 * max(expected_raw, 1.0)

    return {
        "n_windows": len(freqs),
        "window_n": window_n,
        "pooling": pooling,
        "mean_hz": mean_hz,
        "sd_mhz": sd_mhz,
        "successive_difference_mhz": succ_mhz,
        "sd_uncertainty_mhz": sd_mhz / math.sqrt(2.0 * (len(freqs) - 1)),
        "min_hz": min(freqs),
        "max_hz": max(freqs),
        "periods_per_window_mean": statistics.mean(edge_counts),
        "periods_per_window_expected": expected_periods,
        "raw_crossings_per_window": statistics.mean(raw_counts),
        "raw_crossings_expected": expected_raw,
        "looks_like_noise_not_signal": noisy,
        "bias_code": bias,
        "sd_codes": statistics.stdev(codes),
    }


def report(a, noise_mhz):
    W = 74
    print("=" * W)
    print("  OFFLINE CROSSING ESTIMATOR OVER RAW CODES")
    print("=" * W)
    print(f"  Windows:                {a['n_windows']}  of {a['window_n']} samples")
    print(f"  Pooling:                {a['pooling']} of up to {MAX_EDGES} edges per direction")
    print(f"  Periods per window:     {a['periods_per_window_mean']:.1f} "
          f"(a clean fundamental gives {a['periods_per_window_expected']:.1f})")
    if a["looks_like_noise_not_signal"]:
        print()
        print("  ** The crossing count is far above what a 50 Hz fundamental")
        print("     gives. This looks like a QUIET input, not a stimulus.")
        print("     The dispersion below is not a frequency measurement.")
    print()
    print(f"  Mean frequency:         {a['mean_hz']:.5f} Hz")
    print(f"  Standard deviation:     {a['sd_mhz']:.3f} +/- {a['sd_uncertainty_mhz']:.3f} mHz")
    print(f"  Successive-difference:  {a['successive_difference_mhz']:.3f} mHz")
    r = a["successive_difference_mhz"] / a["sd_mhz"] if a["sd_mhz"] else 0.0
    print(f"  ratio:                  {r:.3f}   "
          + ("white, no drift" if r > 0.9 else "DRIFT PRESENT - the source is wandering"))
    print(f"  Range:                  {a['min_hz']:.4f} to {a['max_hz']:.4f} Hz")
    print()
    print("  -- separating the source from the instrument -----------------")
    m = a["sd_mhz"]
    print(f"  This measurement (mode 0: timing ~ 0.01 us)   {m:.3f} mHz")
    print(f"  Converter noise through the same geometry     {noise_mhz:.3f} mHz")
    if m > noise_mhz:
        src = math.sqrt(m * m - noise_mhz * noise_mhz)
        print(f"  -> SOURCE's own frequency dispersion          {src:.3f} mHz")
        print()
        pred = math.sqrt(PREDICTED_CANON_MHZ ** 2 + src * src)
        print(f"  Predicted V2 total on this stimulus:")
        print(f"    sqrt({PREDICTED_CANON_MHZ:.3f}^2 + {src:.3f}^2) = {pred:.3f} mHz")
        print(f"    Run E measured                             {RUN_E_MHZ:.3f} mHz")
        d = (pred - RUN_E_MHZ) / RUN_E_MHZ * 100.0
        print(f"    agreement                                  {d:+.1f} %")
        print()
        if abs(d) <= 12.0:
            print("  ** RUN E IS EXPLAINED. The excess over the archived 5.46 mHz")
            print("     is the stimulus source, not the instrument. The archived")
            print("     figure stands and the ladder closes it.")
        else:
            print("  ** RUN E IS NOT EXPLAINED by the source alone. Something")
            print("     else differs from June. Do not write the attribution")
            print("     into 3.12 until it is found.")
    else:
        print("  -> the source contributes nothing measurable above the noise.")
        print()
        print("  ** RUN E IS NOT EXPLAINED by the stimulus. Look elsewhere.")
    print("=" * W)


# ---------------------------------------------------------------------------
def _synth(f_hz, n, fs, amp_v, noise_mv, seed, jitter_us=0.0):
    rnd = random.Random(seed)
    rows, t = [], 0.0
    for i in range(n):
        t += 1e6 / fs + (rnd.gauss(0.0, jitter_us) if jitter_us else 0.0)
        v = amp_v * math.sin(2.0 * math.pi * f_hz * t * 1e-6)
        v += rnd.gauss(0.0, noise_mv * 1e-3)
        rows.append((int(round(t)), v / (VREF / FULL_SCALE) + 1872.15))
    return rows


def self_test():
    print("SELF-TEST")
    print("-" * 68)
    ok = True

    # 1. unbiased on a clean tone
    a = estimate(_synth(50.0, 40000, 2000.0, 1.0939, 0.0, 1), 200, "mean", 50.0)
    good = abs(a["mean_hz"] - 50.0) * 1000.0 < 1.0
    ok &= good
    print(f"  [{'PASS' if good else 'FAIL'}] unbiased on a noiseless 50 Hz tone   "
          f"{(a['mean_hz'] - 50.0) * 1000:+.3f} mHz")

    # 2. dispersion matches the analytic propagation.
    #    NOTE the pooling law. Consecutive periods in a chain SHARE crossings,
    #    so the mean of k periods is (t_last - t_first)/k and its error falls
    #    as 1/k, not 1/sqrt(k). Two chains (rising and falling) are averaged,
    #    which divides by a further sqrt(2):
    #        sigma_f = f0^2 * sqrt2 * sigma_t / k_chain / sqrt2
    #                = f0^2 * sigma_t / k_chain
    noise_mv, amp = 0.7389, 1.0939
    S = 2.0 * math.pi * 50.0 * amp
    sigma_t = noise_mv * 1e-3 / S
    a = estimate(_synth(50.0, 120000, 2000.0, amp, noise_mv, 2), 200, "mean", 50.0)
    k_chain = a["periods_per_window_mean"] / 2.0
    want = 50.0 ** 2 * sigma_t / k_chain * 1000.0
    good = abs(a["sd_mhz"] - want) / want < 0.15
    ok &= good
    print(f"  [{'PASS' if good else 'FAIL'}] dispersion vs analytic propagation   "
          f"{a['sd_mhz']:.3f} vs {want:.3f} mHz over {k_chain:.2f} periods/chain")

    # 3. a quiet input is flagged, not silently analysed
    a = estimate(_synth(50.0, 40000, 2000.0, 0.0008, 0.7389, 3), 200, "median", 50.0)
    good = a["looks_like_noise_not_signal"]
    ok &= good
    print(f"  [{'PASS' if good else 'FAIL'}] quiet input flagged as not a signal")

    # 4. an injected source wander comes back out
    rows = _synth(50.0, 120000, 2000.0, amp, 0.0, 4)
    a0 = estimate(rows, 200, "mean", 50.0)
    rows = []
    rnd = random.Random(5)
    t = 0.0
    for w in range(600):                       # 600 windows, each at its own f
        fw = 50.0 + rnd.gauss(0.0, 0.006)      # 6.0 mHz of source wander
        for i in range(200):
            t += 1e6 / 2000.0
            rows.append((int(round(t)),
                         amp * math.sin(2 * math.pi * fw * t * 1e-6)
                         / (VREF / FULL_SCALE) + 1872.15))
    a1 = estimate(rows, 200, "mean", 50.0)
    got = math.sqrt(max(a1["sd_mhz"] ** 2 - a0["sd_mhz"] ** 2, 0.0))
    good = abs(got - 6.0) / 6.0 < 0.30
    ok &= good
    print(f"  [{'PASS' if good else 'FAIL'}] recovers 6.0 mHz of injected source "
          f"wander -> {got:.3f} mHz")

    print("-" * 68)
    print("ALL PASS" if ok else "FAILURES PRESENT")
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("log", nargs="?")
    ap.add_argument("--f0", type=float, default=F0_HZ_DEFAULT)
    ap.add_argument("--window", type=int, default=WINDOW_N_DEFAULT)
    ap.add_argument("--pooling", choices=("median", "mean"), default="mean",
                    help="mean reproduces the offline September estimator; "
                         "median reproduces V2 (default: mean)")
    ap.add_argument("--noise-mhz", type=float, default=SEPT_CLEAN_MHZ,
                    help="converter-noise term to remove in quadrature "
                         f"(default {SEPT_CLEAN_MHZ}, the September figure)")
    ap.add_argument("--json", type=str, default=None)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)

    if args.self_test:
        return self_test()
    if not args.log:
        ap.error("give a log file, or --self-test")

    a = estimate(read_log(args.log), args.window, args.pooling, args.f0)
    report(a, args.noise_mhz)
    if args.json:
        a["input"] = {"file": args.log, "f0_hz": args.f0,
                      "noise_mhz": args.noise_mhz}
        with open(args.json, "w", encoding="utf-8", newline="\n") as fh:
            json.dump(a, fh, indent=2)
        print(f"  JSON written: {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
