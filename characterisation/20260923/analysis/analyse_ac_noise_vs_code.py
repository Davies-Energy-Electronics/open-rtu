#!/usr/bin/env python3
"""
analyse_ac_noise_vs_code.py — is the conditioning chain's noise code-dependent?

WHY THIS EXISTS, AND WHY THE STEPPED-DC SWEEP DOES NOT WORK
-----------------------------------------------------------
The first attempt at this question drove the generator through sixteen DC
levels and measured the spread within each. It cannot work, and the reason is
in Section 2.1.2 of the paper: stage two of the conditioning chain is an
alternating-current coupling, C1 = 1 uF into a ~502 kOhm Thevenin resistance,
giving a first-order high-pass with a time constant near half a second. A DC
level change does not survive it. What the sweep captures is the high-pass
transient — an exponential decay back to the bias point after every step — so
the "within-step spread" is the decay curve, not noise, and it is two orders of
magnitude larger than the noise floor.

The front end passes alternating current by design. So the measurement has to
be made on an alternating-current signal, and the noise has to be separated
from the signal locally rather than by holding the input still.

METHOD
------
A Savitzky-Golay local polynomial is fitted in a short sliding window and
subtracted. Over a 5- or 7-sample window at 2 kHz (2.5-3.5 ms) a 50 Hz sinusoid
is very well approximated by a quadratic, so the residual is noise plus a small
signal-leakage term. Each residual is then binned by the LOCAL FITTED LEVEL —
the instantaneous code the converter was working at — and the standard
deviation is computed per bin. That gives noise against code directly, through
the real signal path, with no direct-current requirement anywhere.

Subtracting a local fit removes some of the noise along with the signal, so the
residual standard deviation is smaller than the true noise by a fixed factor
that depends only on the window length and polynomial order. The factor is
computed exactly from the filter's own coefficients (it is sqrt(1 - 2*h0 +
sum h_i^2) for the centre-tap impulse response h) and is applied as a
correction. --calibrate checks it against a quiet capture of known noise.

USAGE
-----
    # 1. check the estimator against a quiet capture whose SD you know:
    python analyse_ac_noise_vs_code.py --calibrate A1_quiet_single.txt

    # 2. measure noise against code on an AC capture:
    python analyse_ac_noise_vs_code.py B1_ac_running.txt --json ac_noise.json
    python analyse_ac_noise_vs_code.py run1.txt run2.txt run3.txt \\
        --baseline-sd 0.9050 --json ac_noise.json --plot Figure_21.png

Standard library plus numpy; matplotlib only if --plot is given.

Jack Davies, September 2026
"""
from __future__ import annotations
import argparse, json, math, sys
import numpy as np

LSB_V    = 3.3 / 4095.0
WIN      = 9          # Savitzky-Golay window, samples
ORDER    = 4          # polynomial order
# WIN and ORDER were chosen by scanning both against synthetic captures whose
# true noise was known: a 50 Hz sine plus a third harmonic carrying, in one
# case, noise independent of code and, in the other, noise rising 2.19x across
# the range. At (9, 4) the estimator returns the code-independent case as
# 0.945 against a true 0.905 (+4.5 %) with a max/min across bins of 1.07, and
# the code-dependent case as 2.04 against a true 2.19 with a rank correlation
# of +1.00. So the method carries a +4.5 % bias and a 1.07x floor on the
# across-bin ratio, both of which are reported rather than hidden: a measured
# ratio under about 1.15 is the method, not the instrument.
LEAK_BIAS  = 1.045    # measured overestimate of the pooled figure
LEAK_RATIO = 1.07     # across-bin ratio the method produces on flat noise
NBINS    = 16         # code bins, to match the sixteen levels the sweep intended
MIN_PER_BIN = 400     # ignore a bin with fewer residuals than this
EDGE_TRIM   = 50      # drop this many samples at each end of the record


def read_capture(path):
    codes, meta = [], {}
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            t = line.strip()
            if not t:
                continue
            if t.startswith("#"):
                b = t.lstrip("#").strip()
                if "," in b:
                    k, _, v = b.partition(",")
                    meta.setdefault(k.strip(), v.strip())
                continue
            try:
                codes.append(int(t))
            except ValueError:
                pass
    if not codes:
        sys.exit(f"{path}: no ADC codes found")
    return np.asarray(codes, dtype=float), meta


def sg_coeffs(win, order):
    """Centre-tap Savitzky-Golay smoothing coefficients."""
    half = win // 2
    x = np.arange(-half, half + 1, dtype=float)
    A = np.vander(x, order + 1, increasing=True)
    # row of (A^T A)^-1 A^T giving the constant term = the smoothed centre value
    return np.linalg.pinv(A)[0]


def residual_scale(h):
    """
    The residual at the centre tap is  r = y0 - sum h_i y_i.
    For independent noise of variance s^2 this has variance
        s^2 * (1 - 2*h0 + sum h_i^2),
    so the residual SD understates the true noise SD by the square root of that.
    """
    half = len(h) // 2
    return math.sqrt(1.0 - 2.0 * h[half] + float(np.sum(h ** 2)))


def local_residuals(codes, win=WIN, order=ORDER):
    """Return (fitted level, residual) for every interior sample."""
    h = sg_coeffs(win, order)
    smooth = np.convolve(codes, h[::-1], mode="same")
    half = win // 2
    s = slice(half + EDGE_TRIM, len(codes) - half - EDGE_TRIM)
    return smooth[s], (codes[s] - smooth[s]), residual_scale(h)


def analyse(path, nbins=NBINS):
    codes, meta = read_capture(path)
    level, resid, scale = local_residuals(codes)
    # Drop samples adjacent to a railed reading: a clipped sample is not noise.
    ok = (codes[WIN // 2 + EDGE_TRIM: len(codes) - WIN // 2 - EDGE_TRIM] > 0) & \
         (codes[WIN // 2 + EDGE_TRIM: len(codes) - WIN // 2 - EDGE_TRIM] < 4095)
    level, resid = level[ok], resid[ok]

    lo, hi = np.percentile(level, [0.5, 99.5])
    edges = np.linspace(lo, hi, nbins + 1)
    idx = np.clip(np.digitize(level, edges) - 1, 0, nbins - 1)

    rows = []
    for b in range(nbins):
        m = idx == b
        n = int(m.sum())
        if n < MIN_PER_BIN:
            continue
        sd = float(resid[m].std(ddof=1)) / scale
        rows.append({"bin": b, "n": n,
                     "mean_code": float(level[m].mean()),
                     "sd_code": sd,
                     "sd_millivolts": sd * LSB_V * 1000})
    return {"file": path, "label": meta.get("label", "?"),
            "mode": meta.get("mode", "?"),
            "fs_hz": float(meta.get("fs_hz", 2000.0)),
            "scale": scale, "code_span": [float(lo), float(hi)],
            "bins": rows}


def report(res, baseline_sd=None):
    print("=" * 78)
    print(f"  {res['file']}   label={res['label']}   mode={res['mode']}"
          f"   fs={res['fs_hz']:.2f} Hz")
    print(f"  code span {res['code_span'][0]:.0f} to {res['code_span'][1]:.0f}"
          f"   ({len(res['bins'])} usable bins of {NBINS})")
    print("=" * 78)
    if not res["bins"]:
        print("  NO USABLE BINS. The record does not traverse a range of codes -")
        print("  is the generator running?")
        return None
    print(f"{'bin':>4}{'n':>8}{'mean code':>11}{'volts':>9}"
          f"{'SD code':>9}{'SD mV':>8}")
    print("-" * 78)
    for r in res["bins"]:
        print(f"{r['bin']:>4}{r['n']:>8}{r['mean_code']:>11.1f}"
              f"{r['mean_code'] * LSB_V:>9.4f}{r['sd_code']:>9.4f}"
              f"{r['sd_millivolts']:>8.4f}")
    print("-" * 78)

    sds = np.array([r["sd_code"] for r in res["bins"]])
    means = np.array([r["mean_code"] for r in res["bins"]])
    ns = np.array([r["n"] for r in res["bins"]], dtype=float)
    pooled = math.sqrt(float(np.sum((ns - 1) * sds ** 2) / np.sum(ns - 1)))
    lo, hi = sds.min(), sds.max()
    ratio = hi / lo if lo > 0 else float("inf")
    rm = np.argsort(np.argsort(means)).astype(float)
    rs = np.argsort(np.argsort(sds)).astype(float)
    rho = float(np.corrcoef(rm, rs)[0, 1]) if len(sds) > 2 else float("nan")

    print(f"  pooled across bins      {pooled:.4f} code "
          f"({pooled * LSB_V * 1000:.4f} mV)")
    print(f"  range across bins       {lo:.4f} to {hi:.4f} code "
          f"(max/min = {ratio:.2f}x)")
    print(f"  rank correlation with code   rho = {rho:+.3f}")
    print(f"  method floor: this estimator returns {LEAK_RATIO:.2f}x and a "
          f"{100*(LEAK_BIAS-1):.1f} % high pooled")
    print(f"                figure on synthetic noise that is genuinely flat.")
    print(f"  leakage-corrected pooled {pooled / LEAK_BIAS:.4f} code "
          f"({pooled / LEAK_BIAS * LSB_V * 1000:.4f} mV)")
    if baseline_sd is not None:
        print(f"  quiet-input baseline    {baseline_sd:.4f} code"
              f"   -> corrected/baseline = "
              f"{pooled / LEAK_BIAS / baseline_sd:.3f}x")
    print()
    if ratio < 1.15 and abs(rho) < 0.50:
        print("  VERDICT: the noise floor is NOT materially code-dependent.")
        print("  The single quiet-input figure used in the error budget applies")
        print("  across the working range.")
    elif ratio >= 1.15 and abs(rho) >= 0.50:
        print("  VERDICT: the noise floor IS code-dependent and trends with code.")
        print("  Quote the pooled figure and say which end of the range is worst.")
    else:
        print("  VERDICT: mixed - the spread across bins and the trend with code")
        print("  disagree. Usually one bin is an outlier; read the table.")
    print()
    return {"pooled_sd_code": pooled, "min_sd_code": float(lo),
            "max_sd_code": float(hi), "ratio_max_min": float(ratio),
            "rank_corr_with_code": rho,
            "pooled_sd_millivolts": pooled * LSB_V * 1000,
            "pooled_sd_code_corrected": pooled / LEAK_BIAS,
            "pooled_sd_millivolts_corrected": pooled / LEAK_BIAS * LSB_V * 1000,
            "method_ratio_floor": LEAK_RATIO, "method_bias": LEAK_BIAS}


def calibrate(path):
    """Check the estimator against a quiet capture whose true SD is known."""
    codes, meta = read_capture(path)
    true_sd = codes.std(ddof=1)
    level, resid, scale = local_residuals(codes)
    est = resid.std(ddof=1) / scale
    print("=" * 78)
    print(f"  CALIBRATION against {path}")
    print("=" * 78)
    print(f"  whole-record SD (the truth, input is quiet)   {true_sd:.4f} code")
    print(f"  Savitzky-Golay window {WIN}, order {ORDER}, "
          f"residual scale {scale:.4f}")
    print(f"  estimator output                             {est:.4f} code")
    print(f"  error                                        "
          f"{100 * (est - true_sd) / true_sd:+.2f} %")
    print()
    if abs(est - true_sd) / true_sd < 0.05:
        print("  PASS - the estimator recovers a known noise SD to better than 5 %.")
    else:
        print("  FAIL - do not trust the code-dependence figures until this passes.")
    print()
    return abs(est - true_sd) / true_sd < 0.05


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("captures", nargs="*")
    ap.add_argument("--calibrate", default=None,
                    help="quiet capture to validate the estimator against")
    ap.add_argument("--bins", type=int, default=NBINS)
    ap.add_argument("--baseline-sd", type=float, default=None)
    ap.add_argument("--json", default=None)
    ap.add_argument("--plot", default=None)
    a = ap.parse_args(argv)

    if a.calibrate:
        ok = calibrate(a.calibrate)
        if not a.captures:
            return 0 if ok else 1

    out = {"captures": [], "summary": {}}
    for p in a.captures:
        res = analyse(p, a.bins)
        res["summary"] = report(res, a.baseline_sd)
        out["captures"].append(res)

    good = [c for c in out["captures"] if c.get("summary")]
    if len(good) > 1:
        pooled = math.sqrt(float(np.mean([c["summary"]["pooled_sd_code"] ** 2
                                          for c in good])))
        ratios = [c["summary"]["ratio_max_min"] for c in good]
        print("=" * 78)
        print(f"  ACROSS {len(good)} CAPTURES: pooled {pooled:.4f} code "
              f"({pooled * LSB_V * 1000:.4f} mV)")
        print(f"  max/min across bins, per capture: "
              f"{', '.join(f'{r:.2f}x' for r in ratios)}")
        print("=" * 78)
        out["summary"] = {"n_captures": len(good), "pooled_sd_code": pooled,
                          "pooled_sd_millivolts": pooled * LSB_V * 1000,
                          "ratio_max_min_per_capture": ratios}
    elif good:
        out["summary"] = good[0]["summary"]

    if a.json:
        with open(a.json, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"  wrote {a.json}")

    if a.plot and good:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.2, 4.2), dpi=160)
        for c in good:
            ax.plot([r["mean_code"] for r in c["bins"]],
                    [r["sd_code"] for r in c["bins"]],
                    "o-", lw=1.2, ms=4, label=c["label"])
        if a.baseline_sd is not None:
            ax.axhline(a.baseline_sd, ls="--", lw=1.0, color="0.4",
                       label=f"quiet input, {a.baseline_sd:.3f} code")
        ax.set_xlabel("converter code")
        ax.set_ylabel("input-referred noise (ADC codes)")
        ax.set_title("Noise against converter code, measured on an AC signal")
        ax.set_ylim(bottom=0)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(a.plot)
        print(f"  wrote {a.plot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
