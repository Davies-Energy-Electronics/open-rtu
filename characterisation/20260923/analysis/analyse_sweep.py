#!/usr/bin/env python3
"""
analyse_sweep.py — noise against converter code, from a SHORT-dwell staircase.

WHAT CHANGED, AND WHY
---------------------
The first version of this script assumed each staircase level was a settled
direct-current level. It is not. Stage two of the conditioning chain is an
alternating-current coupling (C1 = 1 uF; Section 2.1.2), so a direct-current
level never reaches the converter: each step is passed as an edge and then
decays back to the bias point with a time constant measured on this bench at
533.7 +/- 21.8 ms. At the 1000 ms dwell first used, each level decayed 85 % of
the way back before the next arrived, and what the script measured was the
decay curve rather than the noise.

The fix is to make the dwell SHORT against that time constant. At 25 ms the
level falls by only 4.6 % across a dwell, which a linear detrend inside the
dwell removes almost exactly, and the staircase passes through the coupling
essentially intact. Validated in simulation against the measured time constant:
with genuinely flat noise of 0.905 codes the estimator returns 0.913 (+0.9 %)
with a spread across bins of 1.09x, and with noise rising 2.8x across the range
it returns 2.98x with a rank correlation of +1.00. At 250 ms and above it fails,
which is what the bench showed.

A short dwell also removes the objection that sank the alternating-current
method: inside a dwell the signal is FLAT, so there is no waveform curvature for
a local fit to mis-track and no deterministic harmonic content to leak into the
residual and masquerade as code-dependent noise.

Because the whole staircase rides on the slowly decaying coupling, the same
nominal level lands at a slightly different code on each traverse. Segments are
therefore binned by their MEASURED mean code, not by step index, which pools
traverses correctly without assuming they align.

USAGE
-----
    python analyse_sweep.py C1_short_run1.txt --dwell-ms 25 \\
        --baseline-sd 0.9050 --json sweep_noise.json --plot Figure_21.png
    python analyse_sweep.py run1.txt run2.txt run3.txt --dwell-ms 25

Standard library plus numpy; matplotlib only if --plot is given.

Jack Davies, September 2026
"""
from __future__ import annotations
import argparse, json, math, sys
import numpy as np

LSB_V       = 3.3 / 4095.0
TAU_MS      = 533.7      # measured high-pass time constant of the chain
SETTLE_FRAC = 0.35       # discard this fraction at the START of each dwell
NBINS       = 16
MIN_PER_BIN = 300
METHOD_RATIO_FLOOR = 1.09   # what this estimator returns on genuinely flat noise
METHOD_BIAS        = 1.009  # and how much it overstates the pooled figure


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


def segment(codes, dwell_samples):
    """Split at step edges. The threshold is set from the data's own
    sample-to-sample spread, so it adapts to however far the staircase
    actually swings after the coupling."""
    sm = np.convolve(codes, np.ones(5) / 5, mode="same")
    d = np.abs(np.diff(sm))
    thr = max(8.0, 5.0 * float(np.median(np.abs(np.diff(codes)))))
    edges = np.flatnonzero(d > thr)
    bounds = [0] + [int(e) + 1 for e in edges] + [len(codes)]
    keep = int(0.5 * dwell_samples)
    return [(a, b) for a, b in zip(bounds[:-1], bounds[1:]) if b - a >= keep]


def analyse(path, dwell_ms, nbins=NBINS):
    codes, meta = read_capture(path)
    fs = float(meta.get("fs_hz", 2000.0))
    dwell_samples = int(round(dwell_ms * 1e-3 * fs))
    segs = segment(codes, dwell_samples)

    lv, rs = [], []
    for a, b in segs:
        k = int((b - a) * SETTLE_FRAC)
        core = codes[a + k:b]
        if len(core) < 12 or core.min() <= 0 or core.max() >= 4095:
            continue
        t = np.arange(len(core), dtype=float)
        p = np.polyfit(t, core, 1)
        lv.append(np.full(len(core), core.mean()))
        rs.append(core - np.polyval(p, t))
    if not lv:
        return {"file": path, "label": meta.get("label", "?"),
                "mode": meta.get("mode", "?"), "fs_hz": fs,
                "n_segments": len(segs), "bins": []}

    lv = np.concatenate(lv)
    rs = np.concatenate(rs)
    lo, hi = np.percentile(lv, [1, 99])
    edges = np.linspace(lo, hi, nbins + 1)
    idx = np.clip(np.digitize(lv, edges) - 1, 0, nbins - 1)

    rows = []
    for i in range(nbins):
        m = idx == i
        if m.sum() < MIN_PER_BIN:
            continue
        sd = float(rs[m].std(ddof=1))
        rows.append({"bin": i, "n": int(m.sum()),
                     "mean_code": float(lv[m].mean()),
                     "sd_code": sd, "sd_millivolts": sd * LSB_V * 1000})
    return {"file": path, "label": meta.get("label", "?"),
            "mode": meta.get("mode", "?"), "fs_hz": fs,
            "dwell_ms": dwell_ms,
            "decay_in_dwell_pct": 100 * (1 - math.exp(-dwell_ms / TAU_MS)),
            "n_segments": len(segs), "code_span": [float(lo), float(hi)],
            "bins": rows}


def report(res, baseline_sd=None):
    print("=" * 78)
    print(f"  {res['file']}   label={res['label']}   mode={res['mode']}"
          f"   fs={res['fs_hz']:.2f} Hz")
    if "dwell_ms" in res:
        print(f"  dwell {res['dwell_ms']} ms -> level falls "
              f"{res['decay_in_dwell_pct']:.1f} % across a dwell "
              f"(tau = {TAU_MS:.1f} ms)")
        if res["decay_in_dwell_pct"] > 15:
            print("  ** WARNING: that is too much. Shorten the dwell to 25 ms.")
    print(f"  {res['n_segments']} dwells segmented, "
          f"{len(res['bins'])} usable bins of {NBINS}")
    print("=" * 78)
    if not res["bins"]:
        print("  NO USABLE BINS. Either the staircase was not running, or the")
        print("  dwell is so long the levels have collapsed onto each other.")
        return None
    print(f"{'bin':>4}{'n':>8}{'mean code':>11}{'volts':>9}{'SD code':>9}{'SD mV':>8}")
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
    rho = float(np.corrcoef(np.argsort(np.argsort(means)),
                            np.argsort(np.argsort(sds)))[0, 1]) \
        if len(sds) > 2 else float("nan")

    print(f"  pooled across bins      {pooled:.4f} code "
          f"({pooled * LSB_V * 1000:.4f} mV)")
    print(f"  corrected for method    {pooled / METHOD_BIAS:.4f} code "
          f"({pooled / METHOD_BIAS * LSB_V * 1000:.4f} mV)")
    print(f"  range across bins       {lo:.4f} to {hi:.4f} code "
          f"(max/min = {ratio:.2f}x)")
    print(f"  rank correlation with code   rho = {rho:+.3f}")
    print(f"  method floor: on synthetic noise that is genuinely flat this")
    print(f"                estimator returns {METHOD_RATIO_FLOOR:.2f}x and a "
          f"{100 * (METHOD_BIAS - 1):.1f} % high pooled figure.")
    if baseline_sd is not None:
        print(f"  quiet-input baseline    {baseline_sd:.4f} code -> "
              f"corrected/baseline = {pooled / METHOD_BIAS / baseline_sd:.3f}x")
    print()
    # The ratio decides first. When the spread across bins is at the method's
    # own floor there is no dependence to have a direction, so the rank
    # correlation is reading the estimator's small residual trend rather than
    # the instrument, and it must not be allowed to raise a false alarm.
    if ratio < 1.20:
        print("  VERDICT: the noise floor is NOT materially code-dependent.")
        print("  The single quiet-input figure the error budget uses applies")
        print("  across the working range, and the last assumption in that")
        print("  budget becomes a measurement.")
    elif abs(rho) >= 0.50:
        print("  VERDICT: the noise floor IS code-dependent and trends with code.")
        print("  Quote the pooled figure in the budget and state which end of")
        print("  the range is worst.")
    else:
        print("  VERDICT: mixed - spread and trend disagree. Usually one bin is")
        print("  an outlier; read the table before writing anything.")
    print()
    return {"pooled_sd_code": pooled,
            "pooled_sd_code_corrected": pooled / METHOD_BIAS,
            "pooled_sd_millivolts_corrected": pooled / METHOD_BIAS * LSB_V * 1000,
            "min_sd_code": float(lo), "max_sd_code": float(hi),
            "ratio_max_min": float(ratio), "rank_corr_with_code": rho,
            "method_ratio_floor": METHOD_RATIO_FLOOR}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("captures", nargs="+")
    ap.add_argument("--dwell-ms", type=float, default=25.0,
                    help="dwell per level as flashed (default 25)")
    ap.add_argument("--bins", type=int, default=NBINS)
    ap.add_argument("--baseline-sd", type=float, default=None)
    ap.add_argument("--json", default=None)
    ap.add_argument("--plot", default=None)
    a = ap.parse_args(argv)

    out = {"captures": [], "summary": {}}
    for p in a.captures:
        res = analyse(p, a.dwell_ms, a.bins)
        res["summary"] = report(res, a.baseline_sd)
        out["captures"].append(res)

    good = [c for c in out["captures"] if c.get("summary")]
    if len(good) > 1:
        pooled = math.sqrt(float(np.mean(
            [c["summary"]["pooled_sd_code_corrected"] ** 2 for c in good])))
        ratios = [c["summary"]["ratio_max_min"] for c in good]
        print("=" * 78)
        print(f"  ACROSS {len(good)} CAPTURES: pooled {pooled:.4f} code "
              f"({pooled * LSB_V * 1000:.4f} mV)")
        print(f"  max/min per capture: {', '.join(f'{r:.2f}x' for r in ratios)}")
        print(f"  (the method's own floor is {METHOD_RATIO_FLOOR:.2f}x)")
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
        ax.set_title("Noise against converter code")
        ax.set_ylim(bottom=0)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(a.plot)
        print(f"  wrote {a.plot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
