#!/usr/bin/env python3
"""
analyse_ac.py — noise and estimator performance under real operating conditions.

Reads a characterise.ino capture taken with Feather 2 running its normal
50 Hz output, and answers two questions the quiet-input capture cannot:

  1. What is the input-referred noise WITH THE STIMULUS PRESENT? A sinusoid
     is least-squares fitted to each 200-sample window and the residual taken.
     If that residual is much larger than the quiet-input figure, the extra
     is coming from the source, not the measurement chain — Feather 2's 8-bit
     DAC puts one code at roughly 10.6 mV referred to the ADC input, which is
     an order of magnitude above the measured converter noise.

  2. What frequency dispersion does the estimator actually achieve on this
     bench today? The exact firmware algorithm — two-point interpolation,
     same-polarity periods from both edges, median — is run over consecutive
     200-sample windows and the scatter reported. That is directly comparable
     to the realised floor quoted in Section 2.3 (5.46 mHz with n - 1
     normalisation, 5.40 mHz with population normalisation), reproduced from
     raw samples rather than inherited from a 2026-06-27 capture.

Usage:
    python analyse_ac.py B1_ac_running.txt
    python analyse_ac.py B1_ac_running.txt --quiet-sigma-mv 0.687 --json ac.json

Standard library plus numpy.
"""
from __future__ import annotations
import argparse, json, math, sys
import numpy as np

LSB_V     = 3.3 / 4095.0
SIGMA_Q_V = LSB_V / math.sqrt(12)
FREQ_WIN  = 200
F_NOM     = 50.0


def read_capture(path):
    codes, meta = [], {}
    for line in open(path, errors="replace"):
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            b = s.lstrip("#").strip()
            if "," in b:
                k, _, v = b.partition(",")
                meta[k.strip()] = v.strip()
            continue
        if s.lstrip("-").isdigit():
            codes.append(int(s))
    if not codes:
        sys.exit(f"{path}: no ADC codes found")
    return np.asarray(codes, dtype=float), meta


def fit_sine(v, t, f):
    """Least squares for A*sin + B*cos + C at known f. Returns residual."""
    M = np.column_stack([np.sin(2 * math.pi * f * t),
                         np.cos(2 * math.pi * f * t),
                         np.ones_like(t)])
    coef, *_ = np.linalg.lstsq(M, v, rcond=None)
    return v - M @ coef, coef


def refine_frequency(v, t, f0):
    """Coarse scan then golden refine on residual energy. The generator is
    not at exactly 50 Hz and a fit at the wrong f leaks signal into the
    residual, which would masquerade as noise."""
    def cost(f):
        r, _ = fit_sine(v, t, f)
        return float(r @ r)
    grid = np.linspace(f0 - 1.0, f0 + 1.0, 201)
    f = min(grid, key=cost)
    lo, hi = f - 0.01, f + 0.01
    for _ in range(60):
        m1 = lo + (hi - lo) / 3
        m2 = hi - (hi - lo) / 3
        if cost(m1) < cost(m2):
            hi = m2
        else:
            lo = m1
    return 0.5 * (lo + hi)


def estimate_freq(v_volts, t_us):
    """Exact port of measureFreq()/ringMeasureFreqVb(): two-point linear
    interpolation truncated to integer microseconds, same-polarity periods
    from both edges combined, median."""
    pos, neg = [], []
    for i in range(1, len(v_volts)):
        pv, cv = v_volts[i - 1], v_volts[i]
        pt, ct = t_us[i - 1], t_us[i]
        if pv < 0.0 and cv >= 0.0:
            dv = cv - pv
            frac = (-pv / dv) if dv > 1e-6 else 0.5
            pos.append(pt + int(min(max(frac, 0.0), 1.0) * (ct - pt)))
        elif pv >= 0.0 and cv < 0.0:
            dv = pv - cv
            frac = (pv / dv) if dv > 1e-6 else 0.5
            neg.append(pt + int(min(max(frac, 0.0), 1.0) * (ct - pt)))
    per = ([pos[i + 1] - pos[i] for i in range(len(pos) - 1)] +
           [neg[i + 1] - neg[i] for i in range(len(neg) - 1)])
    if len(per) < 3:
        return None
    med = float(np.median(per))
    if med < 100.0:
        return None
    return 1e6 / med


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("capture")
    ap.add_argument("--quiet-sigma-mv", type=float, default=0.687,
                    help="sigma_v measured on the quiet-input capture (default 0.687)")
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)

    codes, meta = read_capture(a.capture)
    fs = float(meta.get("fs_hz", 2000.0))
    n = len(codes)
    v_all = codes * LSB_V
    bias = float(np.mean(v_all))

    nwin = n // FREQ_WIN
    if nwin < 4:
        sys.exit("  Capture too short for windowed analysis.")

    # Global frequency first — the generator is not exactly 50 Hz.
    t0 = np.arange(min(4000, n)) / fs
    f_hat = refine_frequency(v_all[:len(t0)] - bias, t0, F_NOM)

    resids, amps, freqs = [], [], []
    for w in range(nwin):
        sl = slice(w * FREQ_WIN, (w + 1) * FREQ_WIN)
        seg = v_all[sl]
        t = np.arange(FREQ_WIN) / fs
        r, coef = fit_sine(seg, t, f_hat)
        resids.append(r.std(ddof=3))
        amps.append(math.hypot(coef[0], coef[1]))
        t_us = np.round(t * 1e6).astype(np.int64)
        fe = estimate_freq(seg - bias, t_us)
        if fe is not None and 45.0 < fe < 55.0:
            freqs.append(fe)

    resids = np.array(resids)
    amps = np.array(amps)
    freqs = np.array(freqs)

    sigma_ac = float(np.median(resids))
    sigma_quiet = a.quiet_sigma_mv * 1e-3
    extra = math.sqrt(max(sigma_ac ** 2 - sigma_quiet ** 2, 0.0))
    A = float(np.median(amps))

    print(f"\n{'='*68}")
    print(f"  OPERATING-CONDITION ANALYSIS — {meta.get('label', a.capture)}")
    print(f"{'='*68}")
    print(f"  samples {n}   fs {fs:.2f} Hz   windows {nwin} x {FREQ_WIN}")
    print(f"  fitted fundamental        {f_hat:.4f} Hz")
    print(f"  fitted amplitude          {A:.4f} V "
          f"({A/LSB_V:.0f} codes)   bias {bias:.4f} V")
    print(f"  {'-'*64}")
    print(f"  residual after sine fit   {sigma_ac*1e3:.4f} mV"
          f"   = {sigma_ac/LSB_V:.2f} LSB")
    print(f"  quiet-input reference     {sigma_quiet*1e3:.4f} mV")
    print(f"  in quadrature, the extra  {extra*1e3:.4f} mV"
          f"   <-- attributable to the stimulus and signal-dependent terms")
    print(f"  {'-'*64}")
    if len(freqs) >= 4:
        print(f"  ESTIMATOR, run on these samples with the firmware algorithm")
        print(f"    windows yielding an estimate   {len(freqs)} of {nwin}")
        print(f"    mean                           {np.mean(freqs):.5f} Hz")
        print(f"    dispersion (sd)                {np.std(freqs, ddof=1)*1e3:.3f} mHz")
        d = np.diff(freqs)
        print(f"    successive-difference floor    "
              f"{np.std(d, ddof=1)/math.sqrt(2)*1e3:.3f} mHz"
              f"   <-- same statistic as Sec 2.3's 5.46 mHz floor")
    else:
        print("  Estimator produced too few valid windows — is the AC running?")
    print(f"{'='*68}")

    # Interpretation
    r = sigma_ac / sigma_quiet if sigma_quiet else float("inf")
    if r > 2.5:
        print("  >> THE STIMULUS DOMINATES. Noise with the signal present is")
        print(f"     {r:.1f}x the quiet-input figure. The realised floor is a")
        print("     property of the 8-bit DAC source, not of the acquisition")
        print("     chain. The floor cannot be attributed to the chain alone,")
        print("     and that favours the instrument (see Section 3.12).")
    elif r > 1.4:
        print("  >> PARTIALLY. The stimulus contributes but does not dominate.")
        print("     Report both terms; neither alone explains the floor.")
    else:
        print("  >> THE STIMULUS IS NOT THE EXPLANATION. Noise with the signal")
        print("     present is close to the quiet figure, so the gap between")
        print("     1.75 mHz predicted and 5.46 mHz realised lies elsewhere.")
        print("     Do NOT edit the paper until it is found.")
    print()

    if a.json:
        out = {"label": meta.get("label", a.capture), "fs_hz": fs,
               "f_fitted_hz": f_hat, "amplitude_volts": A, "bias_volts": bias,
               "sigma_ac_millivolts": sigma_ac * 1e3,
               "sigma_quiet_millivolts": a.quiet_sigma_mv,
               "stimulus_term_millivolts": extra * 1e3,
               "n_windows": nwin, "n_estimates": len(freqs),
               "freq_mean_hz": float(np.mean(freqs)) if len(freqs) else None,
               "freq_sd_mhz": float(np.std(freqs, ddof=1) * 1e3) if len(freqs) > 1 else None}
        with open(a.json, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"  -> {a.json}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
