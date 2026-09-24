#!/usr/bin/env python3
"""
verify_stimulus_residual_312.py - Section 3.12, paragraph "With the production
stimulus applied ..." (12.28 mV residual, 98.7 % harmonic, 1.57 % THD,
1.399 mV broadband, and the simulated 0.07 mHz / 0.14 mHz effect of the
distortion on the released estimator).

RECORD
  20260913/B1_ac_running.txt  (byte-identical to 20260923/captures/B1_ac_running.txt)
  deadline-scheduled 2 kHz characterisation loop, production 50 Hz stimulus.

METHOD
  1. Residual.  Exactly analyse_ac.py: fundamental refined on the first 4000
     samples, A*sin + B*cos + C fitted in each 200-sample window at that
     frequency, residual SD (ddof=3) per window, median over the 81 windows.
     (This is the released ac.json figure, sigma_ac_millivolts.)
  2. Harmonic share.  The 81 window residuals are concatenated (16,200
     samples) and a one-sided power spectrum formed (rectangular window,
     power-normalised so it sums to the residual mean square).  Power within
     +/-HALF_BW bins of every multiple h*f1, h = 1..19, is "harmonic"; the
     remainder is "broadband".  harmonic share = P_harm / P_total;
     broadband rms = sqrt(P_total - P_harm).  The result is reported for
     HALF_BW = 1..8 bins and for a Hann window as well, because the broadband
     figure depends (weakly) on that choice.
  3. THD.  rms of the harmonic residual (h = 2..19) over the rms of the fitted
     fundamental A/sqrt2; cross-checked against analyse_sinad.py (IEEE 1241,
     7-term Blackman-Harris, --harmonics 19) run on the same record.
  4. Simulation of the released estimator at the measured distortion.
     Harmonic amplitudes and phases (h = 2..19) are measured by a single
     least-squares multi-harmonic fit to the whole B1 record at the refined
     fundamental.  Two synthetic records are then built on the deadline
     loop's integer 500 us grid at the fitted fundamental frequency and
     amplitude, with the SAME Gaussian noise realisation (common random
     numbers, fixed seed) and quantised to integer codes:
         clean      = fundamental + noise
         distorted  = fundamental + measured harmonics + noise
     The exact port of measureFreq() (analyse_ac.estimate_freq: interpolated
     crossings truncated to integer us, both edges, median) is run on
     consecutive 200-sample windows, with the known DC bias as V2 would
     measure it.  Reported: change in dispersion (SD) and in mean, at the
     broadband noise (1.399 mV) and at the quiet-input noise (0.659 mV),
     with a bootstrap standard error on each difference.

Usage:  python3 verify_stimulus_residual_312.py [--windows 20000]
Deterministic (numpy default_rng, fixed seeds).  numpy only.
"""
from __future__ import annotations
import argparse, math, subprocess, sys, os
import numpy as np
import common_312 as c

HARM_MAX = 19          # 19 * 50 Hz = 950 Hz, the last harmonic below Nyquist


def window_residuals(v, fs, f1):
    res, sds, amps = [], [], []
    t = np.arange(c.FREQ_WIN) / fs
    for w in range(len(v) // c.FREQ_WIN):
        seg = v[w * c.FREQ_WIN:(w + 1) * c.FREQ_WIN]
        r, coef = c.fit_sine(seg, t, f1)
        res.append(r)
        sds.append(r.std(ddof=3))
        amps.append(math.hypot(coef[0], coef[1]))
    return np.concatenate(res), np.array(sds), float(np.median(amps))


def power_spectrum(r, fs, hann=False):
    n = len(r)
    w = np.hanning(n) if hann else np.ones(n)
    X = np.fft.rfft((r - r.mean()) * w)
    P = np.abs(X) ** 2 / (n * np.sum(w ** 2))
    P[1:] *= 2.0
    if n % 2 == 0:
        P[-1] /= 2.0
    return np.fft.rfftfreq(n, 1.0 / fs), P


def harmonic_mask(freqs, f1, half_bw, hmin=1, hmax=HARM_MAX):
    df = freqs[1] - freqs[0]
    m = np.zeros(len(freqs), bool)
    for h in range(hmin, hmax + 1):
        k = int(round(h * f1 / df))
        m[max(0, k - half_bw):min(len(m), k + half_bw + 1)] = True
    return m


def multiharmonic_fit(v, t, f1, hmax=HARM_MAX):
    cols = [np.ones_like(t)]
    for h in range(1, hmax + 1):
        cols += [np.sin(2 * np.pi * h * f1 * t), np.cos(2 * np.pi * h * f1 * t)]
    M = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(M, v, rcond=None)
    dc = coef[0]
    ab = coef[1:].reshape(-1, 2)
    amp = np.hypot(ab[:, 0], ab[:, 1])
    ph = np.arctan2(ab[:, 1], ab[:, 0])          # a*sin + b*cos = amp*sin(x+ph)
    return dc, amp, ph


def run_estimator(codes, bias_code, t_us):
    v = (codes - bias_code) * c.LSB_V
    out = []
    for w in range(len(codes) // c.FREQ_WIN):
        s = slice(w * c.FREQ_WIN, (w + 1) * c.FREQ_WIN)
        fe = c.measure_freq_median(v[s], t_us[s])
        out.append(np.nan if fe is None else fe)
    return np.array(out)


def simulate(amp, ph, f1, sigma_mv, nwin, seed, fs_grid_us=500):
    rng = np.random.default_rng(seed)
    n = nwin * c.FREQ_WIN
    t_us = (np.arange(n, dtype=np.int64) * fs_grid_us)
    t = t_us * 1e-6
    bias_code = 1870.0                                   # B1 mean code
    fund = amp[0] * np.sin(2 * np.pi * f1 * t + ph[0])
    harm = np.zeros(n)
    for h in range(2, len(amp) + 1):
        harm += amp[h - 1] * np.sin(2 * np.pi * h * f1 * t + ph[h - 1])
    noise = rng.normal(0.0, sigma_mv * 1e-3, n)
    clean = np.round((fund + noise) / c.LSB_V + bias_code)
    dist = np.round((fund + harm + noise) / c.LSB_V + bias_code)
    fc = run_estimator(clean, bias_code, t_us)
    fd = run_estimator(dist, bias_code, t_us)
    ok = np.isfinite(fc) & np.isfinite(fd)
    fc, fd = fc[ok], fd[ok]
    d_sd = (fd.std(ddof=1) - fc.std(ddof=1)) * 1e3
    d_mean = (fd.mean() - fc.mean()) * 1e3
    # bootstrap SE of the two differences (paired resampling of windows)
    br = np.random.default_rng(seed + 1)
    bs_sd, bs_mean = [], []
    for _ in range(400):
        i = br.integers(0, len(fc), len(fc))
        bs_sd.append((fd[i].std(ddof=1) - fc[i].std(ddof=1)) * 1e3)
        bs_mean.append((fd[i].mean() - fc[i].mean()) * 1e3)
    return {"sigma_mv": sigma_mv, "windows": int(len(fc)),
            "sd_clean_mhz": fc.std(ddof=1) * 1e3, "sd_distorted_mhz": fd.std(ddof=1) * 1e3,
            "mean_clean_offset_mhz": (fc.mean() - f1) * 1e3,
            "mean_distorted_offset_mhz": (fd.mean() - f1) * 1e3,
            "delta_sd_mhz": d_sd, "delta_sd_se_mhz": float(np.std(bs_sd)),
            "delta_mean_mhz": d_mean, "delta_mean_se_mhz": float(np.std(bs_mean))}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--windows", type=int, default=20000)
    a = ap.parse_args(argv)

    codes, meta = c.read_codes(c.B1)
    fs = float(meta["fs_hz"])
    v = codes.astype(float) * c.LSB_V
    bias = v.mean()
    f1 = c.refine_frequency(v[:4000] - bias, np.arange(4000) / fs, 50.0)
    r, sds, A = window_residuals(v, fs, f1)
    sigma_ac = float(np.median(sds))
    print("=" * 74)
    print("  B1_ac_running (13 Sep), production stimulus, 81 windows of 200")
    print("=" * 74)
    print(f"  fitted fundamental {f1:.5f} Hz, amplitude {A:.4f} V")
    print(f"  residual after single-tone fit (median of window SDs)  {sigma_ac*1e3:.3f} mV"
          f"   [paper 12.28]")

    table = []
    for hann in (False, True):
        fr, P = power_spectrum(r, fs, hann)
        for hb in (1, 2, 3, 5, 8):
            m = harmonic_mask(fr, f1, hb)
            share = P[m].sum() / P.sum()
            bb = math.sqrt(P[~m].sum())
            table.append({"window": "hann" if hann else "rect", "half_bw_bins": hb,
                          "harmonic_share_pct": 100 * share, "broadband_mv": bb * 1e3})
    print("\n  harmonic share / broadband rms of the concatenated window residuals")
    for row in table:
        print(f"    {row['window']:4s} +/-{row['half_bw_bins']} bins   "
              f"{row['harmonic_share_pct']:6.3f} %   {row['broadband_mv']:.4f} mV")
    ref = next(x for x in table if x["window"] == "rect" and x["half_bw_bins"] == 5)
    bb_lo = min(x["broadband_mv"] for x in table)
    bb_hi = max(x["broadband_mv"] for x in table)
    print(f"  -> reference choice (rect, +/-5 bins): {ref['harmonic_share_pct']:.2f} %, "
          f"{ref['broadband_mv']:.3f} mV   [paper 98.7 %, 1.399 mV]")
    print(f"     sensitivity over the choices above: broadband {bb_lo:.3f} to {bb_hi:.3f} mV")

    # THD from the residual
    fr, P = power_spectrum(r, fs, False)
    mh = harmonic_mask(fr, f1, 5, hmin=2)
    thd_resid = math.sqrt(P[mh].sum()) / (A / math.sqrt(2)) * 100
    sinad = subprocess.run([sys.executable, os.path.join(c.HERE, "analyse_sinad.py"),
                            c.B1, "--harmonics", str(HARM_MAX)],
                           capture_output=True, text=True).stdout
    thd_sinad = float([l for l in sinad.splitlines() if l.strip().startswith("THD")][0]
                      .split("(")[1].split("%")[0])
    print(f"\n  THD, harmonic residual (h=2..19) / fundamental rms   {thd_resid:.3f} %")
    print(f"  THD, analyse_sinad.py --harmonics 19 (IEEE 1241)     {thd_sinad:.3f} %"
          f"   [paper 1.57 %]")
    print(f"  broadband / quiet-input 0.659 mV                     "
          f"{ref['broadband_mv']/0.659:.2f}x   [paper 'roughly twice']")

    # distortion model from the whole record
    t = np.arange(len(v)) / fs
    fg = c.refine_frequency(v - bias, t, 50.0)
    dc, amp, ph = multiharmonic_fit(v - bias, t, fg)
    print(f"\n  multi-harmonic fit, whole record at {fg:.5f} Hz: A1 = {amp[0]*1e3:.1f} mV;"
          f" harmonic rms {math.sqrt(np.sum(amp[1:]**2)/2)*1e3:.2f} mV"
          f" (THD {math.sqrt(np.sum(amp[1:]**2))/amp[0]*100:.3f} %)")

    sims = []
    print(f"\n  SIMULATION, exact port of measureFreq(), {a.windows} windows, "
          f"common random numbers")
    for k, s in enumerate((1.399, 0.659)):
        res = simulate(amp, ph, fg, s, a.windows, seed=312 + k)
        sims.append(res)
        print(f"    noise {s:.3f} mV: SD clean {res['sd_clean_mhz']:.3f} -> distorted "
              f"{res['sd_distorted_mhz']:.3f} mHz, change {res['delta_sd_mhz']:+.3f} "
              f"+/- {res['delta_sd_se_mhz']:.3f};  mean change {res['delta_mean_mhz']:+.3f} "
              f"+/- {res['delta_mean_se_mhz']:.3f} mHz")
    print("    [paper: dispersion changed by 0.07 mHz, mean by 0.14 mHz]")

    c.save_json("stimulus_residual_312.json", {
        "record": os.path.relpath(c.B1, c.HERE), "f1_hz": f1, "amplitude_v": A,
        "residual_mv": sigma_ac * 1e3, "harmonic_table": table,
        "reference_choice": ref, "thd_pct_from_residual": thd_resid,
        "thd_pct_analyse_sinad_h19": thd_sinad,
        "harmonic_amplitudes_mv": (amp * 1e3).tolist(),
        "simulation": sims})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
