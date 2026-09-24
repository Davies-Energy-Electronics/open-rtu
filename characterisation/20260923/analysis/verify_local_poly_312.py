#!/usr/bin/env python3
"""
verify_local_poly_312.py - Section 3.12, the rejected alternating-current
method: "a local-polynomial estimator that recovers a known noise figure from
synthetic records to better than one per cent returns four times the
quiet-input noise on a real record, and a spectral decomposition of its
residual places the excess at 550, 650, 750, 850 and 950 Hz".

ESTIMATOR
  The released one, unchanged: analyse_ac_noise_vs_code.py (Savitzky-Golay,
  9-sample window, order 4, centre-tap residual divided by the analytic
  residual scale, 16 code bins, pooled SD).

METHOD
  1. Synthetic recovery (fixed seed).  Gaussian noise of SD 0.905 codes
     is added to (a) a constant level, (b) a clean 50 Hz sine of the measured
     amplitude (1356 codes), (c) that sine plus a third harmonic at the level
     measured on B1 (-38.8 dBc), each at 2 kHz, 16,384 samples, not quantised
     and then quantised to integer codes; 20 records per case.  The released
     local_residuals() + binned pooled SD is applied and the error reported
     against the SD of (record - noiseless signal), which for a quantised
     record includes the quantisation error.
  2. Real records.  The same estimator on B1_ac_running.txt (13 Sep) and on
     the quiet A1_quiet_single.txt (--calibrate path).  Ratio of the pooled
     figure (raw and with the released 1.045 leakage correction) to the
     quiet-input 0.659 mV (0.818 codes) and to the warmed 0.905 codes.
  3. Spectral decomposition.  Hann-windowed power spectrum of the B1 residual
     (divided by the residual scale), power within +/-4 bins of each multiple
     of the fitted fundamental; the five largest are listed with their share
     of the residual power.

Usage:  python3 verify_local_poly_312.py
"""
from __future__ import annotations
import math
import numpy as np
import common_312 as c
import analyse_ac_noise_vs_code as sg

QUIET_CODE_0659 = 0.659e-3 / c.LSB_V
WARM_CODE = 0.905


def pooled(codes):
    """The released binned, pooled SD (analyse() + report()'s pooling)."""
    level, resid, scale = sg.local_residuals(codes)
    lo, hi = np.percentile(level, [0.5, 99.5])
    edges = np.linspace(lo, hi, sg.NBINS + 1)
    idx = np.clip(np.digitize(level, edges) - 1, 0, sg.NBINS - 1)
    num = den = 0.0
    for b in range(sg.NBINS):
        m = idx == b
        if m.sum() < sg.MIN_PER_BIN:
            continue
        sd = resid[m].std(ddof=1) / scale
        num += (m.sum() - 1) * sd ** 2
        den += m.sum() - 1
    if den == 0:                                  # constant level: one population
        return float(resid.std(ddof=1) / scale)
    return math.sqrt(num / den)


def main():
    rng = np.random.default_rng(3120)
    n, fs, true = 16384, 2000.0, WARM_CODE
    t = np.arange(n) / fs
    A = 1356.0
    h3 = A * 10 ** (-38.77 / 20)
    cases = {"constant level": np.zeros(n),
             "50 Hz sine": A * np.sin(2 * np.pi * 50 * t),
             "50 Hz sine + H3 (-38.8 dBc)": A * np.sin(2 * np.pi * 50 * t)
             + h3 * np.sin(2 * np.pi * 150 * t + 0.7)}
    print("=" * 78 + "\n  1. SYNTHETIC RECOVERY, true SD 0.905 codes, released SG(9,4) estimator\n" + "=" * 78)
    synth = []
    NREP = 20
    for name, sig in cases.items():
        for q in (False, True):
            errs = []
            for _ in range(NREP):
                x = sig + 1872.0 + rng.normal(0, true, n)
                if q:
                    x = np.round(x)
                # the known figure: SD of everything that is not the signal
                # (for a quantised record that includes the quantisation error)
                ref = float(np.std(x - sig - 1872.0, ddof=1))
                errs.append(100 * (pooled(x) / ref - 1))
            errs = np.array(errs)
            synth.append({"case": name, "quantised": q, "n_records": NREP,
                          "mean_error_pct": float(errs.mean()),
                          "sd_error_pct": float(errs.std(ddof=1)),
                          "max_abs_error_pct": float(np.abs(errs).max())})
            print(f"  {name:30s} {'quantised' if q else 'continuous':10s} "
                  f"mean error {errs.mean():+.2f} %, SD {errs.std(ddof=1):.2f} %, "
                  f"max |error| {np.abs(errs).max():.2f} %  ({NREP} records)")
    print("=" * 78 + "\n  2. REAL RECORDS\n" + "=" * 78)
    qcodes, _ = c.read_codes(c.SEPT_QUIET[0])
    qcodes = qcodes.astype(float)
    lvl, res, scale = sg.local_residuals(qcodes)
    qest = res.std(ddof=1) / scale
    qtrue = qcodes.std(ddof=1)
    print(f"  A1_quiet_single: estimator {qest:.4f} vs whole-record SD {qtrue:.4f} codes "
          f"({100*(qest/qtrue-1):+.2f} %)")
    bcodes, meta = c.read_codes(c.B1)
    bcodes = bcodes.astype(float)
    res_b = sg.analyse(c.B1)
    sds = np.array([r["sd_code"] for r in res_b["bins"]])
    ns = np.array([r["n"] for r in res_b["bins"]], float)
    pb = math.sqrt(float(np.sum((ns - 1) * sds ** 2) / np.sum(ns - 1)))
    pbc = pb / sg.LEAK_BIAS
    print(f"  B1_ac_running pooled {pb:.4f} codes ({pb*c.LSB_V*1e3:.3f} mV); "
          f"leakage-corrected {pbc:.4f} codes")
    print(f"    / quiet 0.659 mV ({QUIET_CODE_0659:.4f} codes): raw {pb/QUIET_CODE_0659:.2f}x, "
          f"corrected {pbc/QUIET_CODE_0659:.2f}x   [paper 'four times']")
    print(f"    / warmed 0.905 codes:                  raw {pb/WARM_CODE:.2f}x, "
          f"corrected {pbc/WARM_CODE:.2f}x")

    print("=" * 78 + "\n  3. SPECTRAL DECOMPOSITION OF THE B1 RESIDUAL\n" + "=" * 78)
    fs_b = float(meta["fs_hz"])
    v = bcodes * c.LSB_V
    f1 = c.refine_frequency(v[:4000] - v.mean(), np.arange(4000) / fs_b, 50.0)
    lvl, r, scale = sg.local_residuals(bcodes)
    r = r / scale
    N = len(r)
    w = np.hanning(N)
    P = np.abs(np.fft.rfft((r - r.mean()) * w)) ** 2 / (N * np.sum(w ** 2))
    P[1:] *= 2
    df = fs_b / N
    tot = P.sum()
    rows = []
    for h in range(1, 20):
        k = int(round(h * f1 / df))
        rows.append((h, h * f1, float(P[max(0, k - 4):k + 5].sum())))
    hp = sum(p for *_, p in rows)
    top = sorted(rows, key=lambda x: -x[2])[:5]
    for h, f, p in sorted(top, key=lambda x: x[1]):
        print(f"  H{h:<3d} {f:7.1f} Hz  {100*p/tot:5.1f} % of residual power")
    print(f"  top five together {100*sum(p for *_, p in top)/tot:.1f} %; all harmonics "
          f"{100*hp/tot:.1f} %; residual rms {math.sqrt(tot):.3f} codes, non-harmonic "
          f"{math.sqrt(tot-hp):.3f} codes")
    c.save_json("local_poly_312.json", {
        "synthetic": synth, "quiet_A1_estimate_codes": qest, "quiet_A1_true_codes": qtrue,
        "B1_pooled_codes": pb, "B1_pooled_corrected_codes": pbc,
        "ratio_to_0659_raw": pb / QUIET_CODE_0659, "ratio_to_0659_corrected": pbc / QUIET_CODE_0659,
        "ratio_to_0905_raw": pb / WARM_CODE, "ratio_to_0905_corrected": pbc / WARM_CODE,
        "top5_harmonics": [{"h": h, "f_hz": f, "share_pct": 100 * p / tot} for h, f, p in top],
        "harmonic_share_pct": 100 * hp / tot})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
