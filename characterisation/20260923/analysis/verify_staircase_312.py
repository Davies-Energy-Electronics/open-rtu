#!/usr/bin/env python3
"""
verify_staircase_312.py - Section 3.12, the 25 ms staircase: "at the bias
point the staircase returns 1.00 codes against the 0.905 of the static quiet
input ... Away from the bias point it returns 1.7 to 2.3 codes", "the measured
spread correlates at +0.64 with distance from the bias point ... and at +0.03
with code", and the "factor of about 2.4" upper bound.

INPUT
  20260923/captures/sweep_noise.json - the released output of analyse_sweep.py
  on C1_short_run1..3.txt.  The per-bin table is re-derived from the raw
  captures with analyse_sweep.analyse() and checked against the JSON first.
  0.905 codes is the quadrature mean of the three WARMED single-channel
  quiet-input records (see verify_quiet_noise_312.py).

METHOD
  All 45 (level, SD) bins of the three records are pooled.  Distance from the
  bias point = |mean code - 1872.2| (the quiet-input mean code).  Pearson r of
  SD against distance and against code (Spearman given alongside).  The bias
  point is the bin whose mean code is within 5 codes of 1872; "away from the
  bias point" excludes that bin and the next one up (2064, the second
  stationary level named in the text), and also shown excluding only the first.

Usage:  python3 verify_staircase_312.py
"""
from __future__ import annotations
import json, math, os
import numpy as np
from scipy import stats
import common_312 as c
import analyse_sweep as sw

BIAS_CODE = 1872.2


def main():
    d = json.load(open(os.path.join(c.CAP_0923, "sweep_noise.json")))
    # re-derive from raw captures and compare
    maxdiff = 0.0
    for cap in d["captures"]:
        fn = os.path.basename(cap["file"].replace("\\", "/"))
        res = sw.analyse(os.path.join(c.CAP_0923, fn), 25.0)
        for a, b in zip(res["bins"], cap["bins"]):
            maxdiff = max(maxdiff, abs(a["sd_code"] - b["sd_code"]))
    print(f"  raw captures re-analysed with analyse_sweep.py: max |dSD| vs sweep_noise.json "
          f"= {maxdiff:.2e} codes")

    M, S, R = [], [], []
    for cap in d["captures"]:
        for b in cap["bins"]:
            M.append(b["mean_code"]); S.append(b["sd_code"]); R.append(cap["label"])
    M, S = np.array(M), np.array(S)
    dist = np.abs(M - BIAS_CODE)
    r_d, p_d = stats.pearsonr(dist, S)
    r_c, p_c = stats.pearsonr(M, S)
    s_d = stats.spearmanr(dist, S)[0]
    s_c = stats.spearmanr(M, S)[0]
    print(f"  {len(S)} bins: Pearson r(SD, distance) = {r_d:+.3f} (p={p_d:.1e}); "
          f"r(SD, code) = {r_c:+.3f} (p={p_c:.2f})   [paper +0.64, +0.03]")
    print(f"           Spearman   (SD, distance) = {s_d:+.3f};  (SD, code) = {s_c:+.3f}")
    print("  per record (Pearson dist / code): " + "; ".join(
        f"{cap['label']} {stats.pearsonr(np.abs(np.array([b['mean_code'] for b in cap['bins']])-BIAS_CODE), [b['sd_code'] for b in cap['bins']])[0]:+.2f}"
        f" / {stats.pearsonr([b['mean_code'] for b in cap['bins']], [b['sd_code'] for b in cap['bins']])[0]:+.2f}"
        for cap in d["captures"]))

    at_bias = S[np.abs(M - 1872) < 5]
    second = S[np.abs(M - 2064) < 5]
    away2 = S[(np.abs(M - 1872) >= 5) & (np.abs(M - 2064) >= 5)]
    away1 = S[np.abs(M - 1872) >= 5]
    print(f"  bias-point bin (1872): {np.round(at_bias,3)} -> mean {at_bias.mean():.3f}, "
          f"rms {c.rms(at_bias):.3f} codes   [paper 1.00]")
    print(f"  second stationary bin (2064): {np.round(second,3)}")
    print(f"  away from both: {away2.min():.3f} to {away2.max():.3f} codes   [paper 1.7 to 2.3]")
    print(f"  away from 1872 only: {away1.min():.3f} to {away1.max():.3f} codes")
    two_lowest = sorted({int(round(x)) for x in M[np.argsort(S)[:6]]})
    print(f"  bins holding the six lowest SDs (codes): {two_lowest}")
    ratios = d["summary"]["ratio_max_min_per_capture"]
    print(f"  max/min per record {', '.join(f'{x:.2f}' for x in ratios)}   [paper 'about 2.4']")
    print(f"  code span of bin means {M.min():.0f} to {M.max():.0f}")
    c.save_json("staircase_312.json", {
        "pearson_sd_vs_distance": r_d, "pearson_sd_vs_code": r_c,
        "spearman_sd_vs_distance": s_d, "spearman_sd_vs_code": s_c,
        "bias_point_sd_codes": at_bias.tolist(), "bias_point_mean": float(at_bias.mean()),
        "away_range_codes": [float(away2.min()), float(away2.max())],
        "ratio_max_min_per_record": ratios, "max_abs_diff_vs_json": maxdiff})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
