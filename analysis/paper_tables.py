#!/usr/bin/env python3
"""
paper_tables.py
===============

Regenerates, from the released captures, every table value in Section 3 of
the manuscript that no other released script produces, and checks each one
against the figure printed in the paper.

Covered here (manuscript numbering as submitted):

    Section 2.3   realised noise floor of the canonical capture
                  (successive-difference statistic, n - 1 normalisation)
    Section 3.1   per-region statistics of the canonical capture on the
                  full-stream basis (n = 379)
    Table 2       TOST with exact Student-t p-values, 90 % t intervals,
                  95 % percentile-bootstrap intervals, and the serial-
                  correlation (effective sample size) sensitivity check
    Table 3       TVE-equivalent statistics for all five regions
                  (read from the tve_metrics.py sidecar)
    Table 4       stationary-/dynamic-region pooled dispersions and the
                  ratios to the Cramer-Rao bound
    Table 5       per-region mean / sd for V1-V4 on the index-aligned basis
    Table 7       three-phase captures on the Section 3.1 partition
    Table 8       V6 parallel-firmware capture, per region
    Section 3.8   run-to-run repeatability of the three 24 July 2026 V6
                  captures on the Section 3.1 partition

Every region below is the Section 3.1 partition of the 9 August 2019
trajectory, by replay index:

    Pre-event nominal        0 - 45
    RoCoF descent           46 - 119
    Nadir + early recovery 120 - 200
    Late recovery          201 - 300
    Post-event settled     301 - 359

Usage (from the repository root or from analysis/):
    python analysis/paper_tables.py            # print tables, write JSON
    python analysis/paper_tables.py --check    # also verify against the paper

Writes paper_tables_output.json next to this script's repository root.
Requires numpy and scipy (both already required by reproduce_all.py).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
import sys

import numpy as np
from scipy import stats as sps

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE) if os.path.isdir(os.path.join(os.path.dirname(HERE), "analysis")) else HERE
sys.path.insert(0, HERE)

import hybrid_freq_analysis as H          # noqa: E402  (released V3 emulation)
import kalman_freq as K                   # noqa: E402  (released V4 filter)
import analyse_v6_repeatability as R      # noqa: E402  (released V6 parser/aligner)
from tost_metrics import bootstrap_ci_mean  # noqa: E402  (released bootstrap, seed 20260628)


def find(name):
    for d in ("captures", "trajectories", "", "analysis"):
        p = os.path.join(ROOT, d, name)
        if os.path.isfile(p):
            return p
    raise FileNotFoundError(name)


CANONICAL = "replay_20260627_102507V2.csv"
V1_CSV = "replay_20260625_180520V1.csv"
V6_CSV = "replay_20260716_152206_parallel.csv"
BAL_CSV = "replay_20260714_164522_3phase_balanced.csv"
THD_CSV = "replay_20260714_174352_3phase_thd_varied.csv"
REPEAT = ["replay_20260724_130842_v6_run1.csv",
          "replay_20260724_131622_v6_run2.csv",
          "replay_20260724_132305_v6_run3.csv"]
TRAJ = "neso_trajectory.csv"
TVE_JSON = "replay_20260627_102507V2_tve.json"

REGIONS = [("Pre-event nominal", 0, 45),
           ("RoCoF descent", 46, 119),
           ("Nadir + early recovery", 120, 200),
           ("Late recovery", 201, 300),
           ("Post-event settled", 301, 359)]
STATIONARY = ("Pre-event nominal", "Post-event settled")
DYNAMIC = ("RoCoF descent", "Nadir + early recovery", "Late recovery")
ALPHA = 0.05
CRB_MHZ = 3.5408423594085563   # crb_analysis.py, conservative 40 dB, N = 200, 1816.5 Hz


# ---------------------------------------------------------------------------
# loaders
# ---------------------------------------------------------------------------

def load_scada(path):
    """(replay_index, f_Vb, f_target) for every replay-tagged row with a target."""
    out = []
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            try:
                i = int(r["replay_index"])
                if i < 0 or r["f_target_hz"] in ("", None):
                    continue
                out.append((i, float(r["f_Vb"]), float(r["f_target_hz"])))
            except (KeyError, ValueError):
                continue
    return out


def region_of(i):
    for name, a, b in REGIONS:
        if a <= i <= b:
            return name
    return None


# ---------------------------------------------------------------------------
# Table 2 - TOST
# ---------------------------------------------------------------------------

def tost_exact(e, delta, df=None):
    n = len(e)
    m = float(np.mean(e))
    s = float(np.std(e, ddof=1))
    if df is None:
        df = n - 1
        se = s / math.sqrt(n)
    else:
        se = s / math.sqrt(df + 1)
    t1 = (m + delta) / se
    t2 = (m - delta) / se
    p1 = float(sps.t.sf(t1, df))
    p2 = float(sps.t.cdf(t2, df))
    return max(p1, p2)


def ci90(e, n_eff=None):
    n = len(e)
    m = float(np.mean(e))
    s = float(np.std(e, ddof=1))
    neff = n if n_eff is None else n_eff
    se = s / math.sqrt(neff)
    tc = float(sps.t.ppf(1 - ALPHA, neff - 1))
    return m - tc * se, m + tc * se


def table2(rows):
    wins = [("Pre-event nominal (i = 0-45)", [e for i, e in rows if 0 <= i <= 45]),
            ("Post-event settled (i = 301-359)", [e for i, e in rows if 301 <= i <= 359]),
            ("Full replay window (i = 0-359)", [e for i, e in rows])]
    out = []
    for name, e in wins:
        e = np.array(e)
        n = len(e)
        lo, hi = ci90(e)
        boot = bootstrap_ci_mean(list(e))                       # released: 95 %, seed 20260628
        x = e - e.mean()
        r1 = float(np.sum(x[1:] * x[:-1]) / np.sum(x * x))
        neff = n * (1 - r1) / (1 + r1)
        lo_e, hi_e = ci90(e, neff)
        out.append({
            "window": name, "n": n,
            "mean_mhz": float(e.mean()), "sd_mhz": float(e.std(ddof=1)),
            "ci90_t_mhz": [lo, hi],
            "ci95_bootstrap_mhz": [boot["bootstrap_ci_lower_mhz"], boot["bootstrap_ci_upper_mhz"]],
            "p_tost_delta10": tost_exact(e, 10.0),
            "p_tost_delta100": tost_exact(e, 100.0),
            "crossover_continuous_mhz": max(abs(lo), abs(hi)),
            "lag1_autocorrelation": r1,
            "n_eff": neff,
            "crossover_neff_mhz": max(abs(lo_e), abs(hi_e)),
            "p_tost_delta100_neff": tost_exact(e, 100.0, df=neff - 1),
        })
    return out


# ---------------------------------------------------------------------------
# Section 2.3 / 3.1 - canonical capture, full-stream basis
# ---------------------------------------------------------------------------

def canonical_stats(sc):
    per = {}
    for name, a, b in REGIONS:
        e = np.array([(f - t) * 1000 for i, f, t in sc if a <= i <= b])
        per[name] = {"n": len(e), "mean_mhz": float(e.mean()), "sd_mhz": float(e.std(ddof=1)),
                     "max_abs_mhz": float(np.abs(e).max())}
    f_pre = np.array([f for i, f, t in sc if 0 <= i <= 45]) * 1000
    f_set = np.array([f for i, f, t in sc if 301 <= i <= 359]) * 1000
    floor = float(np.std(np.diff(f_pre), ddof=1) / math.sqrt(2))
    floor_pop = float(np.std(np.diff(f_pre), ddof=0) / math.sqrt(2))
    tail = float(np.std(np.diff(f_set), ddof=1) / math.sqrt(2))
    ref_pre = np.array([t for i, f, t in sc if 0 <= i <= 45]) * 1000
    e_all = np.array([abs(f - t) * 1000 for i, f, t in sc])
    return {"regions": per, "noise_floor_pre_event_mhz": floor,
            "noise_floor_pre_event_population_norm_mhz": floor_pop,
            "noise_floor_settled_mhz": tail,
            "reference_sd_pre_event_mhz": float(np.std(ref_pre, ddof=1)),
            "coverage": {"tier2_pct": float(100 * np.mean(e_all <= 100)),
                         "tier1_pct": float(100 * np.mean(e_all <= 10))}}


# ---------------------------------------------------------------------------
# Tables 4 and 5 - four estimators on the index-aligned basis
# ---------------------------------------------------------------------------

def four_estimators():
    v1 = H.load_frames(find(V1_CSV))
    v2 = H.load_frames(find(CANONICAL))
    aligned = H.align_v1_v2(v1, v2)
    hyb, _ = H.emulate_state_machine(aligned)
    rows = K.load_capture(find(CANONICAL))
    neso = K.load_neso(find(TRAJ))
    filt = K.kalman_run(rows, K.derive_Q(rows, neso), K.derive_R(rows))
    v4 = {}
    for d, fv in zip(rows, filt):
        v4.setdefault(d["i"], fv)
    series = {"V1": [], "V2": [], "V3": [], "V4": []}
    for o in hyb:
        i, t = o["replay_idx"], o["f_target_hz"]
        series["V1"].append((i, (o["f_v1"] - t) * 1000))
        series["V2"].append((i, (o["f_v2"] - t) * 1000))
        series["V3"].append((i, (o["f_hybrid"] - t) * 1000))
        series["V4"].append((i, (v4[i] - t) * 1000))
    t5 = {}
    for alg, s in series.items():
        t5[alg] = {}
        for name, a, b in REGIONS + [("Aggregate (aligned)", 0, 359)]:
            e = np.array([x for i, x in s if a <= i <= b])
            t5[alg][name] = {"n": len(e), "mean_mhz": float(e.mean()), "sd_mhz": float(e.std(ddof=1))}
    t4 = {}
    for alg in series:
        def pooled(names):
            num = sum((t5[alg][nm]["n"] - 1) * t5[alg][nm]["sd_mhz"] ** 2 for nm in names)
            den = sum(t5[alg][nm]["n"] - 1 for nm in names)
            return math.sqrt(num / den), sum(t5[alg][nm]["n"] for nm in names)
        st, n_st = pooled(STATIONARY)
        dy, n_dy = pooled(DYNAMIC)
        t4[alg] = {"stationary_sd_mhz": st, "n_stationary": n_st,
                   "ratio_to_crb": st / CRB_MHZ,
                   "dynamic_sd_mhz": dy, "n_dynamic": n_dy}
    return t5, t4


# ---------------------------------------------------------------------------
# Table 7 - three-phase captures
# ---------------------------------------------------------------------------

def table7():
    out = {}
    for key, fn in (("balanced", BAL_CSV), ("thd_injected", THD_CSV)):
        sc = load_scada(find(fn))
        d = {}
        for name, a, b in REGIONS + [("Overall", 0, 359)]:
            e = np.array([abs(f - t) * 1000 for i, f, t in sc if a <= i <= b])
            d[name] = {"n": len(e), "max_abs_mhz": float(e.max()), "mean_abs_mhz": float(e.mean())}
        out[key] = d
    return out


# ---------------------------------------------------------------------------
# Table 8 - V6 parallel capture
# ---------------------------------------------------------------------------

def table8(t5):
    all_rows = []
    missing = set()
    with open(find(V6_CSV), newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["replay_index"] == "":
                continue
            i = int(r["replay_index"])
            if r["f_target_hz"] == "":
                missing.add(i)
                continue
            all_rows.append((i, (float(r["f_Vb"]) - float(r["f_target_hz"])) * 1000))
    d = {}
    for name, a, b in REGIONS + [("Overall", 0, 359)]:
        e = np.array([x for i, x in all_rows if a <= i <= b])
        v2sd = t5["V2"]["Aggregate (aligned)" if name == "Overall" else name]["sd_mhz"]
        d[name] = {"n": len(e), "max_abs_mhz": float(np.abs(e).max()),
                   "mean_abs_mhz": float(np.abs(e).mean()),
                   "signed_mean_mhz": float(e.mean()), "sd_mhz": float(e.std(ddof=1)),
                   "v2_sd_mhz": v2sd, "sd_ratio": v2sd / float(e.std(ddof=1))}
    e = np.array([x for _, x in all_rows])
    cov = {"tier2_pct": float(100 * np.mean(np.abs(e) <= 100)),
           "tier1_pct": float(100 * np.mean(np.abs(e) <= 10))}
    idx = sorted(missing)
    # Sensitivity: score the unreferenced records against the published trajectory instead.
    traj = [float(r["frequency_hz"]) for r in csv.DictReader(open(find(TRAJ), newline="", encoding="utf-8"))]
    extra = []
    with open(find(V6_CSV), newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["replay_index"] != "" and r["f_target_hz"] == "":
                i = int(r["replay_index"])
                extra.append((float(r["f_Vb"]) - traj[i]) * 1000)
    e_inc = np.concatenate([e, np.array(extra)]) if extra else e
    sens = {"n": int(len(e_inc)), "mean_abs_mhz": float(np.abs(e_inc).mean()),
            "tier2_pct": float(100 * np.mean(np.abs(e_inc) <= 100))}
    return {"regions": d, "coverage": cov,
            "excluded_indices_without_reference": [idx[0], idx[-1]] if idx else [],
            "sensitivity_including_unreferenced": sens}


def table8_excluded_frames():
    n = 0
    with open(find(V6_CSV), newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["replay_index"] != "" and r["f_target_hz"] == "":
                n += 1
    return n


# ---------------------------------------------------------------------------
# Section 3.8 - repeatability on the Section 3.1 partition
# ---------------------------------------------------------------------------

def repeatability(t8):
    traj = R.load_trajectory(find(CANONICAL))
    tgt = np.repeat(traj, 1000 // R.CADENCE_MS)
    weights = {name: t8["regions"][name]["n"] for name, _, _ in REGIONS}
    wsum = sum(weights.values())
    runs = []
    for fn in REPEAT:
        grid, nrep = R.load_capture(find(fn))
        lag, cc = R.align(grid, tgt)
        err = np.abs(grid[lag:lag + len(tgt)] - tgt) * 1000.0
        reg = {}
        for name, a, b in REGIONS:
            seg = err[a * 50:(b + 1) * 50]
            reg[name] = float(np.nanmean(seg))
        overall = float(np.nanmean(err))
        weighted = sum(reg[n] * weights[n] for n in weights) / wsum
        runs.append({"file": fn, "overall_mean_mhz": overall,
                     "weighted_to_table8_mhz": weighted, "regions_mean_abs_mhz": reg})
    ov = [r["overall_mean_mhz"] for r in runs]
    wt = [r["weighted_to_table8_mhz"] for r in runs]
    deltas = {name: [r["regions_mean_abs_mhz"][name] - t8["regions"][name]["mean_abs_mhz"] for r in runs]
              for name, _, _ in REGIONS}
    return {"runs": runs,
            "overall_range_mhz": [min(ov), max(ov)], "weighted_range_mhz": [min(wt), max(wt)],
            "spread_overall_mhz": max(ov) - min(ov), "spread_weighted_mhz": max(wt) - min(wt),
            "region_minus_table8_mhz": deltas}


# ---------------------------------------------------------------------------
# Table 3 - TVE, all five regions
# ---------------------------------------------------------------------------

def table3():
    d = json.load(open(find(TVE_JSON), encoding="utf-8"))
    out = {}
    for name in [r[0] for r in REGIONS]:
        g = d["regions"][name]
        out[name] = {"frames": g["frame_count"], "tve_max_pct": g["tve_max_pct"],
                     "tve_mean_pct": g["tve_mean_pct"]}
    o = d["overall"]
    out["Aggregate (full replay)"] = {"frames": d["input"]["frames_processed"],
                                      "tve_max_pct": o["tve_max_pct"], "tve_mean_pct": o["tve_mean_pct"]}
    return out


# ---------------------------------------------------------------------------
# values printed in the submitted manuscript and their tolerances
# ---------------------------------------------------------------------------

EXPECTED = [
    ("Sec 2.3 realised floor (mHz)", lambda o: o["canonical"]["noise_floor_pre_event_mhz"], 5.46, 0.005),
    ("Sec 2.8 settled-tail floor (mHz)", lambda o: o["canonical"]["noise_floor_settled_mhz"], 21.0, 0.5),
    ("Sec 3.1 dynamic sd min (mHz)", lambda o: min(o["canonical"]["regions"][n]["sd_mhz"] for n in DYNAMIC), 98.5, 0.5),
    ("Sec 3.1 dynamic sd max (mHz)", lambda o: max(o["canonical"]["regions"][n]["sd_mhz"] for n in DYNAMIC), 259.3, 0.5),
    ("Table 2 pre p(10)", lambda o: o["table2"][0]["p_tost_delta10"], 0.9998, 0.0005),
    ("Table 2 settled p(10)", lambda o: o["table2"][1]["p_tost_delta10"], 0.43, 0.005),
    ("Table 2 full p(10)", lambda o: o["table2"][2]["p_tost_delta10"], 0.22, 0.005),
    ("Table 2 pre log10 p(100)", lambda o: math.log10(o["table2"][0]["p_tost_delta100"]), math.log10(2.3e-30), 0.05),
    ("Table 2 settled log10 p(100)", lambda o: math.log10(o["table2"][1]["p_tost_delta100"]), math.log10(5.3e-27), 0.05),
    ("Table 2 full log10 p(100)", lambda o: math.log10(o["table2"][2]["p_tost_delta100"]), math.log10(1.6e-18), 0.05),
    ("Table 2 pre 90% CI lo", lambda o: o["table2"][0]["ci90_t_mhz"][0], -26.1, 0.06),
    ("Table 2 pre 90% CI hi", lambda o: o["table2"][0]["ci90_t_mhz"][1], -16.2, 0.06),
    ("Table 2 settled 90% CI lo", lambda o: o["table2"][1]["ci90_t_mhz"][0], -17.4, 0.06),
    ("Table 2 settled 90% CI hi", lambda o: o["table2"][1]["ci90_t_mhz"][1], -0.9, 0.06),
    ("Table 2 full 90% CI lo", lambda o: o["table2"][2]["ci90_t_mhz"][0], -16.0, 0.06),
    ("Table 2 full 90% CI hi", lambda o: o["table2"][2]["ci90_t_mhz"][1], 19.4, 0.06),
    ("Sec 3.2 n_eff pre", lambda o: o["table2"][0]["n_eff"], 15.1, 0.1),
    ("Sec 3.2 n_eff settled", lambda o: o["table2"][1]["n_eff"], 11.7, 0.1),
    ("Sec 3.2 n_eff full", lambda o: o["table2"][2]["n_eff"], 93.7, 0.1),
    ("Sec 3.2 crossover(n_eff) pre", lambda o: o["table2"][0]["crossover_neff_mhz"], 30.5, 0.5),
    ("Sec 3.2 crossover(n_eff) settled", lambda o: o["table2"][1]["crossover_neff_mhz"], 29.9, 0.5),
    ("Sec 3.2 crossover(n_eff) full", lambda o: o["table2"][2]["crossover_neff_mhz"], 37.5, 0.5),
    ("Table 3 late recovery TVE mean (%)", lambda o: o["table3"]["Late recovery"]["tve_mean_pct"], 3.54, 0.006),
    ("Table 3 settled TVE mean (%)", lambda o: o["table3"]["Post-event settled"]["tve_mean_pct"], 1.32, 0.006),
    ("Table 4 V2 stationary sd (mHz)", lambda o: o["table4"]["V2"]["stationary_sd_mhz"], 32.8, 0.05),
    ("Table 4 V2 dynamic sd (mHz)", lambda o: o["table4"]["V2"]["dynamic_sd_mhz"], 202.2, 0.05),
    ("Table 4 V3 stationary sd (mHz)", lambda o: o["table4"]["V3"]["stationary_sd_mhz"], 43.2, 0.05),
    ("Table 4 V4 stationary sd (mHz)", lambda o: o["table4"]["V4"]["stationary_sd_mhz"], 30.8, 0.05),
    ("Table 4 V4 dynamic sd (mHz)", lambda o: o["table4"]["V4"]["dynamic_sd_mhz"], 178.4, 0.05),
    ("Table 8 pre-event sd ratio", lambda o: o["table8"]["regions"]["Pre-event nominal"]["sd_ratio"], 0.99, 0.005),
    ("Table 5 V1 aggregate sd (mHz)", lambda o: o["table5"]["V1"]["Aggregate (aligned)"]["sd_mhz"], 131.3, 0.05),
    ("Table 7 bal. late recovery mean (mHz)", lambda o: o["table7"]["balanced"]["Late recovery"]["mean_abs_mhz"], 200.0, 0.05),
    ("Table 7 THD overall mean (mHz)", lambda o: o["table7"]["thd_injected"]["Overall"]["mean_abs_mhz"], 169.1, 0.05),
    ("Table 8 overall mean |e| (mHz)", lambda o: o["table8"]["regions"]["Overall"]["mean_abs_mhz"], 78.3, 0.05),
    ("Table 8 Tier 2 coverage (%)", lambda o: o["table8"]["coverage"]["tier2_pct"], 68.9, 0.05),
    ("Table 8 excluded frames", lambda o: o["table8_excluded_frames"], 400, 0),
    ("Sec 3.8 weighted mean min (mHz)", lambda o: o["repeatability"]["weighted_range_mhz"][0], 71.3, 0.05),
    ("Sec 3.8 weighted mean max (mHz)", lambda o: o["repeatability"]["weighted_range_mhz"][1], 74.7, 0.05),
    ("Sec 3.8 overall spread (mHz)", lambda o: o["repeatability"]["spread_overall_mhz"], 3.0, 0.05),
    ("Sec 3.8 weighted spread (mHz)", lambda o: o["repeatability"]["spread_weighted_mhz"], 3.5, 0.05),
    ("Sec 2.3 reference sd, pre-event (mHz)", lambda o: o["canonical"]["reference_sd_pre_event_mhz"], 20.1, 0.05),
    ("Sec 3.2 canonical Tier 2 frame coverage (%)", lambda o: o["canonical"]["coverage"]["tier2_pct"], 44.3, 0.05),
    ("Sec 3.2 canonical Tier 1 frame coverage (%)", lambda o: o["canonical"]["coverage"]["tier1_pct"], 6.6, 0.05),
    ("Sec 3.8 mean |e| incl. unreferenced (mHz)", lambda o: o["table8"]["sensitivity_including_unreferenced"]["mean_abs_mhz"], 79.1, 0.05),
    ("Sec 3.8 Tier 2 incl. unreferenced (%)", lambda o: o["table8"]["sensitivity_including_unreferenced"]["tier2_pct"], 68.3, 0.05),
    ("Sec 3.1 canonical mean |e| (mHz)", lambda o: o["canonical"]["mean_abs_error_mhz"], 159.4, 0.05),
]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="verify against the manuscript values")
    ap.add_argument("--json", default=os.path.join(ROOT, "paper_tables_output.json"))
    a = ap.parse_args(argv)

    sc = load_scada(find(CANONICAL))
    rows = [(i, (f - t) * 1000) for i, f, t in sc]
    out = {"schema_version": "1.0", "tool": "paper_tables.py",
           "regions": [list(r) for r in REGIONS]}
    out["canonical"] = canonical_stats(sc)
    out["canonical"]["mean_abs_error_mhz"] = float(np.mean([abs(f - t) * 1000 for i, f, t in sc]))
    out["table2"] = table2(rows)
    out["table3"] = table3()
    t5, t4 = four_estimators()
    out["table4"] = t4
    out["table5"] = t5
    out["table7"] = table7()
    out["table8"] = table8(t5)
    out["table8_excluded_frames"] = table8_excluded_frames()
    out["repeatability"] = repeatability(out["table8"])

    with open(a.json, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(out, fh, indent=2)

    print("\nTable 2 (exact Student-t TOST; alpha = 0.05)")
    for r in out["table2"]:
        print(f"  {r['window']:34s} n={r['n']:3d}  mean={r['mean_mhz']:+7.2f}  sd={r['sd_mhz']:6.1f}  "
              f"90%CI=[{r['ci90_t_mhz'][0]:+.1f},{r['ci90_t_mhz'][1]:+.1f}]  "
              f"boot95=[{r['ci95_bootstrap_mhz'][0]:+.1f},{r['ci95_bootstrap_mhz'][1]:+.1f}]  "
              f"p10={r['p_tost_delta10']:.4g}  p100={r['p_tost_delta100']:.2g}  "
              f"rho1={r['lag1_autocorrelation']:.2f}  n_eff={r['n_eff']:.1f}  "
              f"x_neff={r['crossover_neff_mhz']:.1f}  p100_neff={r['p_tost_delta100_neff']:.2g}")
    print("\nTable 3 (TVE-equivalent)")
    for k, v in out["table3"].items():
        print(f"  {k:26s} {v['frames']:4d}  max {v['tve_max_pct']:.2f} %  mean {v['tve_mean_pct']:.2f} %")
    print("\nTable 4 (regional columns)")
    for k, v in t4.items():
        print(f"  {k}: stationary {v['stationary_sd_mhz']:.1f} (n={v['n_stationary']})  "
              f"ratio {v['ratio_to_crb']:.1f}x  dynamic {v['dynamic_sd_mhz']:.1f} (n={v['n_dynamic']})")
    print("\nTable 5 (mean / sd, mHz, index-aligned)")
    for name in [r[0] for r in REGIONS] + ["Aggregate (aligned)"]:
        cells = "  ".join(f"{t5[a][name]['mean_mhz']:+.1f} / {t5[a][name]['sd_mhz']:.1f}" for a in ("V1", "V2", "V3", "V4"))
        print(f"  {name:24s} n={t5['V1'][name]['n']:3d}  {cells}")
    print("\nTable 7 (|error|, mHz)")
    for name in [r[0] for r in REGIONS] + ["Overall"]:
        b = out["table7"]["balanced"][name]
        t = out["table7"]["thd_injected"][name]
        print(f"  {name:24s} n={b['n']:3d}/{t['n']:3d}  bal max {b['max_abs_mhz']:.1f} mean {b['mean_abs_mhz']:.1f}  "
              f"THD max {t['max_abs_mhz']:.1f} mean {t['mean_abs_mhz']:.1f}")
    print("\nTable 8 (V6)")
    for name in [r[0] for r in REGIONS] + ["Overall"]:
        v = out["table8"]["regions"][name]
        print(f"  {name:24s} n={v['n']:6d} max {v['max_abs_mhz']:.1f} mean|e| {v['mean_abs_mhz']:.1f} "
              f"signed {v['signed_mean_mhz']:+.1f} sd {v['sd_mhz']:.1f} V2 sd {v['v2_sd_mhz']:.1f} ratio {v['sd_ratio']:.2f}")
    print(f"  coverage: Tier 2 {out['table8']['coverage']['tier2_pct']:.1f} %, Tier 1 {out['table8']['coverage']['tier1_pct']:.1f} %; "
          f"excluded (no reference value): indices {out['table8']['excluded_indices_without_reference']}, "
          f"{out['table8_excluded_frames']} frames")
    rp = out["repeatability"]
    print("\nSection 3.8 repeatability (Section 3.1 partition)")
    for r in rp["runs"]:
        print(f"  {r['file']}: overall {r['overall_mean_mhz']:.1f}  weighted {r['weighted_to_table8_mhz']:.1f}  "
              + "  ".join(f"{k.split()[0]} {v:.1f}" for k, v in r["regions_mean_abs_mhz"].items()))
    print(f"  spread overall {rp['spread_overall_mhz']:.1f} mHz, weighted {rp['spread_weighted_mhz']:.1f} mHz")
    for k, v in rp["region_minus_table8_mhz"].items():
        print(f"    {k:24s} run minus Table 8: " + ", ".join(f"{x:+.1f}" for x in v))
    c = out["canonical"]
    print(f"\nSection 2.3 floor {c['noise_floor_pre_event_mhz']:.3f} mHz (n-1), "
          f"{c['noise_floor_pre_event_population_norm_mhz']:.3f} mHz (n); settled tail {c['noise_floor_settled_mhz']:.2f} mHz")

    if a.check:
        bad = 0
        print("\nChecks against the manuscript:")
        for label, fn, exp, tol in EXPECTED:
            got = fn(out)
            ok = abs(got - exp) <= tol
            bad += not ok
            print(f"  [{'OK' if ok else 'FAIL'}] {label}: got {got:.4g}, expected {exp:.4g} (tol {tol})")
        print("ALL TABLE CHECKS PASSED" if not bad else f"{bad} TABLE CHECK(S) FAILED")
        return 1 if bad else 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
