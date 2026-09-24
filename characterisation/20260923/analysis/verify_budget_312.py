#!/usr/bin/env python3
"""
verify_budget_312.py - the error-budget numbers of Section 3.12 and Figure 20,
recomputed from the released records.

  (1) Predicted dispersion from measured noise (bars 1, 2: paper 1.69, 3.59 mHz).
      No released script produces these.  Two documented propagations are given:
        analytic  sigma_f = f0^2 * sqrt2 * sigma_t / k,  sigma_t = sigma_v / (2 pi f0 A),
                  k = 4 same-polarity periods per chain in a 200-sample window at
                  2 kHz (the chain of characterise_loop.py, pooled over one chain);
        Monte Carlo  the exact port of measureFreq() (analyse_ac.estimate_freq:
                  median, integer-us truncation) on synthetic 200-sample windows
                  at 2 kHz, A = 1.0928 V, Gaussian noise sigma_v, uniformly random
                  phase per window, 20,000 windows, fixed seed.  Also reported with
                  mean pooling (estimate_from_codes crossings), which is the
                  statistic of bar 3.
  (2) Pooling statistic (bars 3, 4: paper 3.08, 5.06 mHz, penalty 1.64).
      estimate_from_codes.estimate() on 20260923/captures/F_stim_baseline.txt
      (mode 0, production stimulus, 23 Sep), mean and median pooling.  The
      13 Sep B1 figure (analyse_ac.py, exact port, median) is shown alongside.
  (3) Timing walk (paper 0.010/0.684/1.014/1.208 us; 2.42/2.65/2.32 mHz; 4.27;
      4.03 mHz with stimulus).  characterise_loop.analyse() on the raw
      A_b*, B_b*, C_b*, D_b* and G_stim_v2arith logs; the "random jitter"
      (period-2 alternation removed) is averaged over the three records of each
      mode and propagated with characterise_loop's chain.
  (4) Crossing counts (paper 2,648 -> 902).  The crossings_pos / crossings_neg
      counters the mode-3 firmware prints in the footer of D_b1..3 and
      G_stim_v2arith, and an independent recount from the logged codes with the
      logged bias_volts and the same float arithmetic.
  (5) Predicted canonical floor (bar 5: 6.46 mHz; ~6.1 mHz on 1816.5 Hz basis).
  (6) Re-measurement of 23 Sep (bar 7: 7.97 mHz over 305 windows; 7.81 mHz
      successive difference) - parse_v2_serial.py on runE_v2.txt; Gaussian
      range check by simulation.
  (7) Derived percentages (4.53 mHz, 18 %, 19 %, 46 %).

Usage:  python3 verify_budget_312.py [--mc-windows 20000]
Deterministic.  numpy only.
"""
from __future__ import annotations
import argparse, math, os, statistics
import numpy as np
import common_312 as c
import characterise_loop as cl
import parse_v2_serial as pv

ARCHIVED_FLOOR = 5.46       # 27 Jun 2026 replay capture; not in this bundle


def analytic_mhz(sigma_mv, A=c.A_BENCH, k=4):
    return c.mhz_from_sigma_v(sigma_mv * 1e-3, A) / k


def mc_predict(sigma_mv, nwin, seed, A=c.A_BENCH, fs=2000.0):
    rng = np.random.default_rng(seed)
    t = np.arange(c.FREQ_WIN) / fs
    t_us = np.round(t * 1e6).astype(np.int64)
    tl = t_us.tolist()
    med, mean = [], []
    for _ in range(nwin):
        v = A * np.sin(2 * np.pi * c.F0 * t + rng.uniform(0, 2 * np.pi)) \
            + rng.normal(0, sigma_mv * 1e-3, c.FREQ_WIN)
        fe = c.measure_freq_median(v, t_us)
        if fe is not None:
            med.append(fe)
        per = []
        vl = v.tolist()
        for rising in (True, False):
            xs = c._efc.crossings(tl, vl, rising)[:12]
            per += [xs[i] - xs[i - 1] for i in range(1, len(xs))]
        mean.append(1e6 / statistics.mean(per))
    return np.std(med, ddof=1) * 1e3, np.std(mean, ddof=1) * 1e3


def recount_crossings(path):
    """Mode-3 firmware crossing counter, recomputed from the logged codes."""
    rows = c.read_tcodes(path)
    bias = float(c.header(path)["bias_volts"])
    v = [code * (3.3 / 4095.0) - bias for _, code in rows]
    npos = nneg = 0
    for a, b in zip(v, v[1:]):
        if a < 0.0 <= b:
            npos += 1
        elif a >= 0.0 > b:
            nneg += 1
    return npos, nneg, len(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mc-windows", type=int, default=20000)
    a = ap.parse_args(argv)
    out = {}
    W = 76

    # ---- (1) predictions from measured noise ------------------------------
    print("=" * W + "\n  (1) PREDICTED DISPERSION FROM MEASURED NOISE (Fig 20 bars 1-2)\n" + "=" * W)
    pred = {}
    for s, paper in ((0.659, 1.69), (1.399, 3.59)):
        an = analytic_mhz(s)
        mc_med, mc_mean = mc_predict(s, a.mc_windows, seed=int(s * 1000))
        pred[s] = {"analytic_k4": an, "mc_exact_port_median": mc_med,
                   "mc_mean_pooling": mc_mean, "paper": paper}
        print(f"  sigma_v {s:.3f} mV: analytic {an:.3f} | MC exact port (median) {mc_med:.3f}"
              f" | MC mean pooling {mc_mean:.3f} mHz   [paper {paper}]")
    out["predictions"] = pred

    # ---- (2) pooling ------------------------------------------------------
    print("=" * W + "\n  (2) POOLING STATISTIC (Fig 20 bars 3-4)\n" + "=" * W)
    F = os.path.join(c.CAP_0923, "F_stim_baseline.txt")
    rows = c.read_tcodes(F)
    em = c.estimate_windows(rows, 200, "mean", 50.0)
    ed = c.estimate_windows(rows, 200, "median", 50.0)
    codes, meta = c.read_codes(c.B1)
    v = codes.astype(float) * c.LSB_V
    v = v - v.mean()
    t_us = np.round(np.arange(c.FREQ_WIN) / float(meta["fs_hz"]) * 1e6).astype(np.int64)
    fb = [c.measure_freq_median(v[i:i + 200], t_us) for i in range(0, len(v) - 199, 200)]
    b1_med = np.std([x for x in fb if x], ddof=1) * 1e3
    print(f"  F_stim_baseline (23 Sep, mode 0 + stimulus), {em['n_windows']} windows:")
    print(f"    mean pooling   {em['sd_mhz']:.3f} mHz   [paper 3.08]")
    print(f"    median pooling {ed['sd_mhz']:.3f} mHz   [paper 5.06]")
    print(f"    penalty        {ed['sd_mhz']/em['sd_mhz']:.3f}      [paper 1.64]")
    print(f"  B1_ac_running (13 Sep), exact port, median: {b1_med:.3f} mHz  [paper 3.05; ac.json 3.047]")
    out["pooling"] = {"record": "20260923/captures/F_stim_baseline.txt",
                      "mean_mhz": em["sd_mhz"], "median_mhz": ed["sd_mhz"],
                      "penalty": ed["sd_mhz"] / em["sd_mhz"], "b1_13sep_median_mhz": b1_med}

    # ---- (3) timing walk --------------------------------------------------
    print("=" * W + "\n  (3) TIMING WALK\n" + "=" * W)
    modes = {}
    for m, name in zip("ABCD", ("0 deadline, ts after", "1 deadline, ts before",
                                "2 free-run", "3 free-run + V2 arithmetic")):
        sds, rj, mhz = [], [], []
        for i in (1, 2, 3):
            res = cl.analyse(cl.read_log(os.path.join(c.CAP_0923, f"{m}_b{i}.txt")),
                             50.0, 1.0939)
            sds.append(res["interval_us"]["sd"])
            rj.append(res["alternation"]["random_jitter_us"])
            mhz.append(res["alternation"]["mhz_from_random_jitter"])
        modes[m] = {"name": name, "interval_sd_us": sds, "random_jitter_us": rj,
                    "mean_random_jitter_us": float(np.mean(rj)),
                    "mean_interval_sd_us": float(np.mean(sds)),
                    "mean_mhz": float(np.mean(mhz)),
                    "half_range_over_mean_pct": 50 * (max(rj) - min(rj)) / np.mean(rj),
                    "range_over_mean_pct": 100 * (max(rj) - min(rj)) / np.mean(rj)}
        print(f"  mode {name:28s} random jitter {np.mean(rj):.4f} us "
              f"(interval SD {np.mean(sds):.4f}), {np.mean(mhz):.3f} mHz, "
              f"range/mean {modes[m]['range_over_mean_pct']:.1f} %")
    q = lambda x, y: math.sqrt(modes[y]["mean_mhz"] ** 2 - modes[x]["mean_mhz"] ** 2)
    attrib = {"ts_placement": q("A", "B"), "free_run": q("B", "C"), "arithmetic": q("C", "D")}
    print(f"  quadrature: placement {attrib['ts_placement']:.3f}, free-run "
          f"{attrib['free_run']:.3f}, arithmetic {attrib['arithmetic']:.3f} mHz"
          f"   [paper 2.42, 2.65, 2.32]; mode 3 {modes['D']['mean_mhz']:.3f} [4.27]")
    g = cl.analyse(cl.read_log(os.path.join(c.CAP_0923, "G_stim_v2arith.txt")), 50.0, 1.0939)
    g_mhz = g["alternation"]["mhz_from_random_jitter"]
    print(f"  mode 3 with stimulus (G_stim_v2arith): {g['alternation']['random_jitter_us']:.4f} us"
          f" -> {g_mhz:.3f} mHz   [paper 4.03]; rate {g['sampling_rate_hz']['from_mean']:.2f} Hz")
    out["timing_walk"] = {"modes": modes, "attribution_mhz": attrib, "G_stim_mhz": g_mhz}

    # ---- (4) crossing counts ----------------------------------------------
    print("=" * W + "\n  (4) CROSSING COUNTS, mode 3 firmware counters (whole 16,384-sample record)\n" + "=" * W)
    cross = {}
    for f in ("D_b1.txt", "D_b2.txt", "D_b3.txt", "G_stim_v2arith.txt"):
        p = os.path.join(c.CAP_0923, f)
        h = c.header(p)
        rp, rn, n = recount_crossings(p)
        cross[f] = {"fw_pos": int(h["crossings_pos"]), "fw_neg": int(h["crossings_neg"]),
                    "recount_pos": rp, "recount_neg": rn}
        print(f"  {f:20s} firmware pos {h['crossings_pos']:>5} neg {h['crossings_neg']:>5} "
              f"total {int(h['crossings_pos'])+int(h['crossings_neg']):>5} | recount "
              f"pos {rp} neg {rn}")
    qt = [cross[f]["fw_pos"] + cross[f]["fw_neg"] for f in ("D_b1.txt", "D_b2.txt", "D_b3.txt")]
    print(f"  quiet input, total both polarities: {min(qt)}-{max(qt)}, mean {np.mean(qt):.0f}; "
          f"stimulus: {cross['G_stim_v2arith.txt']['fw_pos']+cross['G_stim_v2arith.txt']['fw_neg']}")
    print("  [paper: 2,648 (= D_b1 POSITIVE-going only) vs 902 (= G, BOTH polarities)]")
    out["crossings"] = cross

    # ---- (5) predicted floor ----------------------------------------------
    pf = math.hypot(ed["sd_mhz"], g_mhz)
    pf_can = math.hypot(ed["sd_mhz"] * 1816.5 / 2000.0, g_mhz)
    print("=" * W + "\n  (5) PREDICTED CANONICAL FLOOR\n" + "=" * W)
    print(f"  sqrt({ed['sd_mhz']:.3f}^2 + {g_mhz:.3f}^2) = {pf:.3f} mHz   [paper 6.46]")
    print(f"  noise term scaled by 1816.5/2000: {pf_can:.3f} mHz   [paper ~6.1]")

    # ---- (6) run E --------------------------------------------------------
    print("=" * W + "\n  (6) RE-MEASUREMENT 23 SEP (runE_v2.txt)\n" + "=" * W)
    vals = pv.extract(os.path.join(c.CAP_0923, "runE_v2.txt"))
    e = pv.analyse(vals)
    rng = np.random.default_rng(305)
    sim_rng = [np.ptp(rng.normal(0, 1, e["n_windows"])) for _ in range(20000)]
    obs = (e["max_hz"] - e["min_hz"]) * 1e3 / e["sd_mhz"]
    pct = 100 * np.mean(np.array(sim_rng) <= obs)
    print(f"  windows {e['n_windows']}, SD {e['sd_mhz']:.3f} mHz, successive-difference "
          f"{e['successive_difference_mhz']:.3f} mHz   [paper 305, 7.97, 7.81]")
    print(f"  range {(e['max_hz']-e['min_hz'])*1e3:.1f} mHz = {obs:.2f} SD; Gaussian n=305 "
          f"median range {np.median(sim_rng):.2f} SD, observed at percentile {pct:.0f}")
    out["runE"] = e | {"range_in_sd": obs, "gaussian_range_percentile": pct}

    # ---- (7) derived ------------------------------------------------------
    print("=" * W + "\n  (7) DERIVED\n" + "=" * W)
    d = {"apparent_residual": math.sqrt(ARCHIVED_FLOOR ** 2 - 3.05 ** 2),
         "pred_over_archived_pct": 100 * (pf / ARCHIVED_FLOOR - 1),
         "pred_below_runE_pct": 100 * (1 - pf / e["sd_mhz"]),
         "runE_over_archived_pct": 100 * (e["sd_mhz"] / ARCHIVED_FLOOR - 1)}
    print(f"  sqrt(5.46^2-3.05^2) = {d['apparent_residual']:.3f} [4.53]; prediction "
          f"{d['pred_over_archived_pct']:+.1f} % vs 5.46 [18]; {d['pred_below_runE_pct']:.1f} % "
          f"below 7.97 [19]; floors differ {d['runE_over_archived_pct']:.1f} % [46]")
    out["predicted_floor_mhz"] = pf
    out["predicted_floor_1816_mhz"] = pf_can
    out["derived"] = d
    c.save_json("budget_312.json", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
