#!/usr/bin/env python3
"""
make_figure_20_error_budget.py - Figure 20 (numbered 19 in earlier drafts), the frequency error budget at the 200-sample
estimator basis, regenerated from the released records only.

Every bar is computed here, none is typed in:

  1  converter noise alone      analytic chain (characterise_loop.py) of the
                                13 Sep quiet-input mean sigma_v, pooled over the
                                k = 4 periods of one crossing chain:
                                sigma_f = f0^2 * sqrt2 * (sigma_v / 2 pi f0 A) / 4
  2  with stimulus present      the same chain at the broadband residual of
                                B1_ac_running (verify_stimulus_residual_312.py)
  3  mean pooling               estimate_from_codes.py --pooling mean on
                                20260923/captures/F_stim_baseline.txt
  4  median pooling             the same record, --pooling median
  5  predicted canonical floor  hypot(bar 4, mode-3 timing term measured with the
                                stimulus, G_stim_v2arith.txt, characterise_loop.py)
  6  archived replay 27 Jun     repository paper_tables_output.json,
                                canonical.noise_floor_pre_event_mhz
                                (analysis/paper_tables.py on
                                captures/replay_20260627_102507V2.csv)
  7  re-measured 23 Sep         parse_v2_serial.py on runE_v2.txt

Usage:
  python3 make_figure_20_error_budget.py [--paper-tables PATH] [--out PATH]
PATH defaults to the repository-root paper_tables_output.json
(characterisation/../paper_tables_output.json).
"""
from __future__ import annotations
import argparse, json, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import common_312 as c
import characterise_loop as cl
import parse_v2_serial as pv
import verify_stimulus_residual_312 as vs

COL = {"pred": "#1f6fb4", "meas": "#e0651c", "canon": "#2a8466", "real": "#8c4fc6"}


def compute(paper_tables):
    # bar 1: 13 Sep quiet-input mean (detrended, analyse_noise.py method)
    sig_q = np.mean([c.sd_raw_and_detrended(c.read_codes(p)[0])[1] for p in c.SEPT_QUIET]) * c.LSB_V
    # bar 2: broadband residual of B1
    codes, meta = c.read_codes(c.B1)
    fs = float(meta["fs_hz"])
    v = codes.astype(float) * c.LSB_V
    f1 = c.refine_frequency(v[:4000] - v.mean(), np.arange(4000) / fs, 50.0)
    r, _, _ = vs.window_residuals(v, fs, f1)
    fr, P = vs.power_spectrum(r, fs)
    m = vs.harmonic_mask(fr, f1, 5)
    sig_bb = math.sqrt(P[~m].sum())
    b1 = c.mhz_from_sigma_v(sig_q) / 4
    b2 = c.mhz_from_sigma_v(sig_bb) / 4
    rows = c.read_tcodes(os.path.join(c.CAP_0923, "F_stim_baseline.txt"))
    b3 = c.estimate_windows(rows, 200, "mean", 50.0)["sd_mhz"]
    b4 = c.estimate_windows(rows, 200, "median", 50.0)["sd_mhz"]
    g = cl.analyse(cl.read_log(os.path.join(c.CAP_0923, "G_stim_v2arith.txt")), 50.0, 1.0939)
    tj = g["alternation"]["mhz_from_random_jitter"]
    b5 = math.hypot(b4, tj)
    b6 = json.load(open(paper_tables))["canonical"]["noise_floor_pre_event_mhz"]
    b7 = pv.analyse(pv.extract(os.path.join(c.CAP_0923, "runE_v2.txt")))["sd_mhz"]
    return dict(sig_q=sig_q * 1e3, sig_bb=sig_bb * 1e3, timing=tj,
                bars=[b1, b2, b3, b4, b5, b6, b7])


def draw(d, out):
    b = d["bars"]
    labels = [f"Converter noise alone\n(predicted from measured $\\sigma_v$ = {d['sig_q']:.3f} mV)",
              f"With stimulus present\n(predicted from measured $\\sigma_v$ = {d['sig_bb']:.3f} mV)",
              "Estimator, mean pooling\n(measured, 23 Sep 2026)",
              "Same samples, median pooling as the firmware\n(measured, 23 Sep 2026)",
              f"Plus measured loop timing, {d['timing']:.2f} mHz\n(predicted canonical floor)",
              "Archived replay capture\n(27 Jun 2026)",
              "Same build re-measured\n(23 Sep 2026)"]
    cols = [COL["pred"], COL["pred"], COL["meas"], COL["meas"], COL["canon"], COL["real"], COL["real"]]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    fig, ax = plt.subplots(figsize=(10.0, 4.6), dpi=250)
    y = np.arange(len(b))[::-1]
    ax.barh(y, b, height=0.62, color=cols, edgecolor="white", linewidth=1.0, zorder=3)
    for yi, v in zip(y, b):
        ax.text(v + 0.12, yi, f"{v:.2f}", va="center", ha="left", fontsize=10,
                fontweight="bold", color="#1a1a1a")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.6, color="#333333")
    ax.set_xlim(0, 9.6)
    ax.set_xticks(range(0, 10))
    ax.set_ylim(-1.7, len(b) - 0.5)
    ax.set_xlabel("Single-estimate frequency dispersion (mHz), 200-sample estimator basis")
    ax.grid(axis="x", color="#cccccc", lw=0.8, zorder=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    lo, hi, p = b[5], b[6], b[4]
    ax.annotate("", xy=(lo, -0.9), xytext=(hi, -0.9),
                arrowprops=dict(arrowstyle="<->", color="#555555", lw=1.0))
    ax.plot([p, p], [-1.15, -0.65], color=COL["canon"], lw=2)
    ax.text((lo + hi) / 2, -1.45, f"prediction {p:.2f} mHz lies between the two realised floors",
            ha="center", va="center", fontsize=8.6, style="italic", color="#444444")
    ax.legend(handles=[Patch(color=COL["pred"], label="Predicted from measured noise"),
                       Patch(color=COL["meas"], label="Measured, this work"),
                       Patch(color=COL["canon"], label="Predicted canonical floor"),
                       Patch(color=COL["real"], label="Realised by the unmodified build")],
              loc="upper center", bbox_to_anchor=(0.36, -0.14), ncol=4, frameon=False,
              fontsize=8.6)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    print(f"  wrote {out}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paper-tables", default=os.path.normpath(
        os.path.join(c.HERE, "..", "..", "..", "paper_tables_output.json")))
    ap.add_argument("--out", default=os.path.join(c.FIG_DIR,
                                                  "Figure_20_error_budget.png"))
    a = ap.parse_args(argv)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    d = compute(a.paper_tables)
    paper = [1.69, 3.59, 3.08, 5.06, 6.46, 5.46, 7.97]
    for i, (v, p) in enumerate(zip(d["bars"], paper), 1):
        print(f"  bar {i}: {v:.3f} mHz  (paper {p})")
    draw(d, a.out)
    c.save_json("figure_20_values.json", d)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
