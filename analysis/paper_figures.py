#!/usr/bin/env python3
"""
paper_figures.py
================

Regenerates Figures 11-15, 17 and 18 of the manuscript (as submitted) from
the released captures, with the terminology and region partition the text
uses:

  * the two accuracy thresholds are labelled "Tier 1" (+/-10 mHz, 0.02 %) and
    "Tier 2" (+/-100 mHz, 0.2 %), the paper's own Product Design
    Specification tiers (Section 2.8). Earlier renders labelled them as
    IEC 61400-21 classes, which they are not;
  * region shading on the 9 August 2019 replays follows the Section 3.1
    partition by replay index (0-45, 46-119, 120-200, 201-300, 301-359);
  * the published record is labelled "NESO reference", not "truth": it is a
    measurement, and the paper claims no traceability through it.

Figure 16 (V7.2 DMA replay) is regenerated only when its raw capture, which
is available from the corresponding author on request, is supplied with
--dma PATH.

Usage (from the repository root or from analysis/):
    python analysis/paper_figures.py                 # writes figures/*.png
    python analysis/paper_figures.py --outdir OUT
    python analysis/paper_figures.py --dma replay_20260718_112046_dma_v7_2.txt

Output: 600 dpi at the 140 mm column width used in the manuscript.
Requires numpy and matplotlib.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
from matplotlib.transforms import blended_transform_factory  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE) if os.path.isdir(os.path.join(os.path.dirname(HERE), "analysis")) else HERE
sys.path.insert(0, HERE)

import tost_sweep as TS                  # noqa: E402
import bland_altman as BA                # noqa: E402
import inertia_estimator as IE           # noqa: E402
import analyse_v6_repeatability as R     # noqa: E402
import paper_tables as PT                # noqa: E402

DPI = 360            # figure widths of 10 in -> 3600 px; placed at 140 mm -> ~650 dpi

REGIONS_0809 = [("Pre-event nominal", 0, 45, "#dfeaf8"),
                ("RoCoF descent", 46, 119, "#fbe3d6"),
                ("Nadir + early recovery", 120, 200, "#fbf3d5"),
                ("Late recovery", 201, 300, "#e0f1e2"),
                ("Post-event settled", 301, 359, "#ececec")]
REGIONS_DEC23 = [("Pre-event plateau", 0, 99, "#dfeaf8"),
                 ("Collapse", 100, 119, "#fbe3d6"),
                 ("Nadir and partial recovery", 120, 179, "#fbf3d5"),
                 ("Slow recovery", 180, 359, "#e0f1e2")]


def find(name):
    return PT.find(name)


def load_traj(name="neso_trajectory.csv"):
    out = {}
    with open(find(name), newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            out[int(r["index"])] = float(r["frequency_hz"])
    return np.array([out[i] for i in sorted(out)])


def load_traj_dec():
    with open(find("neso_trajectory_December2023Peak.csv"), newline="", encoding="utf-8") as fh:
        return np.array([float(r["f_neso_Hz"]) for r in csv.DictReader(fh)])


# ---------------------------------------------------------------------------
# replay overlay (Figures 11, 15, 16, 17)
# ---------------------------------------------------------------------------

# The replay figures carry no in-image title: the manuscript caption says what each
# one shows, and a title repeating it reads as a second caption. The working title
# (capture file and summary statistics) is printed to the console instead, and
# --titles puts it back on the figure.
SHOW_TITLES = False


def overlay(x_meas, f_meas, x_err, err_mhz, traj, regions, title, meas_label, out_path,
            ylim=None, lower_ylim=None):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True,
                                   gridspec_kw={"height_ratios": [3, 1.05]})
    ax1.plot(x_meas, f_meas, color="#d62728", linewidth=0.8, linestyle=":",
             label=meas_label, zorder=2)
    ax1.plot(np.arange(len(traj)), traj, color="#1f77b4", linewidth=2.0,
             label="NESO reference", zorder=3)
    ax1.axhline(49.5, color="#888888", linestyle=":", linewidth=0.8, zorder=1)
    ax1.text(2, 49.5, "49.5 Hz statutory limit", color="#444444", fontsize=8,
             va="bottom", ha="left", bbox=dict(boxstyle="round,pad=0.1", facecolor="white",
                                                edgecolor="none", alpha=0.8), zorder=4)
    ax1.set_ylabel("Frequency (Hz)")
    if SHOW_TITLES:
        ax1.set_title(title, fontsize=11)
    else:
        for line in title.splitlines():
            print(f"    {line}")
    ax1.grid(True, alpha=0.3, linestyle=":")
    ax1.legend(loc="center right", fontsize=9, framealpha=0.95)
    if ylim:
        ax1.set_ylim(*ylim)

    for name, a, b, col in regions:
        ax2.axvspan(a - 0.5, b + 0.5, color=col, alpha=0.9, zorder=0, linewidth=0)
    pct = np.asarray(err_mhz) / 50.0 / 10.0          # mHz -> % of 50 Hz
    ax2.plot(x_err, pct, color="#9467bd", linewidth=0.7, zorder=3)
    top = lower_ylim if lower_ylim else max(0.25, float(np.nanmax(pct)) * 1.12)
    ax2.set_ylim(0, top)
    ax2.axhline(0.02, color="#555555", linewidth=0.9, zorder=2)
    ax2.axhline(0.20, color="#555555", linewidth=0.9, zorder=2)
    tr = blended_transform_factory(ax2.transAxes, ax2.transData)
    box = dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.85)
    ax2.text(0.995, 0.20, "Tier 2 (±0.2 %)", transform=tr, fontsize=8, ha="right", va="bottom",
             color="#222222", zorder=5, bbox=box)
    ax2.text(0.995, 0.02, "Tier 1 (±0.02 %)", transform=tr, fontsize=8, ha="right", va="bottom",
             color="#222222", zorder=5, bbox=box)
    for name, a, b, _ in regions:
        ax2.text((a + b) / 2.0, top * 0.95, name, fontsize=7.5, ha="center", va="top",
                 color="#333333", zorder=4)
    ax2.set_xlim(0, 360)
    ax2.set_xlabel("Replay index (1 s per step, aligned)")
    ax2.set_ylabel("|error| (% of 50 Hz)")
    ax2.grid(True, alpha=0.3, linestyle=":", zorder=1)
    fig.tight_layout()
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    print(f"  wrote {out_path}")


def fig11(outdir):
    sc = PT.load_scada(find(PT.CANONICAL))
    sc.sort(key=lambda r: r[0])
    traj = load_traj()
    idx = np.array([i for i, _, _ in sc], dtype=float)
    f = np.array([fv for _, fv, _ in sc])
    e = np.abs(f - traj[idx.astype(int)]) * 1000
    title = (f"NESO replay — {PT.CANONICAL} (V2 canonical firmware, 1 s reporting)\n"
             f"max |error| = {e.max():.1f} mHz, mean |error| = {e.mean():.1f} mHz, n = {len(e)} frames")
    overlay(idx, f, idx, e, traj, REGIONS_0809, title, "Measured $f_{\\mathrm{Vb}}$ (V2 canonical)",
            os.path.join(outdir, "Figure_11_v2_canonical_replay.png"))


def _v6_like(path):
    """SCADA-style 20 ms capture with replay_index / f_target_hz columns."""
    xs, fs, xe, es = [], [], [], []
    count = {}
    traj = load_traj()
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["replay_index"] == "":
                continue
            i = int(r["replay_index"])
            k = count.get(i, 0)
            count[i] = k + 1
            xs.append((i, k, float(r["f_Vb"]), r["f_target_hz"]))
    per = {}
    for i, k, fv, t in xs:
        per[i] = max(per.get(i, 0), k + 1)
    X, F, XE, E = [], [], [], []
    for i, k, fv, t in xs:
        x = i + k / max(per[i], 50)
        X.append(x)
        F.append(fv)
        if t != "":
            XE.append(x)
            E.append(abs(fv - float(t)) * 1000)
        else:
            XE.append(x)
            E.append(np.nan)
    return np.array(X), np.array(F), np.array(XE), np.array(E), traj


def fig15(outdir):
    X, F, XE, E, traj = _v6_like(find(PT.V6_CSV))
    n = int(np.sum(~np.isnan(E)))
    title = (f"NESO replay — replay_20260716_152206_parallel (V6 parallel firmware, 20 ms cadence)\n"
             f"max |error| = {np.nanmax(E):.1f} mHz, mean |error| = {np.nanmean(E):.1f} mHz, "
             f"n = {n:,} scored records of 19,997 captured")
    overlay(X, F, XE, E, traj, REGIONS_0809, title, "Measured $f_{\\mathrm{Vb}}$ (V6, 20 ms)",
            os.path.join(outdir, "Figure_15_v6_parallel_replay.png"))


def _raw_report(path, traj):
    grid, nrep = R.load_capture(path)
    tgt = np.repeat(traj, 1000 // R.CADENCE_MS)
    lag, cc = R.align(grid, tgt)
    seg = grid[lag:lag + len(tgt)]
    x = np.arange(len(tgt)) / (1000 // R.CADENCE_MS)
    err = np.abs(seg - tgt) * 1000.0
    return x, seg, err, nrep, lag, cc


def fig17(outdir):
    traj = load_traj_dec()
    x, seg, err, nrep, lag, cc = _raw_report(find("replay_20260725_195552_dec2023peak_v6.csv"), traj)
    valid = ~np.isnan(err)
    t2 = 100 * np.mean(err[valid] <= 100)
    t1 = 100 * np.mean(err[valid] <= 10)
    title = (f"December 2023 replay — replay_20260725_195552_dec2023peak_v6.csv (V6 parallel firmware, 20 ms)\n"
             f"max |error| = {np.nanmax(err):.1f} mHz, mean |error| = {np.nanmean(err):.1f} mHz, "
             f"n = {int(valid.sum()):,} in-replay records; Tier 2 {t2:.1f} %, Tier 1 {t1:.1f} %")
    overlay(x, seg, x, err, traj, REGIONS_DEC23, title, "Measured $f_{\\mathrm{Vb}}$ (V6, 20 ms)",
            os.path.join(outdir, "Figure_17_december2023_replay.png"))
    print(f"    Dec 2023: lag {lag * R.CADENCE_MS / 1000:.2f} s, corr {cc:.4f}, REPORT {nrep}, "
          f"mean {np.nanmean(err):.1f}, max {np.nanmax(err):.1f}, Tier2 {t2:.1f} %, Tier1 {t1:.1f} %")


def fig16(outdir, dma_path):
    traj = load_traj()
    x, seg, err, nrep, lag, cc = _raw_report(dma_path, traj)
    valid = ~np.isnan(err)
    title = (f"NESO replay — {os.path.basename(dma_path)} (V7.2 DMA firmware, 20 ms cadence)\n"
             f"max |error| = {np.nanmax(err):.1f} mHz, mean |error| = {np.nanmean(err):.1f} mHz, "
             f"n = {int(valid.sum()):,} in-replay records")
    overlay(x, seg, x, err, traj, REGIONS_0809, title, "Measured $f_{\\mathrm{Vb}}$ (V7.2 DMA, 20 ms)",
            os.path.join(outdir, "Figure_16_v72_dma_replay.png"))


# ---------------------------------------------------------------------------
# Figure 12 - TOST sweep
# ---------------------------------------------------------------------------

def fig12(outdir):
    rows = TS.read_canonical_csv(find(PT.CANONICAL))
    grid = [TS.DELTA_MIN_mHz + k * TS.DELTA_STEP_mHz
            for k in range(int((TS.DELTA_MAX_mHz - TS.DELTA_MIN_mHz) / TS.DELTA_STEP_mHz) + 1)]
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = ["#1f77b4", "#2ca02c", "#d62728"]
    markers = ["o", "s", "^"]
    labels = ["Pre-event nominal (n = 48)", "Post-event settled (n = 63)", "Full 360 s replay (n = 379)"]
    for (wname, a, b), c, m, lab in zip(TS.WINDOWS, colors, markers, labels):
        e = TS.extract_window_errors(rows, a, b)
        res, cross = TS.sweep_window(e, grid)
        ax.plot([r["delta_mHz"] for r in res], [r["p_tost"] for r in res], color=c, marker=m,
                markersize=4, linewidth=1.4, label=lab)
        if cross is not None:
            ax.axvline(cross, color=c, linestyle=":", linewidth=1.1, alpha=0.7)
    ax.axhline(TS.ALPHA, color="black", linestyle="--", linewidth=1.0, alpha=0.75, label="α = 0.05")
    for d, lab in ((10.0, "Tier 1 (10 mHz)"), (100.0, "Tier 2 (100 mHz)")):
        ax.axvline(d, color="#7f7f7f", linewidth=1.0, alpha=0.7)
        ax.text(d, 1.035, lab, ha="center", va="bottom", fontsize=9, color="#444444",
                transform=blended_transform_factory(ax.transData, ax.transAxes))
    ax.set_xlim(TS.DELTA_MIN_mHz, TS.DELTA_MAX_mHz)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Equivalence margin δ (mHz)")
    ax.set_ylabel("TOST p-value, max(p₁, p₂)")
    ax.legend(loc="center right", fontsize=9, framealpha=0.95)
    ax.grid(True, linestyle=":", alpha=0.4)
    fig.tight_layout()
    p = os.path.join(outdir, "Figure_12_tost_sensitivity_sweep.png")
    fig.savefig(p, dpi=DPI)
    plt.close(fig)
    print(f"  wrote {p}")


# ---------------------------------------------------------------------------
# Figure 13 - Bland-Altman
# ---------------------------------------------------------------------------

def fig13(outdir):
    frames = BA.load_canonical(find(PT.CANONICAL))
    frames = [f for f in frames if f["replay_idx"] >= 0]
    diffs = [(f["f_measured"] - f["f_target"]) * 1000 for f in frames]
    means = [(f["f_measured"] + f["f_target"]) / 2.0 for f in frames]
    st = BA.stats_for_errors(diffs, [f["f_target"] for f in frames], [f["f_measured"] for f in frames])
    names = {"pre-event nominal": "Pre-event nominal", "descent": "RoCoF descent",
             "nadir": "Nadir + early recovery", "late-recovery": "Late recovery",
             "settled": "Post-event settled"}
    fig, ax = plt.subplots(figsize=(10, 6.6))
    ax.axhspan(-100, 100, alpha=0.10, color="green", zorder=0)
    ax.text(min(means), 96, "Tier 2 band (±100 mHz)", fontsize=9, ha="left", va="top",
            color="#2d6b2d", style="italic")
    for rn, lo, hi in BA.REGIONS:
        xs = [means[i] for i, f in enumerate(frames) if lo <= f["replay_idx"] <= hi]
        ys = [diffs[i] for i, f in enumerate(frames) if lo <= f["replay_idx"] <= hi]
        ax.scatter(xs, ys, c=BA.REGION_COLOURS[rn], s=22, alpha=0.8, edgecolors="none", label=names[rn])
    ax.axhline(st["ba_bias_mHz"], color="#222222", linewidth=1.3, label=f"bias = {st['ba_bias_mHz']:+.2f} mHz")
    ax.axhline(st["ba_loa_upper_mHz"], color="#222222", linestyle="--", linewidth=1.0,
               label=f"bias + 1.96 SD = {st['ba_loa_upper_mHz']:+.1f} mHz")
    ax.axhline(st["ba_loa_lower_mHz"], color="#222222", linestyle="-.", linewidth=1.0,
               label=f"bias − 1.96 SD = {st['ba_loa_lower_mHz']:+.1f} mHz")
    ax.set_xlabel("(measured + reference) / 2 (Hz)")
    ax.set_ylabel("measured − reference (mHz)")
    ax.grid(True, alpha=0.3, linestyle=":")
    ax.legend(loc="lower right", fontsize=9, framealpha=0.95)
    fig.tight_layout()
    p = os.path.join(outdir, "Figure_13_bland_altman.png")
    fig.savefig(p, dpi=DPI)
    plt.close(fig)
    print(f"  wrote {p}  (bias {st['ba_bias_mHz']:+.2f}, LoA {st['ba_loa_lower_mHz']:+.1f} / {st['ba_loa_upper_mHz']:+.1f})")


# ---------------------------------------------------------------------------
# Figure 14 - per-region dispersion of the four estimators
# ---------------------------------------------------------------------------

def fig14(outdir):
    t5, _ = PT.four_estimators()
    regs = [r[0] for r in PT.REGIONS]
    ticks = ["Pre-event\nnominal", "RoCoF\ndescent", "Nadir + early\nrecovery", "Late\nrecovery", "Post-event\nsettled"]
    algs = [("V1", "V1 baseline", "white", "//"), ("V2", "V2 canonical", "#4a4a4a", ""),
            ("V3", "V3 hybrid", "white", ".."), ("V4", "V4 Kalman", "#b0b0b0", "")]
    fig, ax = plt.subplots(figsize=(10, 5.2))
    w = 0.2
    for k in (0, 4):
        ax.axvspan(k - 0.5, k + 0.5, color="#f0f0f0", zorder=0)
        ax.text(k, 297, "quiet", ha="center", va="top", fontsize=9, style="italic", color="#666666")
    ax.text(2, 297, "dynamic", ha="center", va="top", fontsize=9, style="italic", color="#666666")
    for j, (a, lab, fc, h) in enumerate(algs):
        vals = [t5[a][r]["sd_mhz"] for r in regs]
        ax.bar(np.arange(5) + (j - 1.5) * w, vals, w, color=fc, edgecolor="#1a1a1a", hatch=h,
               linewidth=1.0, label=lab, zorder=3)
    ax.set_xticks(np.arange(5))
    ax.set_xticklabels(ticks)
    ax.set_xlim(-0.5, 4.5)
    ax.set_ylim(0, 310)
    ax.set_ylabel("Standard deviation of tracking error (mHz)")
    ax.grid(True, axis="y", alpha=0.35, zorder=1)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.13), ncol=4, frameon=False, fontsize=10)
    fig.tight_layout()
    p = os.path.join(outdir, "Figure_14_region_dispersion.png")
    fig.savefig(p, dpi=DPI)
    plt.close(fig)
    print(f"  wrote {p}")


# ---------------------------------------------------------------------------
# Figure 18 - swing-equation inversion on external inputs (f0 and dP published,
# the initial RoCoF nominal; see Section 3.11)
# ---------------------------------------------------------------------------

def fig18(outdir):
    display = {"Initial trip + embedded loss": "Cumulative after embedded loss (1,481 MW)",
               "Cumulative to LFDD": "Final cumulative total (1,878 MW)",
               "Reference (900 MW)": "Illustrative case (900 MW)"}
    rocof = IE.PUBLISHED_ROCOF
    ref = IE.REF_INERTIA_GVAS
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    axis = [0.05 + 0.01 * k for k in range(36)]
    for label, dp in IE.SCENARIOS:
        ax1.plot(axis, [IE.recover_ek_gvas(dp, r) for r in axis], linewidth=1.6, label=display[label])
    ax1.axhline(ref, color="#333333", linestyle="--", linewidth=1.0,
                label=f"Published pre-event inertia = {ref:.0f} GVA·s")
    rec = IE.recover_ek_gvas(IE.SCENARIOS[0][1], rocof)
    ax1.plot([rocof], [rec], marker="*", markersize=16, color="#d62728", linestyle="none",
             label=f"Recovered at nominal {rocof:.2f} Hz/s = {rec:.1f} GVA·s")
    ax1.set_ylim(0, 700)
    ax1.set_xlabel("|RoCoF| (Hz/s)")
    ax1.set_ylabel("Recovered stored kinetic energy $E_{\\mathrm{k}}$ (GVA·s)")
    ax1.set_title("(a) Swing-equation inversion on external inputs")
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8.5, loc="upper right")
    errs = list(range(-30, 31, 2))
    rec_err = [(IE.recover_ek_gvas(IE.SCENARIOS[0][1], rocof * (1 + ep / 100.0)) - ref) / ref * 100 for ep in errs]
    ax2.plot(errs, rec_err, color="#1f77b4", linewidth=1.6)
    ax2.axhspan(-5, 5, alpha=0.12, color="green")
    ax2.text(-29, 5.5, "±5 % target band", fontsize=9, va="bottom", ha="left", color="#2d6b2d", style="italic")
    base = (rec - ref) / ref * 100
    ax2.plot([0.0], [base], marker="*", markersize=16, color="#d62728", linestyle="none",
             label=f"Nominal {rocof:.2f} Hz/s input: {base:+.1f} %")
    ax2.set_xlabel("RoCoF measurement error (%)")
    ax2.set_ylabel("Inertia recovery error (%)")
    ax2.set_title("(b) RoCoF accuracy requirement")
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9, loc="upper right")
    fig.tight_layout()
    p = os.path.join(outdir, "Figure_18_inertia.png")
    fig.savefig(p, dpi=DPI)
    plt.close(fig)
    print(f"  wrote {p}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=os.path.join(ROOT, "figures"))
    ap.add_argument("--dma", default=None, help="raw V7.2 DMA capture (available on request)")
    ap.add_argument("--titles", action="store_true",
                    help="draw the working title (capture file and statistics) on Figures 11 and 15-17")
    a = ap.parse_args(argv)
    global SHOW_TITLES
    SHOW_TITLES = a.titles
    os.makedirs(a.outdir, exist_ok=True)
    plt.rcParams.update({"font.size": 10.5})
    fig11(a.outdir)
    fig12(a.outdir)
    fig13(a.outdir)
    fig14(a.outdir)
    fig15(a.outdir)
    if a.dma:
        fig16(a.outdir, a.dma)
    else:
        print("  Figure 16 skipped: pass --dma with the raw V7.2 capture to regenerate it")
    fig17(a.outdir)
    fig18(a.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
