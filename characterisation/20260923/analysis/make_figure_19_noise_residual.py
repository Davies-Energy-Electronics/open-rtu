#!/usr/bin/env python3
"""
make_figure_19_noise_residual.py - Figure 19 (numbered 20 in earlier drafts), direct measurement of the acquisition chain,
regenerated from the released records only.

(a) Deviation from the record mean, in ADC codes, for the three quiet-input
    records of the 13 September bundle (20260913/A1_quiet_single.txt,
    A2_quiet_single_run2.txt, A3_quiet_single_run3.txt), as unit-width
    histograms (one bin per code), against the ideal uniform 12-bit quantiser
    (density 1 on [-0.5, 0.5] code, sigma_q = LSB/sqrt12 = 0.233 mV).  The
    legend quotes each record's RAW SD in mV (the values the text quotes,
    0.687 / 0.634 / 0.657; the detrended values are 0.687 / 0.633 / 0.655)
    and the ratio of the three-record mean to sigma_q.
(b) B1_ac_running.txt (13 Sep, production stimulus): residual after the
    single-tone fit in each 200-sample window (analyse_ac.py), the 81 residuals
    concatenated, one-sided power spectrum (rectangular window, power
    normalised); the plotted quantity is the rms amplitude per bin, sqrt(P).
    Harmonics H2-H5 are marked at h * f1.  The dashed line is the level a white
    noise of the quiet-input 0.659 mV would give per bin, 0.659*sqrt(2/N) mV.
    The annotation is the harmonic share of residual power computed here
    (+/-5 bins round every multiple of f1, h = 1..19).

Usage:  python3 make_figure_20.py [--out PATH]
"""
from __future__ import annotations
import argparse, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common_312 as c
import verify_stimulus_residual_312 as vs

COLS = ["#1f6fb4", "#e0651c", "#2a8466"]


def panel_a(ax):
    sds = []
    edges = np.arange(-4.5, 4.51, 1.0)
    for p, col, lab in zip(c.SEPT_QUIET, COLS, ("A1", "A2", "A3")):
        codes, _ = c.read_codes(p)
        x = codes.astype(float)
        raw, _ = c.sd_raw_and_detrended(x)
        sds.append(raw * c.LSB_V * 1e3)
        ax.hist(x - x.mean(), bins=edges, density=True, histtype="step", lw=2.2,
                color=col, label=f"{lab}  $\\sigma_v$ = {sds[-1]:.3f} mV")
    ax.plot([-0.5, -0.5, 0.5, 0.5], [0, 1, 1, 0], ls="--", lw=2, color="#555555",
            label=f"Ideal 12-bit quantisation\n($\\sigma_v$ = {c.SIGMA_Q_V*1e3:.3f} mV)")
    ratio = np.mean(sds) * 1e-3 / c.SIGMA_Q_V
    ax.text(0.03, 0.95, f"measured / ideal = {ratio:.1f}x\n"
            f"({12 - math.log2(ratio):.1f} effective bits)", transform=ax.transAxes,
            va="top", fontsize=8.5, color="#333333")
    ax.set_xlim(-5, 5)
    ax.set_ylim(0, 1.55)
    ax.set_xlabel("Deviation from mean (ADC codes)")
    ax.set_ylabel("Probability density")
    ax.set_title("(a)  Quiet input, three independent runs", fontsize=10.5)
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    return sds, ratio


def panel_b(ax):
    codes, meta = c.read_codes(c.B1)
    fs = float(meta["fs_hz"])
    v = codes.astype(float) * c.LSB_V
    f1 = c.refine_frequency(v[:4000] - v.mean(), np.arange(4000) / fs, 50.0)
    r, _, _ = vs.window_residuals(v, fs, f1)
    fr, P = vs.power_spectrum(r, fs)
    m = vs.harmonic_mask(fr, f1, 5)
    share = 100 * P[m].sum() / P.sum()
    amp = np.sqrt(P) * 1e3
    ax.semilogy(fr[1:], np.maximum(amp[1:], 1e-5), lw=0.6, color="#2f73b8")
    floor = 0.659 * math.sqrt(2.0 / len(r))
    ax.axhline(floor, ls="--", lw=1.4, color="#2a8466")
    ax.text(990, floor * 0.45, "quiet-input (0.659 mV) white-noise level", ha="right",
            va="top", fontsize=8, color="#2a8466",
            bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none", alpha=0.85))
    df = fr[1] - fr[0]
    for h in (2, 3, 4, 5):
        k = int(round(h * f1 / df))
        pk = amp[max(0, k - 5):k + 6].max()
        ax.plot(h * f1, pk * 1.6, marker="v", ms=8, color="#e0651c")
        ax.text(h * f1, pk * 3.2, f"H{h}", ha="center", va="bottom", fontsize=8.5,
                fontweight="bold", color="#e0651c")
    ax.text(0.97, 0.06, f"{share:.1f} % of residual power lies in\n"
            "harmonic bins: distortion, not noise", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8.5,
            bbox=dict(boxstyle="square,pad=0.4", fc="#f2f2f2", ec="none"))
    ax.set_xlim(0, 1000)
    ax.set_ylim(5e-5, 2e2)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Residual amplitude, rms per bin (mV)")
    ax.set_title("(b)  Stimulus present: residual after single-tone fit", fontsize=10.5)
    return share


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(c.FIG_DIR,
                                                  "Figure_19_noise_and_residual.png"))
    a = ap.parse_args(argv)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9.5})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 3.9), dpi=250)
    sds, ratio = panel_a(ax1)
    share = panel_b(ax2)
    for ax in (ax1, ax2):
        ax.grid(color="#dddddd", lw=0.7)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(a.out, bbox_inches="tight")
    print(f"  (a) raw SDs {', '.join(f'{s:.3f}' for s in sds)} mV; mean/ideal {ratio:.2f}x")
    print(f"  (b) harmonic share {share:.2f} %")
    print(f"  wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
