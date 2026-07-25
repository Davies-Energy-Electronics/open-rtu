"""
inertia_estimator.py
====================

Grid-inertia estimation from Rate of Change of Frequency (RoCoF) via the
swing equation, reproducing Table 15 and Figure 12 of Section 5.3.10 of the
open-RTU MDPI Sensors paper.

The swing equation in per-unit form relates the initial RoCoF following a
step power imbalance dP to the system's stored kinetic energy E_k:

    df/dt = -(f0 / (2 * E_k)) * dP_pu

Rearranged for E_k in GVA*s, with dP in MW and RoCoF in Hz/s:

    E_k [GVA*s] = | f0 * dP[MW] | / (2 * |df/dt[Hz/s]|) * 1e-3

For the NESO 9 August 2019 GB blackout, the published initial-descent RoCoF
was approximately 0.16 Hz/s and the published initial generation loss was
1378 MW (Hornsea 737 MW + Little Barford steam 244 MW + distributed-
generation trip 397 MW). Substituting recovers E_k = 215.3 GVA*s against a
published pre-event inertia of ~210 GVA*s, a +2.5% match.

IMPORTANT (matches the paper's own methodological caveat, Section 5.3.10):
The canonical bench replay is TIME-DILATED relative to the real event, so the
input RoCoF used here is the *published real-event* value (0.16 Hz/s), NOT the
replay's measured RoCoF. This script reproduces the algorithmic pathway on the
published inputs; a live-field deployment (Section 6.3 Gate 3) would supply a
directly-measured RoCoF instead. Pass --measured-rocof to see what the replay's
own descent-region RoCoF would recover, for comparison only.

Usage:
    python inertia_estimator.py
    python inertia_estimator.py --rocof 0.16 --ref-inertia 210
    python inertia_estimator.py --figure Figure_12_inertia.png
    python inertia_estimator.py --measured-rocof replay_20260627_102507V2.csv

Writes inertia_output.json (+ optional Figure 12 PNG if matplotlib present).
Standard library only for the numerics.

Author: Jack Davies
"""

from __future__ import annotations
import argparse, csv, json, math, os


# ============================================================================
# Constants (paper defaults)
# ============================================================================

F0_HZ            = 50.0
PUBLISHED_ROCOF  = 0.16      # Hz/s, published initial-descent RoCoF [45]
REF_INERTIA_GVAS = 210.0     # GVA*s, published pre-event NESO inertia [46]

# Three generation-loss scenarios (Table 15)
SCENARIOS = [
    ("Initial trip",       1378.0),   # Hornsea + Little Barford + DG trip
    ("Cumulative to LFDD",  1878.0),   # cumulative loss at LFDD activation
    ("Reference (900 MW)",  900.0),    # smaller reference case
]


# ============================================================================
# Swing equation
# ============================================================================

def recover_ek_gvas(dp_mw: float, rocof_hz_s: float, f0: float = F0_HZ) -> float:
    """Recovered stored kinetic energy E_k in GVA*s."""
    return abs(f0 * dp_mw) / (2.0 * abs(rocof_hz_s)) * 1.0e-3


def descent_rocof_from_csv(path: str,
                           i_lo: int = 46, i_hi: int = 119) -> float | None:
    """Least-squares slope of f_Vb over the descent region (indices 46-119),
    in Hz per replay-index-second. Returned as a positive magnitude.
    Reference/diagnostic only - NOT used for the Table 15 recovery."""
    xs, ys = [], []
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            try:
                i = int(r["replay_index"])
                if i_lo <= i <= i_hi:
                    xs.append(float(i))
                    ys.append(float(r["f_Vb"]))
            except (KeyError, ValueError):
                continue
    n = len(xs)
    if n < 2:
        return None
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    if den == 0:
        return None
    return abs(num / den)


# ============================================================================
# Figure 12
# ============================================================================

def render_figure(rocof: float, ref_inertia: float, recovered_primary: float,
                  out_path: str) -> bool:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[figure] matplotlib not installed - skipping Figure 12.")
        return False

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Panel (a): |RoCoF| -> E_k relationship at three loss levels
    rocof_axis = [r / 100.0 for r in range(5, 41)]  # 0.05 .. 0.40 Hz/s
    for label, dp in SCENARIOS:
        eks = [recover_ek_gvas(dp, r) for r in rocof_axis]
        ax1.plot(rocof_axis, eks, linewidth=1.6, label=f"{label} ({dp:.0f} MW)")
    ax1.axhline(ref_inertia, color="#333333", linestyle="--", linewidth=1.0,
                label=f"Published NESO = {ref_inertia:.0f} GVA\u00b7s")
    ax1.plot([rocof], [recovered_primary], marker="*", markersize=16,
             color="#d62728", linestyle="none",
             label=f"Recovered = {recovered_primary:.1f} GVA\u00b7s")
    ax1.set_xlabel("|RoCoF|  (Hz/s)")
    ax1.set_ylabel("Recovered stored kinetic energy  E$_k$  (GVA\u00b7s)")
    ax1.set_title("(a) Swing-equation recovery")
    ax1.set_ylim(0, 700)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8, loc="upper right")

    # Panel (b): recovery-error sensitivity to RoCoF measurement error
    dp_primary = SCENARIOS[0][1]
    err_pct_axis = [e for e in range(-30, 31, 2)]
    rec_err = []
    for ep in err_pct_axis:
        r_meas = rocof * (1.0 + ep / 100.0)
        ek = recover_ek_gvas(dp_primary, r_meas)
        rec_err.append((ek - ref_inertia) / ref_inertia * 100.0)
    ax2.plot(err_pct_axis, rec_err, color="#1f77b4", linewidth=1.6)
    ax2.axhspan(-5, 5, alpha=0.12, color="green")
    ax2.text(err_pct_axis[0], 5, " \u00b15% target band",
             fontsize=8, va="bottom", ha="left", color="#2d6b2d", style="italic")
    base_err = (recover_ek_gvas(dp_primary, rocof) - ref_inertia) / ref_inertia * 100.0
    ax2.plot([0.0], [base_err], marker="*", markersize=16, color="#d62728",
             linestyle="none", label=f"Achieved = {base_err:+.1f}%")
    ax2.set_xlabel("RoCoF measurement error  (%)")
    ax2.set_ylabel("Inertia recovery error  (%)")
    ax2.set_title("(b) Recovery-error sensitivity")
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=8, loc="upper left")

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return True


# ============================================================================
# Main
# ============================================================================

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rocof", type=float, default=PUBLISHED_ROCOF,
                    help=f"published initial-descent RoCoF Hz/s (default {PUBLISHED_ROCOF})")
    ap.add_argument("--ref-inertia", type=float, default=REF_INERTIA_GVAS,
                    help=f"published pre-event inertia GVA*s (default {REF_INERTIA_GVAS})")
    ap.add_argument("--measured-rocof", type=str, default=None,
                    help="Optional replay CSV: also print the replay's own "
                         "descent-region RoCoF for comparison (diagnostic only)")
    ap.add_argument("--figure", type=str, default="Figure_12_inertia.png")
    ap.add_argument("--json", type=str, default="inertia_output.json")
    ap.add_argument("--no-figure", action="store_true")
    args = ap.parse_args(argv)

    print()
    print("=" * 74)
    print("  GRID INERTIA ESTIMATION FROM ROCOF (Table 15)")
    print("=" * 74)
    print(f"  f0 = {F0_HZ:.1f} Hz,  published RoCoF = {args.rocof:.3f} Hz/s,  "
          f"reference inertia = {args.ref_inertia:.0f} GVA\u00b7s")
    print()
    print(f"  {'Scenario':<22}{'dP [MW]':>10}{'RoCoF [Hz/s]':>14}"
          f"{'E_k [GVA*s]':>14}{'Error':>10}")
    print(f"  {'-'*22}{'-'*10}{'-'*14}{'-'*14}{'-'*10}")

    rows = []
    for label, dp in SCENARIOS:
        ek = recover_ek_gvas(dp, args.rocof)
        err = (ek - args.ref_inertia) / args.ref_inertia * 100.0
        print(f"  {label:<22}{dp:>10.0f}{-args.rocof:>14.2f}"
              f"{ek:>14.1f}{err:>+9.1f}%")
        rows.append({"scenario": label, "dp_mw": dp,
                     "rocof_hz_s": -args.rocof,
                     "recovered_ek_gvas": round(ek, 1),
                     "error_vs_ref_pct": round(err, 1)})
    print("=" * 74)

    recovered_primary = recover_ek_gvas(SCENARIOS[0][1], args.rocof)

    measured = None
    if args.measured_rocof and os.path.isfile(args.measured_rocof):
        m = descent_rocof_from_csv(args.measured_rocof)
        if m is not None:
            measured = m
            print(f"\n  [diagnostic] replay descent-region measured RoCoF: "
                  f"{m:.4f} Hz/s  (NOT used for Table 15 - replay is "
                  f"time-dilated; see Section 5.3.10 caveat)")

    out = {
        "schema_version": "1.0",
        "tool": "inertia_estimator.py",
        "config": {"f0_hz": F0_HZ, "published_rocof_hz_s": args.rocof,
                   "ref_inertia_gvas": args.ref_inertia},
        "primary_recovered_ek_gvas": round(recovered_primary, 1),
        "primary_error_pct": round((recovered_primary - args.ref_inertia)
                                   / args.ref_inertia * 100.0, 1),
        "table_15_rows": rows,
        "replay_measured_descent_rocof_hz_s": (round(measured, 4)
                                               if measured is not None else None),
    }
    with open(args.json, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print(f"\n  -> {args.json}")

    if not args.no_figure:
        if render_figure(args.rocof, args.ref_inertia, recovered_primary, args.figure):
            print(f"  -> {args.figure}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
