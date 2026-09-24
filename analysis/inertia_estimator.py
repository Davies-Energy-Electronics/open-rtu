"""
inertia_estimator.py
====================

Swing-equation inversion for grid inertia, reproducing Figure 18 of
Section 3.11 (Grid-Inertia Recovery) of the open-RTU MDPI Sensors paper.

The swing equation in per-unit form relates the initial RoCoF following a
step power imbalance dP to the system's stored kinetic energy E_k:

    df/dt = -(f0 / (2 * E_k)) * dP_pu

Rearranged for E_k in GVA*s, with dP in MW and RoCoF in Hz/s:

    E_k [GVA*s] = | f0 * dP[MW] | / (2 * |df/dt[Hz/s]|) * 1e-3

WHAT THIS SCRIPT DOES AND DOES NOT ESTABLISH
--------------------------------------------
No term in the recovery (f0, dP, RoCoF) is a measurement made by this
instrument. f0 and dP are published figures; the initial RoCoF of 0.16 Hz/s
is a nominal figure for the event, not a published measurement (see the
provenance block below). The 231.4 GVA*s recovered at 1,481 MW, against the
published 210 GVA*s pre-event inertia, is therefore a worked instance of the
sensitivity in panel (b): the 0.16 Hz/s input sits about 10 % below the
0.176 Hz/s implied by the two published figures, and the recovered inertia
sits about 10 % above 210 GVA*s. It is NOT evidence about the instrument's
accuracy and the paper does not present it as such. The operative output is
the sensitivity analysis in panel (b) of Figure 18: it establishes the RoCoF
measurement accuracy a field deployment must achieve, which is a
specification for Gate 2 (sub-cycle cadence) and Gate 4 (real disturbance)
of the roadmap in Section 4.5.

PROVENANCE OF THE INPUTS - READ THIS BEFORE CITING ANY NUMBER BELOW
-------------------------------------------------------------------
The National Grid ESO Technical Report on the Events of 9 August 2019
(reference [2] of the manuscript as submitted; Table 2 of the report) gives
the generation losses as a running cumulative total:

    Hornsea One deload (799 MW -> 62 MW)            737 MW
    + Little Barford steam turbine trip             244 MW   -> 981 MW
    + embedded generation lost on vector shift     ~150 MW   -> 1,131 MW
    + embedded generation lost on RoCoF protection ~350 MW   -> 1,481 MW
    + Little Barford GT1A                           210 MW   -> 1,691 MW
    + Little Barford GT1B                           187 MW   -> 1,878 MW

  * dP. The primary scenario uses the published cumulative total of
    1,481 MW (after the embedded generation lost on RoCoF protection) and
    the second the published final total of 1,878 MW. Revisions before
    23 September 2026 used 1378 MW, which does not appear in the report: it
    was justified as 737 + 244 + 397 MW, whereas the report gives ~150 MW
    on vector shift and ~350 MW on RoCoF protection. For reference, the
    inversion at 0.16 Hz/s gives 153.3 GVA*s at 981 MW, 176.7 at 1,131 MW,
    231.4 at 1,481 MW and 293.4 at 1,878 MW, against the published
    210 GVA*s pre-event inertia.

  * RoCoF = 0.16 Hz/s - flagged in the code as PAPER INPUT. The report
    states no measured system-wide initial RoCoF for the event. Its only
    Hz/s figure is a protection threshold - "some parts of the system may
    have experienced a rate of change of frequency of 0.125 Hz/s or above"
    - which is not the same quantity. Section 3.11 of the manuscript
    therefore identifies 0.16 Hz/s as a nominal figure for the event and
    makes no claim that it is retrievable from the public record.

The reference inertia of 210 GVA*s IS stated in the report (Table 4) and
needs no flag. Manuscript reference numbers change between revisions (this
file previously cited the report as [10], [45] or [46]), so the documents
are named here as well as numbered: [2] is the ESO technical report and [3]
the ESO interim report in the manuscript as submitted.

WHY THE REPLAY'S OWN SLOPE IS NOT USED
---------------------------------------
The reason is temporal resolution and window length, not time dilation. The
canonical bench replay runs in real time (0.9997 s per replay index, verified
against replay_20260627_102507V2.csv). What differs is the averaging window:

  * The nominal 0.16 Hz/s is an INITIAL RoCoF, conventionally evaluated
    over a sub-second window at the instant of the step.
  * The diagnostic this script prints is a least-squares MEAN slope over the
    whole 74-second descent region (indices 46-119). Measured on the
    canonical capture it is 0.00913 Hz/s - a factor of 17.5 below the
    nominal initial value, not the "factor of nearly nine" stated in
    earlier revisions of this file and in Section 3.11 of the manuscript.
    That earlier figure of 0.018 Hz/s is not what this script produces and
    never was; inertia_output.json has reported 0.0091 since the deposit.

  * Most of that factor is window length rather than record smoothing. Over
    the steepest 10-second window the same capture gives 0.1156 Hz/s on the
    measured channel (at index 50) and 0.0839 Hz/s on the reference
    trajectory (at index 43) - within a factor of 1.4 of the nominal
    initial value. The one-second resolution of the source record still
    prevents a true sub-second initial RoCoF from being formed, which is why
    the inversion uses external inputs; but the 17.5x figure is an artefact
    of comparing a 74-second mean with an instantaneous value and should not
    be read as instrument error.

OUTPUT FILES
------------
  inertia_output.json   Written with LF line endings explicitly, so that its
                        SHA-256 is identical on Windows, macOS and Linux.
                        Before this change the file inherited the platform's
                        line endings and the hash quoted in the Data
                        Availability Statement reproduced on Windows only.
  Figure_12_inertia.png The default filename is historical: this file renders
                        FIGURE 18 of the paper. The name is retained because
                        it is the name under which the artefact is deposited
                        and hashed. Note that a rendered PNG is NOT
                        byte-reproducible across matplotlib and freetype
                        versions - bbox_inches="tight" crops to font metrics,
                        so the image dimensions themselves differ between
                        environments. Any SHA-256 quoted for it can only be
                        reproduced on the machine that made it.

The JSON key "table_15_rows" is retained for wire compatibility with the
deposited sidecar. It refers to the three-scenario table printed on the
console; Section 3.11 of the current manuscript carries no numbered table.

Usage:
    python inertia_estimator.py
    python inertia_estimator.py --rocof 0.16 --ref-inertia 210
    python inertia_estimator.py --figure Figure_12_inertia.png
    python inertia_estimator.py --measured-rocof replay_20260627_102507V2.csv

Reproducibility (re-run against the canonical capture, release v1.2.0):
    Recovered E_k at 1481 MW, 0.16 Hz/s:      231.4 GVA*s  (+10.2% vs 210)
    Recovered E_k at 1878 MW, 0.16 Hz/s:      293.4 GVA*s  (+39.7%)
    Recovered E_k at  900 MW, 0.16 Hz/s:      140.6 GVA*s  (-33.0%)
    Replay descent-region mean slope:         0.0091 Hz/s
    inertia_output.json SHA-256 (LF): quoted in the Data Availability
    Statement and checked by reproduce_all.py (DERIVED_HASHES).

Standard library only for the numerics; matplotlib only for --figure.

Author: Jack Davies
"""

from __future__ import annotations
import argparse, csv, json, math, os


# ============================================================================
# Constants (paper defaults)
# ============================================================================

F0_HZ            = 50.0

# PAPER INPUT - nominal initial-descent RoCoF for the event. The ESO
# technical report (manuscript ref. [2]) states no measured system-wide RoCoF;
# its only Hz/s figure is the 0.125 Hz/s protection threshold, a different
# quantity. Section 3.11 identifies 0.16 Hz/s as nominal. (The constant's
# name is historical and is kept because paper_figures.py imports it.)
PUBLISHED_ROCOF  = 0.16      # Hz/s

# Stated in the ESO technical report (manuscript ref. [2], Table 4). Sourced.
REF_INERTIA_GVAS = 210.0     # GVA*s, published pre-event inertia

# Three generation-loss scenarios (Figure 18 panel (a); see the provenance
# block in the module docstring before citing any of these).
SCENARIOS = [
    # CORRECTED 23 September 2026. Was 1378 MW, which does not appear in the
    # ESO technical report; its published cumulative totals are 981 / 1,131 /
    # 1,481 / 1,691 / 1,878 MW. 1,481 MW is the cumulative loss after the embedded generation
    # lost on RoCoF protection - the same quantity the old figure named, and
    # sourced. It is also consistent with Section 1, which already reports
    # ~500 MW of distribution-connected loss and 1,878 MW cumulative.
    ("Initial trip + embedded loss", 1481.0),
    # 1,878 MW is the published final cumulative loss in the report (after
    # Little Barford GT1A and GT1B). The label "to LFDD" is the paper's.
    ("Cumulative to LFDD", 1878.0),
    ("Reference (900 MW)",  900.0),    # smaller reference case, illustrative
]


# ============================================================================
# Swing equation
# ============================================================================

def recover_ek_gvas(dp_mw: float, rocof_hz_s: float, f0: float = F0_HZ) -> float:
    """Recovered stored kinetic energy E_k in GVA*s."""
    return abs(f0 * dp_mw) / (2.0 * abs(rocof_hz_s)) * 1.0e-3


def descent_rocof_from_csv(path: str,
                           i_lo: int = 46, i_hi: int = 119) -> float | None:
    """Least-squares MEAN slope of f_Vb over the descent region
    (indices 46-119), in Hz per replay-index-second, returned as a positive
    magnitude. Diagnostic only - NOT used for the recovery.

    This is a 74-second mean, not an initial RoCoF, so it is not comparable
    with the published 0.16 Hz/s without that caveat. See the module
    docstring."""
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
# Figure 18 (default file name Figure_12_inertia.png is historical)
# ============================================================================

def render_figure(rocof: float, ref_inertia: float, recovered_primary: float,
                  out_path: str, dpi: int = 600) -> bool:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[figure] matplotlib not installed - skipping Figure 18.")
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
    ax1.set_title("(a) Swing-equation inversion on published inputs")
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
    ax2.set_title("(b) RoCoF accuracy requirement")
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=8, loc="upper left")

    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
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
    ap.add_argument("--dpi", type=int, default=600,
                    help="figure resolution; 600 is the MDPI minimum for line art (default 600)")
    ap.add_argument("--json", type=str, default="inertia_output.json")
    ap.add_argument("--no-figure", action="store_true")
    args = ap.parse_args(argv)

    print()
    print("=" * 74)
    print("  SWING-EQUATION INERTIA INVERSION (Section 3.11, Figure 18)")
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
            print(f"\n  [diagnostic] replay descent-region MEAN slope over "
                  f"indices 46-119: {m:.4f} Hz/s")
            print(f"               Not used for the recovery, and not "
                  f"comparable with the published {args.rocof:.2f} Hz/s "
                  f"initial RoCoF:")
            print(f"               this is a 74 s mean, that is a sub-second "
                  f"instantaneous value. See Section 3.11.")

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
        "rocof_provenance": ("nominal initial RoCoF for the event; NOT measured by "
                             "the instrument and not a published measurement: the "
                             "ESO technical report states no measured system RoCoF "
                             "(its only Hz/s figure is the 0.125 Hz/s protection "
                             "threshold). See Section 3.11"),
        "dp_provenance": ("published cumulative losses in the ESO technical report: "
                          "981 / 1131 / 1481 / 1691 / 1878 MW; the primary scenario "
                          "uses 1481 MW (after the embedded generation lost on "
                          "RoCoF protection)"),
        "claim_scope": ("consistency check on published figures and on this "
                        "implementation; not evidence of instrument accuracy"),
    }
    # newline="\n" is deliberate: without it Python translates to CRLF on
    # Windows and the file's SHA-256 differs between platforms, so the hash
    # quoted in the Data Availability Statement could only ever reproduce on
    # the machine that wrote it. Verified 22 September 2026.
    with open(args.json, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(out, fh, indent=2)
    print(f"\n  -> {args.json}")

    if not args.no_figure:
        if render_figure(args.rocof, args.ref_inertia, recovered_primary, args.figure, args.dpi):
            print(f"  -> {args.figure}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
