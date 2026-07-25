"""
kalman_freq.py
==============

V4 Kalman frequency estimator for Section 5.3.11 (Table 16, V4 row) of the
open-RTU MDPI Sensors paper.

The V4 algorithm is a two-state linear Kalman filter run over the V2
processed-frequency stream (the f_Vb column of the canonical replay). It is
NOT a merge of V1 and V2 - that is the V3 hybrid (hybrid_freq_analysis.py).
V4 smooths the single V2 stream:

    state    x = [ f, RoCoF ]^T
    process  constant-velocity:  f_{k+1} = f_k + RoCoF_k * dt
                                 RoCoF_{k+1} = RoCoF_k    (+ process noise)
    measure  z = f  (the V2 estimate),  H = [1, 0]

    F = [[1, dt],        Q = q * [[dt^4/4, dt^3/2],
         [0,  1]]                 [dt^3/2, dt^2  ]]     (continuous-white-
                                                         noise-acceleration
                                                         discretisation)
    R = [r]

-------------------------------------------------------------------------
PROVENANCE NOTE - read before quoting these numbers.

The Table 16 V4 row in Draft V8 was an ESTIMATE inserted ahead of running
this filter. This script REPLACES that estimate with the filter's true
output on the real V2 data. The tuning constants Q and R are NOT recovered
from a prior run (none existed); they are DERIVED FROM FIRST PRINCIPLES so
the result is defensible and reproducible rather than fitted to a target:

  R (measurement-noise variance) is estimated directly from the V2
  detector's own steady-state scatter: the pooled variance of the V2
  frequency error over the two quiet regions of the replay (pre-event
  nominal, indices 0-45, and post-event settled, indices 301-359). This is
  the honest measurement noise the filter actually sees.

  Q (process-noise intensity) is set from the plausible frame-to-frame
  change in true grid RoCoF. We use the variance of the first-difference of
  the NESO reference RoCoF as the process-noise intensity q. This encodes
  "how much can the true RoCoF change between frames" with no free knob.

Both are printed and stored in the JSON so a reviewer can see exactly where
they came from. Whatever std/RMSE the filter then produces IS the V4 result;
Table 16 should be updated to match this script's output, not the reverse.
-------------------------------------------------------------------------

Usage:
    python kalman_freq.py replay_20260627_102507V2.csv
    python kalman_freq.py replay_20260627_102507V2.csv --neso neso_trajectory.csv
    python kalman_freq.py replay_20260627_102507V2.csv --q 1e-4 --r 4e-4  (manual override)

Writes kalman_output.json (+ optional overlay PNG).
Standard library only for the numerics; matplotlib optional for the figure.

Author: Jack Davies
"""

from __future__ import annotations
import argparse, csv, json, math, os, statistics


# ============================================================================
# Constants
# ============================================================================

DT_S = 1.0        # one replay frame == one NESO index-second
QUIET_REGIONS = [(0, 45), (301, 359)]   # pre-event nominal + post-event settled


# ============================================================================
# Data loading
# ============================================================================

def load_capture(path: str):
    """Return list of dicts sorted by replay_index (>=0 only)."""
    rows = []
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            try:
                i = int(r["replay_index"])
                if i < 0:
                    continue
                rows.append({"i": i,
                             "f_Vb": float(r["f_Vb"]),
                             "f_target": float(r["f_target_hz"])})
            except (KeyError, ValueError):
                continue
    rows.sort(key=lambda d: d["i"])
    return rows


def load_neso(path: str):
    out = {}
    if not os.path.isfile(path):
        return out
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            try:
                out[int(r["index"])] = float(r["frequency_hz"])
            except (KeyError, ValueError):
                continue
    return out


# ============================================================================
# First-principles Q and R derivation
# ============================================================================

def derive_R(rows) -> float:
    """Measurement-noise variance from V2 steady-state scatter (Hz^2).
    Pooled variance of (f_Vb - f_target) over the two quiet regions."""
    resid = []
    for d in rows:
        if any(lo <= d["i"] <= hi for lo, hi in QUIET_REGIONS):
            resid.append(d["f_Vb"] - d["f_target"])
    if len(resid) < 2:
        # Fall back to full-window residual variance
        resid = [d["f_Vb"] - d["f_target"] for d in rows]
    return statistics.pvariance(resid)


def derive_Q(rows, neso) -> float:
    """Process-noise intensity q from the variance of the frame-to-frame
    change in true RoCoF (Hz^2/s^4-ish intensity, applied via the CWNA Q
    matrix). Uses the NESO reference if available, else the emitted target."""
    # Build an index-ordered true-frequency series
    idx = [d["i"] for d in rows]
    if neso:
        f_true = [neso.get(d["i"], d["f_target"]) for d in rows]
    else:
        f_true = [d["f_target"] for d in rows]
    # First difference -> RoCoF; second difference -> RoCoF change
    rocof = [ (f_true[k] - f_true[k-1]) / DT_S for k in range(1, len(f_true)) ]
    d_rocof = [ (rocof[k] - rocof[k-1]) / DT_S for k in range(1, len(rocof)) ]
    if len(d_rocof) < 2:
        return 1.0e-6
    return statistics.pvariance(d_rocof)


# ============================================================================
# Two-state Kalman filter (constant-velocity / CWNA)
# ============================================================================

def kalman_run(rows, q: float, r: float):
    """Run the two-state KF over the f_Vb stream. Returns list of filtered f."""
    dt = DT_S
    # CWNA process-noise matrix scaled by intensity q
    q11 = q * dt**4 / 4.0
    q12 = q * dt**3 / 2.0
    q22 = q * dt**2
    # Initialise
    f0 = rows[0]["f_Vb"]
    x0, x1 = f0, 0.0                      # [f, rocof]
    # Covariance init - generous
    p00, p01, p10, p11 = 1.0, 0.0, 0.0, 1.0
    out = []
    for d in rows:
        z = d["f_Vb"]
        # --- Predict ---
        # x = F x
        x0p = x0 + dt * x1
        x1p = x1
        # P = F P F^T + Q
        # F = [[1, dt],[0,1]]
        fp00 = p00 + dt*p10 + dt*(p01 + dt*p11)
        fp01 = p01 + dt*p11
        fp10 = p10 + dt*p11
        fp11 = p11
        p00p = fp00 + q11
        p01p = fp01 + q12
        p10p = fp10 + q12
        p11p = fp11 + q22
        # --- Update (H = [1,0], scalar measurement) ---
        y = z - x0p                       # innovation
        s = p00p + r                      # innovation covariance
        k0 = p00p / s
        k1 = p10p / s
        x0 = x0p + k0 * y
        x1 = x1p + k1 * y
        # P = (I - K H) P
        p00 = (1 - k0) * p00p
        p01 = (1 - k0) * p01p
        p10 = p10p - k1 * p00p
        p11 = p11p - k1 * p01p
        out.append(x0)
    return out


def error_stats(filtered, rows, ref="f_target"):
    errs = [(filtered[k] - rows[k][ref]) * 1000.0 for k in range(len(rows))]
    n = len(errs)
    mean = statistics.mean(errs)
    std = statistics.stdev(errs)
    rmse = math.sqrt(sum(e*e for e in errs) / n)
    mx = max(abs(e) for e in errs)
    return {"n": n, "mean_mhz": mean, "std_mhz": std,
            "rmse_mhz": rmse, "max_mhz": mx}, errs


# ============================================================================
# Main
# ============================================================================

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", help="canonical V2 replay CSV")
    ap.add_argument("--neso", default="neso_trajectory.csv",
                    help="NESO reference CSV (for Q derivation & error baseline)")
    ap.add_argument("--q", type=float, default=None,
                    help="manual process-noise intensity override")
    ap.add_argument("--r", type=float, default=None,
                    help="manual measurement-noise variance override")
    ap.add_argument("--json", default="kalman_output.json")
    ap.add_argument("--figure", default="Figure_kalman_v4.png")
    ap.add_argument("--no-figure", action="store_true")
    args = ap.parse_args(argv)

    rows = load_capture(args.csv)
    if not rows:
        raise SystemExit(f"No replay rows in {args.csv}")
    neso = load_neso(args.neso)

    # Reference column for error: prefer NESO if aligned, else emitted target
    if neso:
        for d in rows:
            if d["i"] in neso:
                d["f_target"] = neso[d["i"]]

    r_val = args.r if args.r is not None else derive_R(rows)
    q_val = args.q if args.q is not None else derive_Q(rows, neso)

    filtered = kalman_run(rows, q_val, r_val)
    v4_stats, _ = error_stats(filtered, rows)

    # V2 baseline (unfiltered) for the improvement figure
    v2_filtered = [d["f_Vb"] for d in rows]
    v2_stats, _ = error_stats(v2_filtered, rows)

    improvement = (1.0 - v4_stats["std_mhz"] / v2_stats["std_mhz"]) * 100.0

    print()
    print("=" * 72)
    print("  V4 KALMAN FREQUENCY ESTIMATOR (Table 16, V4 row)")
    print("=" * 72)
    print(f"  Input CSV:            {os.path.basename(args.csv)}")
    print(f"  Frames:               {v4_stats['n']}")
    print()
    print("  First-principles tuning (derived, not fitted):")
    print(f"    R (meas-noise var,  from V2 quiet-region scatter): {r_val:.6e} Hz^2"
          f"  ({math.sqrt(r_val)*1000:.2f} mHz RMS)")
    print(f"    Q (process-noise intensity, from d(RoCoF) var):    {q_val:.6e}")
    print()
    print(f"  {'Algorithm':<14}{'mean [mHz]':>12}{'std [mHz]':>12}"
          f"{'RMSE [mHz]':>12}{'max [mHz]':>12}")
    print(f"  {'-'*14}{'-'*12}{'-'*12}{'-'*12}{'-'*12}")
    print(f"  {'V2 canonical':<14}{v2_stats['mean_mhz']:>+12.2f}"
          f"{v2_stats['std_mhz']:>12.2f}{v2_stats['rmse_mhz']:>12.2f}"
          f"{v2_stats['max_mhz']:>12.1f}")
    print(f"  {'V4 Kalman':<14}{v4_stats['mean_mhz']:>+12.2f}"
          f"{v4_stats['std_mhz']:>12.2f}{v4_stats['rmse_mhz']:>12.2f}"
          f"{v4_stats['max_mhz']:>12.1f}")
    print()
    print(f"  V4 improvement over V2 (std): {improvement:+.1f}%")
    print("=" * 72)
    print("  NOTE: this is the TRUE filter output. Update Table 16's V4 row to")
    print("        match these numbers (the previous row was an estimate).")
    print("=" * 72)

    out = {
        "schema_version": "1.0",
        "tool": "kalman_freq.py",
        "source_file": os.path.basename(args.csv),
        "provenance": "V4 row was an estimate in Draft V8; this is the true "
                      "filter output. Q and R derived from first principles.",
        "tuning": {
            "R_measurement_noise_var_hz2": r_val,
            "R_rms_mhz": math.sqrt(r_val) * 1000.0,
            "R_source": "pooled variance of V2 (f_Vb - reference) over quiet "
                        "regions 0-45 and 301-359",
            "Q_process_noise_intensity": q_val,
            "Q_source": "variance of frame-to-frame change in reference RoCoF",
            "dt_s": DT_S,
            "manual_override": (args.q is not None or args.r is not None),
        },
        "v2_stats": v2_stats,
        "v4_stats": v4_stats,
        "v4_improvement_over_v2_std_pct": improvement,
        # convenience block for crb_analysis.py --from-json chaining
        "achieved_std_mhz": {"V4 Kalman": v4_stats["std_mhz"]},
    }
    with open(args.json, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print(f"  -> {args.json}")

    if not args.no_figure:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            idx = [d["i"] for d in rows]
            ref = [d["f_target"] for d in rows]
            fig, ax = plt.subplots(figsize=(10, 5))
            ax.plot(idx, ref, color="#1f77b4", linewidth=1.6, label="Reference")
            ax.plot(idx, v2_filtered, color="#d62728", linewidth=0.8,
                    alpha=0.6, label="V2 canonical")
            ax.plot(idx, filtered, color="#2ca02c", linewidth=1.1,
                    label="V4 Kalman")
            ax.set_xlabel("Replay index (1 s per step)")
            ax.set_ylabel("Frequency (Hz)")
            ax.set_title("V4 Kalman filter over V2 stream vs reference")
            ax.grid(True, alpha=0.3)
            ax.legend(loc="lower right", fontsize=9)
            fig.tight_layout()
            fig.savefig(args.figure, dpi=200, bbox_inches="tight")
            plt.close(fig)
            print(f"  -> {args.figure}")
        except ImportError:
            pass

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
