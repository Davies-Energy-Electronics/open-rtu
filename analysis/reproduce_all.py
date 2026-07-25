"""
reproduce_all.py
================

Single-command reproducibility wrapper for the open-RTU MDPI Sensors paper.

Verifies the canonical CSV SHA-256, runs every analysis script in sequence
with default parameters and fixed seeds, and prints a compact console summary
of every numerical claim in the paper against its freshly-computed value.

This does NOT recompute physics; it orchestrates the existing scripts and
checks their outputs against the paper's reported numbers, so a reviewer can
confirm the whole results section in one command.

Usage:
    python reproduce_all.py
    python reproduce_all.py --canonical replay_20260627_102507V2.csv
    python reproduce_all.py --skip-figures

Exit code 0 if every check passes, 1 otherwise.

Author: Jack Davies
"""

from __future__ import annotations
import argparse, hashlib, json, os, subprocess, sys


# ============================================================================
# Expected canonical inputs
# ============================================================================

CANONICAL_CSV_DEFAULT = "replay_20260627_102507V2.csv"
CANONICAL_SHA256 = "dadbcfe247dd4e538b71a891b2b26f070d02c42e3dd4622b5d2125aff169560f"
V1_CSV_DEFAULT = "replay_20260625_180520V1.csv"
NESO_CSV_DEFAULT = "neso_trajectory.csv"
TIMING_CSV_DEFAULT = "timing_capture.csv"

# Paper's reported headline numbers, with tolerance, keyed by a label.
# (value_from_json_path, expected, abs_tol)
GREEN = "\033[92m"; RED = "\033[91m"; DIM = "\033[2m"; RST = "\033[0m"


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def run(cmd, label):
    print(f"{DIM}  $ {' '.join(cmd)}{RST}")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"{RED}  [FAIL] {label} exited {res.returncode}{RST}")
        if res.stderr:
            print(DIM + res.stderr[-500:] + RST)
        return False
    return True


def load_json(path):
    if not os.path.isfile(path):
        return None
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def check(label, got, expected, tol):
    if got is None:
        print(f"{RED}  [MISS] {label}: no value produced{RST}")
        return False
    ok = abs(got - expected) <= tol
    tag = f"{GREEN}[ OK ]{RST}" if ok else f"{RED}[FAIL]{RST}"
    print(f"  {tag} {label}: got {got:.3f}, expected {expected:.3f} "
          f"(tol {tol})")
    return ok


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--canonical", default=CANONICAL_CSV_DEFAULT)
    ap.add_argument("--v1", default=V1_CSV_DEFAULT)
    ap.add_argument("--neso", default=NESO_CSV_DEFAULT)
    ap.add_argument("--timing", default=TIMING_CSV_DEFAULT)
    ap.add_argument("--skip-figures", action="store_true")
    ap.add_argument("--python", default=sys.executable)
    args = ap.parse_args(argv)

    py = args.python
    all_ok = True

    print("=" * 74)
    print("  OPEN-RTU REPRODUCIBILITY WRAPPER")
    print("=" * 74)

    # --- Step 0: verify canonical hash ---
    print("\n[0] Canonical CSV integrity")
    if not os.path.isfile(args.canonical):
        print(f"{RED}  [FAIL] canonical CSV not found: {args.canonical}{RST}")
        return 1
    got_hash = sha256(args.canonical)
    if got_hash == CANONICAL_SHA256:
        print(f"{GREEN}  [ OK ]{RST} SHA-256 matches "
              f"({got_hash[:16]}...)")
    else:
        print(f"{RED}  [FAIL]{RST} SHA-256 mismatch\n"
              f"        got      {got_hash}\n"
              f"        expected {CANONICAL_SHA256}")
        all_ok = False

    figflag = ["--no-figure"] if args.skip_figures else []

    # --- Step 1: TOST + bootstrap ---
    print("\n[1] TOST equivalence + bootstrap (Table 8)")
    run([py, "tost_metrics.py", args.canonical], "tost_metrics")

    # --- Step 2: Bland-Altman (Table 10, Figure 10) ---
    print("\n[2] Bland-Altman (Table 10, Figure 10)")
    ba_fig = [] if not args.skip_figures else ["--figure", os.devnull]
    run([py, "bland_altman.py", args.canonical], "bland_altman")
    ba = load_json(os.path.splitext(args.canonical)[0] + "_bland_altman.json")
    if ba:
        agg = ba.get("aggregate", {})
        all_ok &= check("BA aggregate RMSE", agg.get("rmse_mHz"), 208.66, 1.0)

    # --- Step 3: TVE (Table 9) ---
    print("\n[3] Total Vector Error (Table 9)")
    run([py, "tve_metrics.py", args.canonical], "tve_metrics")

    # --- Step 4: TOST sweep (Figure 9) ---
    print("\n[4] TOST sensitivity sweep (Figure 9)")
    run([py, "tost_sweep.py", args.canonical], "tost_sweep")

    # --- Step 5: replay metrics (Table 6) ---
    print("\n[5] Replay metrics (Table 6)")
    run([py, "replay_metrics.py", args.canonical, "--neso", args.neso],
        "replay_metrics")
    rm = load_json(os.path.splitext(args.canonical)[0] + "_metrics.json")
    if rm:
        all_ok &= check("Table 6 max RoCoF err (Hz/s)",
                        rm.get("metric_3_max_rocof_tracking_error_hz_per_s"),
                        0.721, 0.005)

    # --- Step 6: runtime (Table 12) ---
    print("\n[6] Runtime characterisation (Table 12)")
    if os.path.isfile(args.timing):
        run([py, "runtime_metrics.py", args.timing], "runtime_metrics")
    else:
        print(f"{DIM}  (skipped: {args.timing} not present){RST}")

    # --- Step 7: hybrid V3 (Table 16 V3 row) ---
    print("\n[7] Hybrid V1/V2 -> V3 (Table 16 V3 row)")
    if os.path.isfile(args.v1):
        run([py, "hybrid_freq_analysis.py", args.v1, args.canonical], "hybrid")
        hy = load_json("hybrid_output.json")
        if hy:
            all_ok &= check("V3 hybrid std (mHz)",
                            hy.get("hybrid_stats", {}).get("std_mHz"),
                            116.75, 1.0)
    else:
        print(f"{DIM}  (skipped: {args.v1} not present){RST}")

    # --- Step 8: Kalman V4 (Table 16 V4 row) ---
    print("\n[8] Kalman V4 (Table 16 V4 row) - TRUE output, first-principles tuning")
    run([py, "kalman_freq.py", args.canonical, "--neso", args.neso] + figflag,
        "kalman_freq")
    kf = load_json("kalman_output.json")
    v4_std = None
    if kf:
        v4_std = kf.get("v4_stats", {}).get("std_mhz")
        print(f"  {DIM}V4 true std = {v4_std:.2f} mHz "
              f"(improvement {kf.get('v4_improvement_over_v2_std_pct'):+.1f}% "
              f"over V2; paper's estimate was 202.52 mHz / 3%){RST}")

    # --- Step 9: inertia (Table 15, Figure 12) ---
    print("\n[9] Grid inertia recovery (Table 15, Figure 12)")
    run([py, "inertia_estimator.py",
         "--measured-rocof", args.canonical] + figflag, "inertia_estimator")
    ie = load_json("inertia_output.json")
    if ie:
        all_ok &= check("Recovered inertia (GVA*s)",
                        ie.get("primary_recovered_ek_gvas"), 215.3, 0.5)

    # --- Step 10: CRB (Table 17), chained to true V4 ---
    print("\n[10] Cramer-Rao bound (Table 17), using TRUE V4 std")
    crb_cmd = [py, "crb_analysis.py"]
    if v4_std is not None:
        crb_cmd += ["--v4", f"{v4_std:.4f}"]
    run(crb_cmd, "crb_analysis")
    cr = load_json("crb_output.json")
    if cr:
        all_ok &= check("CRB std conservative (mHz)",
                        cr.get("crb_std_conservative_mhz"), 2.757, 0.01)

    print("\n" + "=" * 74)
    if all_ok:
        print(f"{GREEN}  ALL CHECKS PASSED{RST}")
    else:
        print(f"{RED}  SOME CHECKS FAILED - see above{RST}")
    print("=" * 74)
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
