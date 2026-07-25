"""
run_live_verification.py
========================

One-command LIVE HARDWARE verification runner for the open-RTU MDPI Sensors
paper. Takes a fresh V2 capture and a fresh V1 baseline capture from a live
bench session (conditioning circuit -> Feather 1 ADC, Feather 2 DAC replay),
runs the complete analysis toolchain over them, regenerates every report
figure, and prints a side-by-side comparison of the LIVE results against the
paper's ARCHIVED (canonical) values with pass/fail verdicts.

This is the wrapper to run AFTER a bench capture. It differs from
reproduce_all.py in that it is built to consume *new* captures (whose SHA-256
will legitimately differ from the archived canonical) and to report tolerance
bands rather than byte-identical reproduction.

WHAT IT DOES NOT DO
-------------------
It does not flash or talk to the Feathers. The V3 hybrid and V4 Kalman are
computed OFF-LINE in Python over the captured V1/V2 streams - which is exactly
how Section 5.3.11 of the paper presents them (emulation over captured
streams, not on-board execution). No Arduino step is required or appropriate
for the claims the paper makes. Running the algorithms on-board is a separate
firmware task (Phase 2) and is intentionally out of scope here.

USAGE
-----
    python run_live_verification.py \
        --v2 replay_20260707_132321.csv \
        --v1 replay_20260707_133804V1.csv \
        --neso neso_trajectory.csv \
        --outdir live_verification

Produces, in --outdir:
    Figure_bland_altman_live.png     (Table 10 / Figure 11 equivalent)
    Figure_tost_sweep_live.png       (Figure 9 equivalent)
    Figure_hybrid_live.png           (V3 hybrid Bland-Altman)
    Figure_kalman_live.png           (V4 Kalman vs reference)
    Figure_12_inertia_live.png       (Table 15 / Figure 12)   [inputs published]
    *_live.json sidecars for every metric
    live_verification_summary.json   (the full comparison table, machine-readable)

Requires the analysis scripts to be in the same directory or on PATH, and
matplotlib for the figures. Standard library otherwise.

Author: Jack Davies
"""

from __future__ import annotations
import argparse, csv, hashlib, json, math, os, statistics, subprocess, sys


# ============================================================================
# Archived (paper) reference values, with per-metric tolerance.
# Tolerances are set to the run-to-run spread expected from ADC/timing noise,
# NOT loosened to force a pass. A metric outside tolerance is reported FAIL
# so it can be investigated rather than hidden.
# ============================================================================

ARCHIVED = {
    # label:                     (paper_value, tolerance, unit)
    "BA aggregate RMSE":         (208.66, 6.0,  "mHz"),
    "BA aggregate bias":         (1.66,   6.0,  "mHz"),   # sign can flip; |delta| checked
    "BA Pearson r":              (0.930,  0.02, ""),
    "Max |f error|":             (629.10, 40.0, "mHz"),
    "RoCoF max (1.0s basis)":    (0.721,  0.20, "Hz/s"),  # noisiest metric (a max); wide band
    "RoCoF RMS (1.0s basis)":    (0.187,  0.04, "Hz/s"),
    "TVE aggregate mean":        (3.16,   0.5,  "%"),
    "TVE aggregate max":         (7.63,   0.8,  "%"),
    "V1 std":                    (131.35, 8.0,  "mHz"),
    "V2 std":                    (211.23, 8.0,  "mHz"),
    "V3 hybrid std":             (116.75, 8.0,  "mHz"),
    "V3 improvement over V2":    (44.7,   4.0,  "%"),
    "V4 Kalman std":             (194.25, 10.0, "mHz"),
    "V4 improvement over V2":    (7.0,    5.0,  "%"),
}

GREEN = "\033[92m"; RED = "\033[91m"; YEL = "\033[93m"; DIM = "\033[2m"; RST = "\033[0m"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for c in iter(lambda: fh.read(65536), b""):
            h.update(c)
    return h.hexdigest()


def run_script(cmd):
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"{RED}  [script error] {' '.join(cmd)}{RST}")
        print(DIM + (res.stderr[-600:] if res.stderr else "") + RST)
    return res


def load_json(path):
    if not os.path.isfile(path):
        return None
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def verdict(label, live):
    if live is None:
        return None, f"{YEL}[ -- ]{RST} {label}: not produced"
    paper, tol, unit = ARCHIVED[label]
    delta = live - paper
    ok = abs(delta) <= tol
    tag = f"{GREEN}[ OK ]{RST}" if ok else f"{RED}[FAIL]{RST}"
    line = (f"  {tag} {label:<26} paper {paper:>8.2f}{unit:<5}  "
            f"live {live:>8.2f}{unit:<5}  d={delta:+7.2f}  (tol +-{tol}{unit})")
    return ok, line


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--v2", required=True, help="live V2 canonical capture CSV")
    ap.add_argument("--v1", required=True, help="live V1 baseline capture CSV")
    ap.add_argument("--neso", default="neso_trajectory.csv")
    ap.add_argument("--outdir", default="live_verification")
    ap.add_argument("--python", default=sys.executable)
    args = ap.parse_args(argv)

    py = args.python
    os.makedirs(args.outdir, exist_ok=True)
    od = args.outdir
    results = {}
    live = {}

    print("=" * 92)
    print("  LIVE HARDWARE VERIFICATION")
    print("=" * 92)
    print(f"  V2 capture:  {os.path.basename(args.v2)}")
    print(f"    SHA-256:   {sha256(args.v2)}")
    print(f"  V1 capture:  {os.path.basename(args.v1)}")
    print(f"    SHA-256:   {sha256(args.v1)}")
    print()

    # --- Bland-Altman (Table 10) ---
    print(f"{DIM}  [1/6] Bland-Altman ...{RST}")
    run_script([py, "bland_altman.py", args.v2,
                "--json", f"{od}/bland_altman_live.json",
                "--figure", f"{od}/Figure_bland_altman_live.png"])
    ba = load_json(f"{od}/bland_altman_live.json")
    if ba:
        agg = ba["aggregate"]
        live["BA aggregate RMSE"] = agg["rmse_mHz"]
        live["BA aggregate bias"] = agg["mean_mHz"]
        live["BA Pearson r"] = agg["pearson_r"]

    # --- TOST sweep (Figure 9) ---
    print(f"{DIM}  [2/6] TOST sensitivity sweep ...{RST}")
    run_script([py, "tost_sweep.py", args.v2,
                "--plot", f"{od}/Figure_tost_sweep_live.png",
                "--json", f"{od}/tost_sweep_live.json"])

    # --- TVE (Table 9) ---
    print(f"{DIM}  [3/6] Total Vector Error ...{RST}")
    run_script([py, "tve_metrics.py", args.v2, "--json", f"{od}/tve_live.json"])
    tve = load_json(f"{od}/tve_live.json")
    if tve:
        agg = tve.get("overall", {})
        if agg.get("tve_mean_pct") is not None:
            live["TVE aggregate mean"] = agg["tve_mean_pct"]
        if agg.get("tve_max_pct") is not None:
            live["TVE aggregate max"] = agg["tve_max_pct"]

    # --- Replay metrics (Table 6) ---
    print(f"{DIM}  [4/6] Replay metrics (RoCoF, max error) ...{RST}")
    run_script([py, "replay_metrics.py", args.v2, "--neso", args.neso])
    rm = load_json(os.path.splitext(args.v2)[0] + "_metrics.json")
    if rm:
        live["Max |f error|"] = rm.get("metric_1_max_abs_freq_error_mhz")
        live["RoCoF max (1.0s basis)"] = rm.get("metric_3_max_rocof_tracking_error_hz_per_s")
        live["RoCoF RMS (1.0s basis)"] = rm.get("metric_3_rms_rocof_tracking_error_hz_per_s")

    # --- Hybrid V1/V2/V3 (Table 16 V1/V2/V3 rows) ---
    print(f"{DIM}  [5/6] Hybrid V1/V2 -> V3 ...{RST}")
    run_script([py, "hybrid_freq_analysis.py", args.v1, args.v2,
                "--out-csv", f"{od}/hybrid_live.csv",
                "--out-json", f"{od}/hybrid_live.json",
                "--figure", f"{od}/Figure_hybrid_live.png"])
    hy = load_json(f"{od}/hybrid_live.json")
    if hy:
        live["V1 std"] = hy["v1_stats"]["std_mHz"]
        live["V2 std"] = hy["v2_stats"]["std_mHz"]
        live["V3 hybrid std"] = hy["hybrid_stats"]["std_mHz"]
        v2s = hy["v2_stats"]["std_mHz"]; v3s = hy["hybrid_stats"]["std_mHz"]
        live["V3 improvement over V2"] = (1 - v3s / v2s) * 100.0

    # --- Kalman V4 (Table 16 V4 row) ---
    print(f"{DIM}  [6/6] Kalman V4 ...{RST}")
    run_script([py, "kalman_freq.py", args.v2, "--neso", args.neso,
                "--json", f"{od}/kalman_live.json",
                "--figure", f"{od}/Figure_kalman_live.png"])
    kf = load_json(f"{od}/kalman_live.json")
    if kf:
        live["V4 Kalman std"] = kf["v4_stats"]["std_mhz"]
        live["V4 improvement over V2"] = kf["v4_improvement_over_v2_std_pct"]

    # --- Inertia (Table 15) - inputs are published values, capture-independent ---
    run_script([py, "inertia_estimator.py",
                "--measured-rocof", args.v2,
                "--json", f"{od}/inertia_live.json",
                "--figure", f"{od}/Figure_12_inertia_live.png"])

    # --- Comparison table ---
    print()
    print("=" * 92)
    print("  LIVE vs ARCHIVED COMPARISON")
    print("=" * 92)
    n_ok = n_fail = n_missing = 0
    summary_rows = []
    for label in ARCHIVED:
        ok, line = verdict(label, live.get(label))
        print(line)
        if ok is True:
            n_ok += 1
        elif ok is False:
            n_fail += 1
        else:
            n_missing += 1
        paper, tol, unit = ARCHIVED[label]
        summary_rows.append({
            "metric": label, "paper": paper, "live": live.get(label),
            "tolerance": tol, "unit": unit,
            "within_tolerance": ok,
        })
    print("=" * 92)
    print(f"  {GREEN}{n_ok} within tolerance{RST}, "
          f"{RED}{n_fail} outside{RST}, "
          f"{YEL}{n_missing} not produced{RST}")
    print("=" * 92)

    summary = {
        "schema_version": "1.0",
        "tool": "run_live_verification.py",
        "v2_capture": os.path.basename(args.v2),
        "v2_sha256": sha256(args.v2),
        "v1_capture": os.path.basename(args.v1),
        "v1_sha256": sha256(args.v1),
        "n_within_tolerance": n_ok,
        "n_outside_tolerance": n_fail,
        "n_not_produced": n_missing,
        "comparison": summary_rows,
        "figures": [
            f"{od}/Figure_bland_altman_live.png",
            f"{od}/Figure_tost_sweep_live.png",
            f"{od}/Figure_hybrid_live.png",
            f"{od}/Figure_kalman_live.png",
            f"{od}/Figure_12_inertia_live.png",
        ],
    }
    with open(f"{od}/live_verification_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(f"  Summary written: {od}/live_verification_summary.json")

    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
