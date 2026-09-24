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
TIMING_CSV_DEFAULT = "runtime_capture_20260725.csv"

# Derived analysis outputs whose SHA-256 is quoted in the Data Availability
# Statement. Regenerating them must reproduce these hashes exactly; if it does
# not, the deposit and the paper have diverged.
#
# CORRECTED 22 September 2026. The previous value for inertia_output.json,
# 96A91FE8...61DB4, was produced on Windows and could only ever reproduce
# there: the file inherited the platform's line endings. inertia_estimator.py
# now writes LF explicitly, so the hash below is the same on Windows, macOS
# and Linux.
# UPDATED for v1.2.0 (24 September 2026): the provenance strings inside
# inertia_output.json were corrected to match Section 3.11 (the 0.16 Hz/s
# input is nominal, not published; the 1378 MW note was stale). The numbers in
# the file are unchanged; the hash changes because the text does. The value
# below is the one quoted in the Data Availability Statement of the
# manuscript as submitted. v1.1.0 quoted 23EEF956...CAF3CC5.
DERIVED_HASHES = {
    "inertia_output.json":
        "91AFA54072CB60C38444D0223FC3AE2C41960CDB8703BD84F90D5DB055FCF9B2",
}

# Figure_12_inertia.png (which renders Figure 18 of the paper) is DELIBERATELY
# NOT hashed here, and the Data Availability Statement should not hash it
# either. A matplotlib PNG is not byte-reproducible across matplotlib and
# freetype versions: bbox_inches="tight" crops the canvas to measured font
# extents, so the image DIMENSIONS differ between environments - 2935x7134
# pixels on the authoring machine against 2926x7117 on a current Linux
# matplotlib 3.10. A hash promise that cannot be kept fails the run for a
# reader who has done nothing wrong, which is worse than making no promise.
# The figure is drawn entirely from inertia_output.json, which IS hashed
# above, so the underlying result is still pinned.

# Paper's reported headline numbers, with tolerance, keyed by a label.
# (value_from_json_path, expected, abs_tol)
GREEN = "\033[92m"; RED = "\033[91m"; DIM = "\033[2m"; RST = "\033[0m"


# ============================================================================
# Location independence  (24 September 2026)
# ============================================================================
# This wrapper used to resolve every script and data file as a bare name in the
# current working directory. That worked in the flat authoring folder and could
# never work in the published repository, whose scripts live in analysis/ and
# whose data lives in captures/ and trajectories/ - so the Data Availability
# Statement's claim that the released harness regenerates the paper was true of
# the working copy and had never been true of the deposit.
#
# Copying one file over the other would have fixed it once and let it diverge
# again on the next edit. Resolving paths instead makes ONE file correct in
# BOTH layouts, which is what actually closes the gap.
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = (os.path.dirname(_HERE)
         if os.path.isdir(os.path.join(os.path.dirname(_HERE), "analysis"))
         else _HERE)
_DATA_DIRS   = ("", "captures", "trajectories", "analysis", "figures")
_SCRIPT_DIRS = ("", "analysis")


def _search(name, dirs, bases):
    """First existing match, else the name unchanged so the caller's own
    error message still fires with something a reader recognises."""
    if os.path.isabs(name):
        return name
    for base in bases:
        for d in dirs:
            cand = os.path.normpath(os.path.join(base, d, name))
            if os.path.isfile(cand):
                return cand
    return name


def data(name):
    """Resolve an input or a derived output."""
    return _search(name, _DATA_DIRS, (os.getcwd(), _ROOT, _HERE))


def script(name):
    """Resolve a sibling analysis script."""
    return _search(name, _SCRIPT_DIRS, (_HERE, os.getcwd(), _ROOT))


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _child_env():
    """Force UTF-8 in child processes.

    Several analysis scripts print box-drawing characters. On Windows the
    default console encoding is cp1252, which cannot represent them, so the
    child raises UnicodeEncodeError and exits non-zero even though its
    numerics completed correctly. Setting PYTHONIOENCODING makes the run
    behave identically on Windows, macOS and Linux.
    """
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONUTF8"] = "1"
    return env


def run(cmd, label):
    print(f"{DIM}  $ {' '.join(cmd)}{RST}")
    res = subprocess.run(cmd, capture_output=True, text=True,
                         encoding="utf-8", errors="replace", env=_child_env())
    if res.returncode != 0:
        print(f"{RED}  [FAIL] {label} exited {res.returncode}{RST}")
        if res.stderr:
            print(DIM + res.stderr[-800:] + RST)
        return False
    return True


def load_json(path):
    if not os.path.isfile(path):
        return None
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def need_json(path, label):
    """Load a JSON output that the paper depends on.

    Returns (data, ok). A missing file is a FAILURE, not a skip: if the
    analysis script did not write its output, the checks that depend on it
    cannot run, and a wrapper that stayed silent would report success while
    verifying nothing. This was a real defect found on 8 August 2026.
    """
    data = load_json(path)
    if data is None:
        print(f"{RED}  [FAIL]{RST} {label}: expected output not found -> {path}\n"
              f"        the checks that depend on it did NOT run")
        return None, False
    return data, True


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
    # The wrapper prints box rules and a degree sign; make its own stream safe
    # on consoles that default to a legacy code page.
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--canonical", default=CANONICAL_CSV_DEFAULT)
    ap.add_argument("--v1", default=V1_CSV_DEFAULT)
    ap.add_argument("--neso", default=NESO_CSV_DEFAULT)
    ap.add_argument("--timing", default=TIMING_CSV_DEFAULT)
    ap.add_argument("--skip-figures", action="store_true")
    ap.add_argument("--python", default=sys.executable)
    args = ap.parse_args(argv)

    # Resolve inputs against both layouts before anything reads them.
    args.canonical = data(args.canonical)
    args.v1        = data(args.v1)
    args.neso      = data(args.neso)
    args.timing    = data(args.timing)

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

    # Delete the outputs this wrapper checks, so that a stale file from an
    # earlier run cannot be mistaken for a fresh computation. Found on
    # 8 August 2026: two scripts write JSON only when --json is passed, so the
    # wrapper had been validating files on disk rather than new results.
    stem = os.path.splitext(args.canonical)[0]
    for stale in (stem + "_tost.json", stem + "_bland_altman.json",
                  stem + "_tve.json", stem + "_tost_sweep.json",
                  stem + "_metrics.json", "kalman_output.json",
                  "inertia_output.json", "crb_output.json",
                  "hybrid_output.json"):
        stale = data(stale)
        if os.path.isfile(stale):
            os.remove(stale)
    print(f"{DIM}  (cleared previous analysis outputs so checks run on fresh results){RST}")

    figflag = ["--no-figure"] if args.skip_figures else []

    # --- Step 1: TOST + bootstrap (paper Table 2) ---
    print("\n[1] TOST equivalence + bootstrap (paper Section 3.2, Table 2)")
    # tost_metrics.py writes its sidecar automatically, named from the input
    # CSV, when --json is not given; this step relies on that default.
    # Verified 8 Aug 2026 -- its JSON self-reports the correct source
    # filename, unlike the two scripts that needed --json adding above.
    all_ok &= run([py, script("tost_metrics.py"), args.canonical], "tost_metrics")
    ts, _ok = need_json(os.path.splitext(args.canonical)[0] + "_tost.json",
                        "TOST metrics"); all_ok &= _ok
    if ts is None:
        print(f"{DIM}        tost_metrics.py writes <csv>_tost.json by default; it did not\n"
              f"        appear, so this step could not be verified.{RST}")
    if ts:
        CLASS_I = "Class I (\u00b110 mHz)"
        CLASS_II = "Class II (\u00b1100 mHz)"
        wins = ts.get("windows", {})
        for key, label, mean_mhz, sd_mhz in (
                ("pre_event_nominal", "pre-event nominal", -21.20, 20.4),
                ("post_event_settled", "post-event settled", -9.15, 39.4),
                ("full_replay_window", "full replay", 1.66, 208.9)):
            w2 = wins.get(key, {}).get(CLASS_II, {})
            w1 = wins.get(key, {}).get(CLASS_I, {})
            all_ok &= check(f"Sec 3.2 {label} mean (mHz)", w2.get("mean"), mean_mhz, 0.02)
            all_ok &= check(f"Sec 3.2 {label} sd (mHz)", w2.get("sd"), sd_mhz, 0.1)
            # The paper's claim is Class II accepted on every window, Class I
            # rejected on every window. Both directions must hold.
            ok = (w2.get("equivalent") is True) and (w1.get("equivalent") is False)
            tag = f"{GREEN}[ OK ]{RST}" if ok else f"{RED}[FAIL]{RST}"
            print(f"  {tag} Sec 3.2 {label}: Class II equivalent="
                  f"{w2.get('equivalent')}, Class I equivalent={w1.get('equivalent')}")
            all_ok &= ok

    # --- Step 2: Bland-Altman (paper Figure 13) ---
    # Under --skip-figures the script's own plot is sent to the null device, so
    # no Figure_10_Bland_Altman.png is written. Without --skip-figures the
    # script keeps its default and writes that file into the current working
    # directory. (In earlier revisions ba_fig was built but never passed, so the plot
    # was written on every run.) The manuscript's Figure 13 is drawn by
    # paper_figures.py in step 13.
    print("\n[2] Bland-Altman (paper Section 3.3, Figure 13)")
    ba_fig = [] if not args.skip_figures else ["--figure", os.devnull]
    all_ok &= run([py, script("bland_altman.py"), args.canonical,
                   "--json", os.path.splitext(args.canonical)[0] + "_bland_altman.json"]
                  + ba_fig,
                  "bland_altman")
    ba, _ok = need_json(os.path.splitext(args.canonical)[0] + "_bland_altman.json",
                        "Bland-Altman"); all_ok &= _ok
    if ba:
        agg = ba.get("aggregate", {})
        all_ok &= check("Sec 3.3 BA aggregate RMSE", agg.get("rmse_mHz"), 208.66, 1.0)

    # --- Step 3: TVE-equivalent (paper Table 3) ---
    print("\n[3] TVE-equivalent statistics (paper Section 3.4, Table 3)")
    all_ok &= run([py, script("tve_metrics.py"), args.canonical,
                   "--json", os.path.splitext(args.canonical)[0] + "_tve.json",
                   "--verify-hash", CANONICAL_SHA256], "tve_metrics")
    tv, _ok = need_json(os.path.splitext(args.canonical)[0] + "_tve.json",
                        "TVE metrics"); all_ok &= _ok
    if tv:
        ov = tv.get("overall", {})
        all_ok &= check("Table 3 aggregate TVE mean (%)", ov.get("tve_mean_pct"), 3.18, 0.02)
        all_ok &= check("Table 3 aggregate TVE max (%)", ov.get("tve_max_pct"), 8.72, 0.02)
        reg = tv.get("regions", {})
        for name, mean_pct in (("Pre-event nominal", 2.43),
                               ("RoCoF descent", 3.85),
                               ("Nadir + early recovery", 3.94)):
            all_ok &= check(f"Table 3 {name} TVE mean (%)",
                            reg.get(name, {}).get("tve_mean_pct"), mean_pct, 0.02)
        # The aggregate maximum cannot be smaller than any region maximum it
        # contains. This caught a stale figure in an earlier draft.
        region_max = [r.get("tve_max_pct") for r in reg.values()
                      if r.get("tve_max_pct") is not None]
        if region_max and ov.get("tve_max_pct") is not None:
            if ov["tve_max_pct"] + 1e-9 < max(region_max):
                print(f"{RED}  [FAIL]{RST} aggregate TVE max {ov['tve_max_pct']:.2f} %% is "
                      f"below the largest region max {max(region_max):.2f} %% - impossible")
                all_ok = False
            else:
                print(f"{GREEN}  [ OK ]{RST} aggregate TVE max is consistent with its regions")

    # --- Step 4: TOST sweep (paper Figure 12) ---
    print("\n[4] TOST sensitivity sweep (paper Section 3.2, Figure 12)")
    all_ok &= run([py, script("tost_sweep.py"), args.canonical,
                   "--json", os.path.splitext(args.canonical)[0] + "_tost_sweep.json",
                   "--verify-hash", CANONICAL_SHA256], "tost_sweep")
    sw, _ok = need_json(os.path.splitext(args.canonical)[0] + "_tost_sweep.json",
                        "TOST sweep"); all_ok &= _ok
    if sw:
        wins = sw.get("windows", {})
        for name, expected_delta in (("Pre-event nominal", 30.0),
                                     ("Post-event settled", 20.0),
                                     ("Full 360 s replay", 20.0)):
            all_ok &= check(f"Sec 3.2 {name} equivalence crossover (mHz)",
                            wins.get(name, {}).get("crossover_delta_mHz"),
                            expected_delta, 0.01)
        # The sweep must be monotonic in delta: once equivalence is accepted,
        # widening the margin cannot reject it again. A non-monotonic sweep
        # would indicate a numerical fault rather than a physical result.
        # NOTE: this was a for/else, whose else clause runs whenever the loop
        # finishes without a break - i.e. always. It printed "all sweeps
        # monotonic" immediately after printing a FAIL. Fixed 22 September
        # 2026 to track the result explicitly.
        monotonic = True
        for name, w in wins.items():
            seq = [row.get("equivalent") for row in w.get("sweep", [])]
            if True in seq and False in seq[seq.index(True):]:
                print(f"{RED}  [FAIL]{RST} {name}: sweep is non-monotonic in delta")
                all_ok = False
                monotonic = False
        if monotonic:
            print(f"{GREEN}  [ OK ]{RST} all sweeps monotonic in delta")

    # --- Step 5: replay metrics (paper Figure 11) ---
    print("\n[5] Replay metrics (paper Section 3.1, Figure 11)")
    all_ok &= run([py, script("replay_metrics.py"), args.canonical, "--neso", args.neso],
                  "replay_metrics")
    rm, _ok = need_json(os.path.splitext(args.canonical)[0] + "_metrics.json",
                        "replay metrics"); all_ok &= _ok
    if rm:
        all_ok &= check("Sec 3.1 max RoCoF err (Hz/s)",
                        rm.get("metric_3_max_rocof_tracking_error_hz_per_s"),
                        0.721, 0.005)

    # --- Step 6: runtime (paper Section 3.6) ---
    print("\n[6] Runtime characterisation (paper Section 3.6)")
    if os.path.isfile(args.timing):
        all_ok &= run([py, script("runtime_metrics.py"), args.timing], "runtime_metrics")
    else:
        print(f"{RED}  [SKIP]{RST} {args.timing} not present - Section 3.6 runtime "
              f"figures are NOT verified by this run")
        all_ok = False

    # --- Step 7: hybrid V3 (paper Table 4, V3 row) ---
    print("\n[7] Hybrid V1/V2 -> V3 (paper Section 3.5, Table 4, V3 row)")
    if os.path.isfile(args.v1):
        all_ok &= run([py, script("hybrid_freq_analysis.py"), args.v1, args.canonical], "hybrid")
        hy, _ok = need_json(data("hybrid_output.json"), "hybrid V3"); all_ok &= _ok
        if hy:
            all_ok &= check("V3 hybrid std (mHz)",
                            hy.get("hybrid_stats", {}).get("std_mHz"),
                            116.75, 1.0)
    else:
        print(f"{DIM}  (skipped: {args.v1} not present){RST}")

    # --- Step 8: Kalman V4 (paper Table 4, V4 row) ---
    print("\n[8] Kalman V4 (paper Section 3.5, Table 4, V4 row) - TRUE output, first-principles tuning")
    all_ok &= run([py, script("kalman_freq.py"), args.canonical, "--neso", args.neso] + figflag,
                  "kalman_freq")
    kf, _ok = need_json(data("kalman_output.json"), "Kalman V4"); all_ok &= _ok
    v4_std = None
    if kf:
        v4_std = kf.get("v4_stats", {}).get("std_mhz")
        all_ok &= check("Table 4 V4 Kalman std (mHz)", v4_std, 194.25, 1.0)
        print(f"  {DIM}V4 true std = {v4_std:.2f} mHz "
              f"(improvement {kf.get('v4_improvement_over_v2_std_pct'):+.1f}% "
              f"over V2; matches Table 4){RST}")

    # --- Step 9: swing-equation inversion on external inputs (paper Figure 18) ---
    # f0 and dP are published; the 0.16 Hz/s initial RoCoF is a nominal figure
    # (Section 3.11). None of the three is measured by the instrument.
    print("\n[9] Swing-equation inversion on external inputs (paper Section 3.11, Figure 18)")
    all_ok &= run([py, script("inertia_estimator.py"),
                   "--measured-rocof", args.canonical] + figflag, "inertia_estimator")
    ie, _ok = need_json(data("inertia_output.json"), "inertia estimator"); all_ok &= _ok
    if ie:
        # 231.4 GVA*s at the published 1,481 MW cumulative loss. The previous
        # expected value, 215.3, was obtained at 1378 MW, which does not appear
        # in the NESO technical report. See inertia_estimator.py.
        all_ok &= check("Inversion on external inputs (GVA*s)",
                        ie.get("primary_recovered_ek_gvas"), 231.4, 0.5)
        prov = ie.get("rocof_provenance")
        if prov:
            print(f"  {DIM}scope: {prov}{RST}")

    # --- Step 10: CRB (paper Section 3.5), chained to true V4 ---
    print("\n[10] Cramer-Rao bound (paper Section 3.5, Table 4 ratios), using TRUE V4 std")
    crb_cmd = [py, script("crb_analysis.py")]
    if v4_std is not None:
        crb_cmd += ["--v4", f"{v4_std:.4f}"]
    all_ok &= run(crb_cmd, "crb_analysis")
    cr, _ok = need_json(data("crb_output.json"), "Cramer-Rao"); all_ok &= _ok
    if cr:
        # CORRECTED 23 September 2026 together with crb_analysis.py and
        # Table 4 of the manuscript, in one change, as the note below required.
        # The superseded value was 2.757 mHz. crb_analysis.py evaluated the
        # Rife-Boorstyn expression with the coefficient 6, which is the
        # complex-exponential form; a real sampled sinusoid takes 12, and the
        # bound is therefore larger by exactly sqrt(2): 3.899 mHz. Confirmed
        # 22 September 2026 by closed form and by inverting the 3x3 Fisher
        # information matrix for (A, f, phi), which gives 3.901 mHz. The
        # whole-replay ratios (quoted in Table 4 of an earlier draft; the
        # published Table 4 quotes stationary-block ratios computed by
        # paper_tables.py) move from 47.5 / 75.8 / 42.4 / 70.5x to
        # 33.6 / 53.6 / 29.9 / 49.8x.
        # SAMPLING RATE CORRECTED 23 Sep 2026: the canonical M2
        # build's per-channel rate is a MEASURED 1816.5 Hz, not the assumed
        # 2000 Hz. The bound is linear in f_s at fixed N, so it falls to
        # 3.541 mHz and the ratios rise to 37.0 / 59.0 / 33.0 / 54.9x.
        # All three were changed together on 23 September 2026.
        all_ok &= check("CRB std conservative (mHz)",
                        cr.get("crb_std_conservative_mhz"), 3.541, 0.01)

    # --- Step 11: derived-output hashes quoted in the Data Availability Statement ---
    print("\n[11] Derived-output hashes (Data Availability Statement)")
    for fname, expected in DERIVED_HASHES.items():
        fpath = data(fname)
        if not os.path.isfile(fpath):
            print(f"{DIM}  (skipped: {fname} not present){RST}")
            continue
        got = sha256(fpath)
        if got.lower() == expected.lower():
            print(f"{GREEN}  [ OK ]{RST} {fname} ({got[:16]}...)")
        else:
            print(f"{RED}  [FAIL]{RST} {fname}\n        got      {got}\n"
                  f"        expected {expected}")
            all_ok = False

    # --- Step 12: per-region tables and statistics (added in release v1.2.0) ---
    # Tables 2-8 and the per-region statistics of Sections 2.3, 3.1, 3.2 and 3.8
    # were computed by no released script before v1.2.0. paper_tables.py
    # regenerates every one of them from the released captures and checks each
    # against the manuscript; it exits non-zero if any value moves.
    print("\n[12] Per-region tables and statistics (paper_tables.py --check; Tables 2-8)")
    stale_pt = data("paper_tables_output.json")
    if os.path.isfile(stale_pt):
        os.remove(stale_pt)
    all_ok &= run([py, script("paper_tables.py"), "--check"], "paper_tables")
    pt, _ok = need_json(data("paper_tables_output.json"), "paper tables"); all_ok &= _ok
    if pt:
        all_ok &= check("Sec 2.3 realised noise floor (mHz)",
                        pt.get("canonical", {}).get("noise_floor_pre_event_mhz"), 5.46, 0.005)
        all_ok &= check("Table 8 overall mean |e| (mHz)",
                        pt.get("table8", {}).get("regions", {}).get("Overall", {}).get("mean_abs_mhz"),
                        78.3, 0.05)

    # --- Step 13: figures (added in release v1.2.0) ---
    if args.skip_figures:
        print(f"\n{DIM}[13] Figures skipped (--skip-figures){RST}")
    else:
        print("\n[13] Figures 11-15, 17 and 18 (paper_figures.py)")
        # Written to figures_regenerated/, not figures/: a rendered PNG is not byte-
        # reproducible across matplotlib versions, and overwriting the released,
        # manifested copies would make `make_hashes.py --check` fail for a reader who
        # had done nothing wrong. Compare the two folders by eye.
        all_ok &= run([py, script("paper_figures.py"), "--outdir",
                       os.path.join(_ROOT, "figures_regenerated")], "paper_figures")

    print("\n" + "=" * 74)
    if all_ok:
        print(f"{GREEN}  ALL CHECKS PASSED{RST}")
    else:
        print(f"{RED}  SOME CHECKS FAILED - see above{RST}")
    print("=" * 74)
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())