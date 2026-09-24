"""
run_repeatability_check.py
==========================

Aggregate N independent V2 captures of the SAME NESO trajectory (default:
9 August 2019) and report the run-to-run spread of every headline metric.
This supports a repeatability claim of the form "independent captures
agree to within run-to-run tolerance" for the canonical V2 build. The
repeatability result reported in the paper (Section 3.8) is for the V6
parallel build and is reproduced by analyse_v6_repeatability.py.

Each capture is passed through run_live_verification.py individually,
then this wrapper aggregates the per-run summary JSONs and reports:

  metric               mean        stdev       min         max         n
  --------------------------------------------------------------------
  BA aggregate RMSE    209.02      1.24        207.30      210.55      3
  ...

USAGE

    python run_repeatability_check.py \
        --v2 replay_20260722_run1.csv replay_20260722_run2.csv replay_20260722_run3.csv \
        --v1 replay_20260722_v1.csv \
        --neso neso_trajectory.csv \
        --outroot repeatability

Produces:
    repeatability/run_1/                 (individual run summary)
    repeatability/run_2/
    repeatability/run_3/
    repeatability/repeatability_summary.json      (the aggregate)
    repeatability/repeatability_table.txt         (human-readable)

Author: Jack Davies
"""

from __future__ import annotations
import argparse, hashlib, json, math, os, statistics, subprocess, sys

GREEN = "\033[92m"; DIM = "\033[2m"; RST = "\033[0m"; YEL = "\033[93m"; RED = "\033[91m"


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for c in iter(lambda: fh.read(65536), b''):
            h.update(c)
    return h.hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--v2', nargs='+', required=True,
        help='N independent V2 captures of the same trajectory')
    ap.add_argument('--v1', required=True,
        help='Single V1 baseline (for the V3 hybrid step in each run)')
    ap.add_argument('--neso', default='neso_trajectory.csv')
    ap.add_argument('--outroot', default='repeatability')
    ap.add_argument('--python', default=sys.executable)
    args = ap.parse_args(argv)

    if len(args.v2) < 2:
        print(f"{RED}Need at least 2 captures for a spread report.{RST}", file=sys.stderr)
        return 2

    os.makedirs(args.outroot, exist_ok=True)
    per_run_summaries = []

    print("=" * 92)
    print(f"  REPEATABILITY CHECK — {len(args.v2)} independent V2 captures")
    print("=" * 92)
    for i, v2 in enumerate(args.v2, 1):
        run_dir = os.path.join(args.outroot, f'run_{i}')
        os.makedirs(run_dir, exist_ok=True)
        print(f"\n  --- Run {i}: {os.path.basename(v2)} (sha256 {sha256(v2)[:12]}...) ---")
        cmd = [args.python, 'run_live_verification.py',
               '--v2', v2, '--v1', args.v1,
               '--neso', args.neso, '--outdir', run_dir]
        r = subprocess.run(cmd)
        if r.returncode not in (0, 1):
            print(f"{YEL}  run {i} returned {r.returncode}; continuing{RST}")
        summary_path = os.path.join(run_dir, 'live_verification_summary.json')
        if os.path.isfile(summary_path):
            with open(summary_path, encoding='utf-8') as fh:
                per_run_summaries.append(json.load(fh))

    if not per_run_summaries:
        print(f"{RED}No per-run summaries produced.{RST}")
        return 2

    # Aggregate per-metric across runs
    metrics_by_name = {}
    for s in per_run_summaries:
        for row in s['comparison']:
            m = row['metric']
            v = row['live']
            if v is None:
                continue
            metrics_by_name.setdefault(m, []).append(v)

    print("\n" + "=" * 96)
    print(f"  RUN-TO-RUN SPREAD (n = {len(per_run_summaries)})")
    print("=" * 96)
    header = f"  {'metric':<28} {'mean':>10} {'stdev':>10} {'min':>10} {'max':>10} {'spread/mean %':>14}"
    print(header)
    print("  " + "-" * (len(header) - 2))
    rows_out = []
    for m, vals in metrics_by_name.items():
        if len(vals) < 2:
            continue
        mean = statistics.mean(vals)
        sd = statistics.stdev(vals)
        spread_pct = (max(vals) - min(vals)) / abs(mean) * 100 if mean != 0 else float('nan')
        print(f"  {m:<28} {mean:>10.3f} {sd:>10.3f} {min(vals):>10.3f} {max(vals):>10.3f} "
              f"{spread_pct:>14.2f}")
        rows_out.append({
            'metric': m, 'n': len(vals), 'mean': mean, 'stdev': sd,
            'min': min(vals), 'max': max(vals), 'spread_pct_of_mean': spread_pct,
        })
    print("=" * 96)

    out = {
        'schema_version': '1.0',
        'tool': 'run_repeatability_check.py',
        'n_runs': len(per_run_summaries),
        'v2_captures': [os.path.basename(v) for v in args.v2],
        'v2_sha256':   [sha256(v) for v in args.v2],
        'v1_capture':  os.path.basename(args.v1),
        'v1_sha256':   sha256(args.v1),
        'spread': rows_out,
    }
    with open(os.path.join(args.outroot, 'repeatability_summary.json'), 'w') as fh:
        json.dump(out, fh, indent=2)
    print(f"  Wrote {args.outroot}/repeatability_summary.json")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
