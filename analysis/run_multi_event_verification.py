"""
run_multi_event_verification.py
===============================

Multi-event variant of run_live_verification.py for the Gate 3 follow-on
groundwork. Runs the analysis toolchain over a V2 capture whose input
trajectory is NOT the 9 August 2019 canonical event — i.e. one of the
prepared multi-event trajectories (January2019Excursions, December2023Peak,
etc.). Reports per-event metrics against the event's own NESO reference
without asserting Aug-2019 tolerance bands (those don't apply — a
different event has a different dynamic signature).

WHY THIS IS SEPARATE FROM run_live_verification.py

run_live_verification.py encodes the 9 Aug 2019 canonical result as
tolerance bands (BA RMSE ~209 mHz, max |f error| ~629 mHz, etc.) so a
fresh bench capture can be judged "within run-to-run tolerance" of the
paper's headline numbers. Those bands are event-specific — a benign
2 mHz-dominated day would produce a much smaller RMSE and would
"fail" the tolerance check even though the instrument is working
correctly. This runner drops the tolerance check and just reports.

USAGE

    python run_multi_event_verification.py \
        --v2 replay_20260722_dec2023peak.csv \
        --neso neso_trajectory_December2023Peak.csv \
        --label December2023Peak \
        --outdir multi_event/December2023Peak

The output multi_event/<label>/multi_event_summary.json is
machine-readable and holds every headline metric plus SHA-256s of the
input files, so it feeds straight into the Gate 3 evidence pack.

Requires bland_altman.py, tost_sweep.py, tve_metrics.py, and
replay_metrics.py on PATH — same as run_live_verification.py.

Author: Jack Davies
"""

from __future__ import annotations
import argparse, hashlib, json, os, subprocess, sys


GREEN = "\033[92m"; DIM = "\033[2m"; RST = "\033[0m"; YEL = "\033[93m"; RED = "\033[91m"


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
    with open(path, encoding='utf-8') as fh:
        return json.load(fh)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--v2', required=True, help='V2 capture CSV from this event')
    ap.add_argument('--neso', required=True,
        help='NESO reference CSV for THIS event (from extract_neso_window.py)')
    ap.add_argument('--label', required=True,
        help='Short event label; used in filenames and the summary JSON')
    ap.add_argument('--outdir', default=None,
        help='Output dir; default is multi_event/<label>')
    ap.add_argument('--python', default=sys.executable)
    args = ap.parse_args(argv)

    od = args.outdir or os.path.join('multi_event', args.label)
    os.makedirs(od, exist_ok=True)
    py = args.python
    metrics = {}

    print("=" * 92)
    print(f"  MULTI-EVENT VERIFICATION — {args.label}")
    print("=" * 92)
    print(f"  V2 capture:  {os.path.basename(args.v2)}")
    print(f"    SHA-256:   {sha256(args.v2)}")
    print(f"  NESO ref:    {os.path.basename(args.neso)}")
    print(f"    SHA-256:   {sha256(args.neso)}")
    print()

    # Bland-Altman
    print(f"{DIM}  [1/4] Bland-Altman ...{RST}")
    run_script([py, 'bland_altman.py', args.v2,
                '--json', f'{od}/bland_altman.json',
                '--figure', f'{od}/Figure_BA_{args.label}.png'])
    ba = load_json(f'{od}/bland_altman.json')
    if ba:
        agg = ba['aggregate']
        metrics['BA aggregate RMSE mHz'] = agg['rmse_mHz']
        metrics['BA aggregate bias mHz'] = agg['mean_mHz']
        metrics['BA Pearson r']          = agg['pearson_r']

    # TOST sweep — reports each event's own equivalence crossover
    print(f"{DIM}  [2/4] TOST sensitivity sweep ...{RST}")
    run_script([py, 'tost_sweep.py', args.v2,
                '--plot', f'{od}/Figure_TOST_{args.label}.png',
                '--json', f'{od}/tost_sweep.json'])
    ts = load_json(f'{od}/tost_sweep.json')
    if ts:
        metrics['TOST crossover full-window mHz'] = ts.get('crossover_full_window_mHz')

    # TVE
    print(f"{DIM}  [3/4] Total Vector Error ...{RST}")
    run_script([py, 'tve_metrics.py', args.v2, '--json', f'{od}/tve.json'])
    tve = load_json(f'{od}/tve.json')
    if tve:
        overall = tve.get('overall', {})
        metrics['TVE aggregate mean %'] = overall.get('tve_mean_pct')
        metrics['TVE aggregate max %']  = overall.get('tve_max_pct')

    # Replay metrics
    print(f"{DIM}  [4/4] Replay metrics ...{RST}")
    run_script([py, 'replay_metrics.py', args.v2, '--neso', args.neso])
    rm = load_json(os.path.splitext(args.v2)[0] + '_metrics.json')
    if rm:
        metrics['Max |f error| mHz']     = rm.get('metric_1_max_abs_freq_error_mhz')
        metrics['Max RoCoF error Hz/s']  = rm.get('metric_3_max_rocof_tracking_error_hz_per_s')
        metrics['RMS RoCoF error Hz/s']  = rm.get('metric_3_rms_rocof_tracking_error_hz_per_s')

    # Report
    print()
    print("=" * 92)
    print(f"  {args.label} — per-event metrics (no archived comparison)")
    print("=" * 92)
    for k, v in metrics.items():
        if v is None:
            tag = f'{YEL}[ -- ]{RST}'
            val = 'not produced'
        else:
            tag = f'{GREEN}[ OK ]{RST}'
            val = f'{v:.3f}' if isinstance(v, float) else str(v)
        print(f'  {tag} {k:<40} {val:>12}')
    print('=' * 92)

    summary = {
        'schema_version': '1.0',
        'tool': 'run_multi_event_verification.py',
        'event_label': args.label,
        'v2_capture': os.path.basename(args.v2),
        'v2_sha256': sha256(args.v2),
        'neso_reference': os.path.basename(args.neso),
        'neso_sha256': sha256(args.neso),
        'metrics': metrics,
        'figures': [
            f'{od}/Figure_BA_{args.label}.png',
            f'{od}/Figure_TOST_{args.label}.png',
        ],
    }
    with open(f'{od}/multi_event_summary.json', 'w', encoding='utf-8') as fh:
        json.dump(summary, fh, indent=2)
    print(f'  Summary written: {od}/multi_event_summary.json')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
