"""
bland_altman.py
===============

Computes Bland-Altman agreement statistics and 95 % confidence intervals for
the V2 M2 frequency-detector output against the NESO reference trajectory
from the canonical replay CSV, and renders Figure 10 of the paper.

Bland-Altman analysis is a standard measurement-agreement methodology
(Bland & Altman, 1986) that decomposes the disagreement between two
measurement techniques into (i) the mean bias between them, and (ii) the
95 % limits of agreement (LoA), defined as bias ± 1.96 × SD_difference.

Usage:
    python bland_altman.py replay_20260627_102507V2.csv
    python bland_altman.py replay_20260627_102507V2.csv --json ba_out.json --figure Figure_10.png

Standard library plus matplotlib only.

Reproducibility:
    Expected input  : replay_20260627_102507V2.csv (canonical V2 capture)
                      SHA-256 dadbcfe247dd4e538b71a891b2b26f070d02c42e3dd4622b5d2125aff169560f
    Expected outputs: <csv>_bland_altman.json  and  Figure_10_Bland_Altman.png

Author: Jack Davies
"""

from __future__ import annotations
import argparse, csv, json, math, os, statistics, sys


REGIONS = [
    ('pre-event nominal', 0, 45),
    ('descent',           46, 119),
    ('nadir',             120, 200),
    ('late-recovery',     201, 300),
    ('settled',           301, 359),
]

REGION_COLOURS = {
    'pre-event nominal': '#2b8cbe',
    'descent':           '#e34a33',
    'nadir':             '#a50026',
    'late-recovery':     '#fd8d3c',
    'settled':           '#2ca02c',
}


def load_canonical(path):
    frames = []
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        r = csv.DictReader(fh)
        for row in r:
            try:
                frames.append({
                    'replay_idx': int(row['replay_index']),
                    'f_target':   float(row['f_target_hz']),
                    'f_measured': float(row['f_Vb']),
                })
            except (ValueError, KeyError):
                continue
    return frames


def region_for(idx):
    for name, lo, hi in REGIONS:
        if lo <= idx <= hi:
            return name
    return None


def pearson_r(x, y):
    n = len(x)
    if n < 2:
        return None
    mx = statistics.mean(x)
    my = statistics.mean(y)
    num = sum((xi-mx)*(yi-my) for xi, yi in zip(x, y))
    dx = math.sqrt(sum((xi-mx)**2 for xi in x))
    dy = math.sqrt(sum((yi-my)**2 for yi in y))
    if dx == 0 or dy == 0:
        return None
    return num / (dx * dy)


def stats_for_errors(errors_mHz, ref, meas):
    n = len(errors_mHz)
    if n < 2:
        return None
    mean = statistics.mean(errors_mHz)
    std = statistics.stdev(errors_mHz)
    se = std / math.sqrt(n)
    return {
        'n':               n,
        'mean_mHz':        mean,
        'std_mHz':         std,
        'ci_95_lo_mHz':    mean - 1.96 * se,
        'ci_95_hi_mHz':    mean + 1.96 * se,
        'rmse_mHz':        math.sqrt(sum(e*e for e in errors_mHz) / n),
        'max_abs_mHz':     max(abs(e) for e in errors_mHz),
        'pearson_r':       pearson_r(ref, meas),
        # Bland-Altman
        'ba_bias_mHz':     mean,
        'ba_sd_diff_mHz':  std,
        'ba_loa_upper_mHz': mean + 1.96 * std,
        'ba_loa_lower_mHz': mean - 1.96 * std,
    }


def render_figure(frames, aggregate_stats, out_path):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    except ImportError:
        print("Figure rendering requires matplotlib. Install with: pip install matplotlib",
              file=sys.stderr)
        return False

    diffs_mHz = [(f['f_measured'] - f['f_target']) * 1000 for f in frames]
    means_Hz = [(f['f_measured'] + f['f_target']) / 2.0 for f in frames]
    regions = [region_for(f['replay_idx']) for f in frames]

    bias = aggregate_stats['ba_bias_mHz']
    loa_upper = aggregate_stats['ba_loa_upper_mHz']
    loa_lower = aggregate_stats['ba_loa_lower_mHz']

    fig, ax = plt.subplots(figsize=(9, 6))
    for region_name in [r[0] for r in REGIONS]:
        xs = [means_Hz[i] for i in range(len(diffs_mHz)) if regions[i] == region_name]
        ys = [diffs_mHz[i] for i in range(len(diffs_mHz)) if regions[i] == region_name]
        ax.scatter(xs, ys, c=REGION_COLOURS[region_name], s=18, alpha=0.75,
                   edgecolors='none', label=region_name)

    ax.axhline(bias, color='#333333', linestyle='-', linewidth=1.2,
               label=f'bias = {bias:+.1f} mHz')
    ax.axhline(loa_upper, color='#333333', linestyle='--', linewidth=1.0,
               label=f'+1.96 SD = {loa_upper:+.1f} mHz')
    ax.axhline(loa_lower, color='#333333', linestyle='--', linewidth=1.0,
               label=f'-1.96 SD = {loa_lower:+.1f} mHz')
    ax.axhspan(-100, 100, alpha=0.08, color='green', zorder=0)
    ax.text(min(means_Hz), 90, 'IEC 61400-21 Class II band (±100 mHz)',
            fontsize=8, ha='left', va='top', color='#2d6b2d', style='italic')

    ax.set_xlabel('(measured + reference) / 2  [Hz]', fontsize=11)
    ax.set_ylabel('measured - reference  [mHz]', fontsize=11)
    ax.set_title('Bland-Altman agreement plot: V2 M2 measurement vs NESO reference',
                 fontsize=11)
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.legend(loc='lower right', fontsize=8, framealpha=0.9)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close()
    return True


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('csv_path')
    ap.add_argument('--json', default=None)
    ap.add_argument('--figure', default=None)
    args = ap.parse_args(argv)

    if not os.path.isfile(args.csv_path):
        print(f'ERROR: input not found: {args.csv_path}', file=sys.stderr)
        return 2

    frames = load_canonical(args.csv_path)
    if not frames:
        print(f'ERROR: no usable frames', file=sys.stderr)
        return 3

    # Per-region + aggregate
    results = {}
    for region_name, lo, hi in REGIONS:
        region_frames = [f for f in frames if lo <= f['replay_idx'] <= hi]
        if not region_frames:
            continue
        errors = [(f['f_measured'] - f['f_target']) * 1000 for f in region_frames]
        ref = [f['f_target'] for f in region_frames]
        meas = [f['f_measured'] for f in region_frames]
        results[region_name] = stats_for_errors(errors, ref, meas)

    agg_errors = [(f['f_measured'] - f['f_target']) * 1000 for f in frames]
    agg_ref = [f['f_target'] for f in frames]
    agg_meas = [f['f_measured'] for f in frames]
    results['aggregate'] = stats_for_errors(agg_errors, agg_ref, agg_meas)

    # Console report
    print('=' * 100)
    print(' BLAND-ALTMAN AGREEMENT STATISTICS: V2 M2 vs NESO Reference')
    print('=' * 100)
    print(f'{"Region":<20} {"n":>4} {"mean":>8} {"CI_lo":>8} {"CI_hi":>8} '
          f'{"std":>8} {"RMSE":>8} {"max|e|":>8} {"LoA-":>8} {"LoA+":>8} {"r":>7}')
    print('-' * 100)
    for name in [r[0] for r in REGIONS] + ['aggregate']:
        s = results.get(name)
        if s is None:
            continue
        r = s['pearson_r'] if s['pearson_r'] is not None else float('nan')
        print(f'{name:<20} {s["n"]:>4d} {s["mean_mHz"]:>8.2f} '
              f'{s["ci_95_lo_mHz"]:>8.2f} {s["ci_95_hi_mHz"]:>8.2f} '
              f'{s["std_mHz"]:>8.2f} {s["rmse_mHz"]:>8.2f} '
              f'{s["max_abs_mHz"]:>8.2f} {s["ba_loa_lower_mHz"]:>8.1f} '
              f'{s["ba_loa_upper_mHz"]:>8.1f} {r:>7.4f}')
    print('=' * 100)

    # JSON output
    out_json = args.json or (os.path.splitext(args.csv_path)[0] + '_bland_altman.json')
    with open(out_json, 'w', encoding='utf-8') as fh:
        json.dump({
            'schema_version': '1.0',
            'tool': 'bland_altman.py',
            'source_file': os.path.basename(args.csv_path),
            'per_region': {k: v for k, v in results.items() if k != 'aggregate'},
            'aggregate': results['aggregate'],
        }, fh, indent=2)
    print(f'JSON sidecar written: {out_json}')

    # Figure
    fig_path = args.figure or 'Figure_10_Bland_Altman.png'
    if render_figure(frames, results['aggregate'], fig_path):
        print(f'Figure written: {fig_path}')

    return 0


if __name__ == '__main__':
    sys.exit(main())
