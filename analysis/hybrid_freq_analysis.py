"""
hybrid_freq_analysis.py

Off-line Python characterisation of the hybrid V1/V2 frequency-detector
selection strategy, applied to the existing V1 baseline
(replay_20260625_180520.csv) and V2 canonical (replay_20260627_102507V2.csv)
captures. The script emulates the state machine of hybrid_freq_detector.cpp
in software, using each captured frame's fVb estimate from V1 and V2 as
the two candidate estimates. It produces a hybrid_output.csv, a
hybrid_output.json sidecar, and a Figure that mirrors the paper's Bland-
Altman analysis but against the hybrid estimate rather than V2 alone.

The purpose of this script is to validate the hybrid design choice before
committing the firmware change to bench validation. If the projected
aggregate hybrid RMSE against the NESO reference is materially better than
either V1 or V2 alone, the firmware implementation is worth pursuing.

Compatible with Python 3.10+ standard library only for the numerical
computation; matplotlib is used for the figure but is optional.

Author: Jack Davies
"""

import csv, math, statistics, json, sys, os, argparse

# Configuration constants — mirror hybrid_freq_detector.h
F_OP_LO_HZ      = 49.5
F_OP_HI_HZ      = 50.5
F_HYST_LO_HZ    = 49.7
F_HYST_HI_HZ    = 50.3
ROCOF_ENTER_HZ  = 0.5
ROCOF_EXIT_HZ   = 0.1
HYST_FRAMES     = 5

MODE_STEADY_V2  = 0
MODE_TRANSIENT_V1 = 1


def load_frames(path):
    """Load a canonical replay CSV into a list of dicts."""
    frames = []
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                frames.append({
                    'monotonic_s':   float(row['monotonic_s']),
                    'frame_count':   int(row['frame_count']),
                    'replay_idx':    int(row['replay_index']),
                    'f_target_hz':   float(row['f_target_hz']),
                    'f_Vb':          float(row['f_Vb']),
                })
            except (ValueError, KeyError):
                continue
    return frames


def align_v1_v2(v1_frames, v2_frames):
    """
    Align V1 and V2 captures by replay index. Both captures are against the
    same 360-sample NESO trajectory, so replay index matches. Where multiple
    frames map to the same replay index, use the first for V1 and the first
    for V2 — this is deterministic and matches how the firmware would
    behave.
    """
    v1_by_idx = {}
    v2_by_idx = {}
    for f in v1_frames:
        v1_by_idx.setdefault(f['replay_idx'], f)
    for f in v2_frames:
        v2_by_idx.setdefault(f['replay_idx'], f)
    common = sorted(set(v1_by_idx.keys()) & set(v2_by_idx.keys()))
    return [(v1_by_idx[i], v2_by_idx[i]) for i in common]


def emulate_state_machine(aligned):
    """
    Run the hybrid state machine over aligned (V1, V2) frame pairs.
    Returns a list of dicts with the hybrid decision on each frame.
    """
    mode = MODE_STEADY_V2
    steady_streak = 0
    last_f_hz = 50.0
    dt_s = 0.1  # one measurement window at 200 samples * 500 us
    transitions = 0
    out = []

    for v1, v2 in aligned:
        # Candidate estimates
        f_v1 = v1['f_Vb']
        f_v2 = v2['f_Vb']
        # Current-mode estimate for RoCoF calculation
        f_current = f_v2 if mode == MODE_STEADY_V2 else f_v1
        rocof = (f_current - last_f_hz) / dt_s
        rocof_abs = abs(rocof)

        # State machine
        if mode == MODE_STEADY_V2:
            band_exit  = (f_current < F_OP_LO_HZ) or (f_current > F_OP_HI_HZ)
            rocof_high = (rocof_abs > ROCOF_ENTER_HZ)
            if band_exit or rocof_high:
                mode = MODE_TRANSIENT_V1
                steady_streak = 0
                transitions += 1
        else:
            in_band  = (F_HYST_LO_HZ <= f_current <= F_HYST_HI_HZ)
            rocof_lo = (rocof_abs < ROCOF_EXIT_HZ)
            if in_band and rocof_lo:
                steady_streak += 1
                if steady_streak >= HYST_FRAMES:
                    mode = MODE_STEADY_V2
                    steady_streak = 0
                    transitions += 1
            else:
                steady_streak = 0

        # Primary output based on updated mode
        f_primary = f_v2 if mode == MODE_STEADY_V2 else f_v1
        last_f_hz = f_primary

        out.append({
            'replay_idx':      v1['replay_idx'],
            'f_target_hz':     v1['f_target_hz'],
            'f_v1':            f_v1,
            'f_v2':            f_v2,
            'f_hybrid':        f_primary,
            'rocof_hz_s':      rocof,
            'mode':            mode,
        })

    return out, transitions


def stats_of_errors(errors_mHz):
    """Standard summary stats of an error list."""
    n = len(errors_mHz)
    if n < 2:
        return {'n': n}
    return {
        'n':         n,
        'mean_mHz':  statistics.mean(errors_mHz),
        'std_mHz':   statistics.stdev(errors_mHz),
        'rmse_mHz':  math.sqrt(sum(e * e for e in errors_mHz) / n),
        'max_mHz':   max(abs(e) for e in errors_mHz),
    }


def pearson_r(xs, ys):
    n = len(xs)
    if n < 2:
        return None
    mx = statistics.mean(xs)
    my = statistics.mean(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if dx * dy == 0:
        return None
    return num / (dx * dy)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('v1_csv')
    ap.add_argument('v2_csv')
    ap.add_argument('--out-csv',  default='hybrid_output.csv')
    ap.add_argument('--out-json', default='hybrid_output.json')
    ap.add_argument('--figure',   default='Figure_hybrid_bland_altman.png')
    args = ap.parse_args(argv)

    v1_frames = load_frames(args.v1_csv)
    v2_frames = load_frames(args.v2_csv)
    print(f'Loaded {len(v1_frames)} V1 frames and {len(v2_frames)} V2 frames')

    aligned = align_v1_v2(v1_frames, v2_frames)
    print(f'Aligned {len(aligned)} paired frames by replay index')

    hybrid, transitions = emulate_state_machine(aligned)
    print(f'State machine executed {transitions} mode transitions')

    # Statistics
    hybrid_errors = [(h['f_hybrid'] - h['f_target_hz']) * 1000 for h in hybrid]
    v1_errors     = [(h['f_v1']     - h['f_target_hz']) * 1000 for h in hybrid]
    v2_errors     = [(h['f_v2']     - h['f_target_hz']) * 1000 for h in hybrid]
    hybrid_stats  = stats_of_errors(hybrid_errors)
    v1_stats      = stats_of_errors(v1_errors)
    v2_stats      = stats_of_errors(v2_errors)

    # Bland-Altman aggregate (hybrid vs reference)
    ba_bias   = hybrid_stats['mean_mHz']
    ba_sd     = hybrid_stats['std_mHz']
    ba_upper  = ba_bias + 1.96 * ba_sd
    ba_lower  = ba_bias - 1.96 * ba_sd
    pearson   = pearson_r([h['f_target_hz'] for h in hybrid],
                          [h['f_hybrid']    for h in hybrid])

    # Report
    print()
    print('=' * 66)
    print('HYBRID vs V1 vs V2 aggregate error against NESO reference')
    print('=' * 66)
    hdr = f'{"":<10s} {"n":>5s} {"mean":>10s} {"std":>10s} {"RMSE":>10s} {"max":>10s}'
    print(hdr)
    for name, s in [('V1', v1_stats), ('V2', v2_stats), ('HYBRID', hybrid_stats)]:
        print(f'{name:<10s} {s["n"]:>5d} {s["mean_mHz"]:>+10.2f} '
              f'{s["std_mHz"]:>10.2f} {s["rmse_mHz"]:>10.2f} {s["max_mHz"]:>10.1f}')
    print()
    print(f'  Bland-Altman bias (hybrid): {ba_bias:+.2f} mHz')
    print(f'  Bland-Altman LoA:           [{ba_lower:+.1f}, {ba_upper:+.1f}] mHz')
    print(f'  Pearson r (hybrid vs ref):  {pearson:.4f}' if pearson else '')

    # Write CSV
    with open(args.out_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(hybrid[0].keys()))
        writer.writeheader()
        writer.writerows(hybrid)
    print(f'\nWrote {args.out_csv}')

    # Write JSON sidecar
    with open(args.out_json, 'w', encoding='utf-8', newline='\n') as f:
        json.dump({
            'tool':          'hybrid_freq_analysis.py',
            'schema':        '1.0',
            'v1_csv':        os.path.basename(args.v1_csv),
            'v2_csv':        os.path.basename(args.v2_csv),
            'n_frames':      len(hybrid),
            'transitions':   transitions,
            'v1_stats':      v1_stats,
            'v2_stats':      v2_stats,
            'hybrid_stats':  hybrid_stats,
            'bland_altman': {
                'bias_mHz':      ba_bias,
                'sd_diff_mHz':   ba_sd,
                'loa_upper_mHz': ba_upper,
                'loa_lower_mHz': ba_lower,
                'pearson_r':     pearson,
            },
            'config': {
                'F_OP_LO_HZ':      F_OP_LO_HZ,
                'F_OP_HI_HZ':      F_OP_HI_HZ,
                'F_HYST_LO_HZ':    F_HYST_LO_HZ,
                'F_HYST_HI_HZ':    F_HYST_HI_HZ,
                'ROCOF_ENTER_HZ':  ROCOF_ENTER_HZ,
                'ROCOF_EXIT_HZ':   ROCOF_EXIT_HZ,
                'HYST_FRAMES':     HYST_FRAMES,
            },
        }, f, indent=2)
    print(f'Wrote {args.out_json}')

    # Optional figure
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 6))
        means = [(h['f_hybrid'] + h['f_target_hz']) / 2 for h in hybrid]
        diffs = hybrid_errors
        modes = [h['mode'] for h in hybrid]
        colors = ['#2b8cbe' if m == MODE_STEADY_V2 else '#e34a33' for m in modes]
        ax.scatter(means, diffs, c=colors, s=18, alpha=0.75, edgecolors='none')
        ax.axhline(ba_bias, color='#333333', linestyle='-',
                   label=f'bias = {ba_bias:+.1f} mHz')
        ax.axhline(ba_upper, color='#333333', linestyle='--',
                   label=f'+1.96 SD = {ba_upper:+.1f} mHz')
        ax.axhline(ba_lower, color='#333333', linestyle='--',
                   label=f'-1.96 SD = {ba_lower:+.1f} mHz')
        ax.axhspan(-100, 100, alpha=0.08, color='green')
        ax.text(min(means), 90, 'IEC 61400-21 Class II band (±100 mHz)',
                fontsize=8, ha='left', va='top', color='#2d6b2d', style='italic')
        # Legend for mode colours
        ax.scatter([], [], c='#2b8cbe', label='steady mode (V2)')
        ax.scatter([], [], c='#e34a33', label='transient mode (V1)')
        ax.set_xlabel('(hybrid + reference) / 2  [Hz]')
        ax.set_ylabel('hybrid - reference  [mHz]')
        ax.set_title('Bland-Altman: Hybrid V1/V2 estimate vs NESO reference')
        ax.legend(loc='lower right', fontsize=8)
        ax.grid(True, alpha=0.3, linestyle=':')
        plt.tight_layout()
        plt.savefig(args.figure, dpi=200, bbox_inches='tight')
        plt.close()
        print(f'Wrote {args.figure}')
    except ImportError:
        pass

    return 0


if __name__ == '__main__':
    sys.exit(main())
