"""
task_1_4_metrics.py
===================

Computes the σ (standard deviation) of the M1, M3, M4, M5 module outputs
across a sliced window of a replay CSV: the four module-precision figures
of a table in an earlier draft of the Open RTU paper (the published
manuscript carries no module-precision table).

Default slice: replay indices 0-45 (pre-event nominal region of the V2
canonical capture). The signal in this region is the NESO 50.045 Hz
baseline +/- < 100 mHz - operationally equivalent to nominal 50 Hz for
the purposes of M1/M3/M4/M5 sigma characterisation. Cross-validation
against the post-event settled region (i=301-359) is recommended.

Usage:
    python task_1_4_metrics.py replay_20260627_102507.csv
    python task_1_4_metrics.py replay_20260627_102507.csv --i-start 0 --i-end 45
    python task_1_4_metrics.py replay_20260627_102507.csv --i-start 301 --i-end 359

Writes a JSON sidecar <csv>_task_1_4.json alongside the input CSV.
"""

from __future__ import annotations
import argparse, csv, json, os, statistics


def load_rows(path, i_start, i_end):
    rows = []
    with open(path, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            try:
                i = int(r['replay_index'])
                if i < i_start or i > i_end:
                    continue
                rows.append(r)
            except (KeyError, ValueError):
                continue
    return rows


def col_stats(rows, key):
    vals = []
    for r in rows:
        try:
            vals.append(float(r[key]))
        except (KeyError, ValueError):
            continue
    if len(vals) < 2:
        return None
    return {
        'n':     len(vals),
        'mean':  statistics.mean(vals),
        'stdev': statistics.stdev(vals),
        'min':   min(vals),
        'max':   max(vals),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('csv', help="Replay CSV from SCADA_HMI.py")
    ap.add_argument('--i-start', type=int, default=0,
                    help='First replay index inclusive (default 0)')
    ap.add_argument('--i-end',   type=int, default=45,
                    help='Last  replay index inclusive (default 45)')
    args = ap.parse_args()

    rows = load_rows(args.csv, args.i_start, args.i_end)
    if not rows:
        raise SystemExit(f"No rows in window i={args.i_start}..{args.i_end}")

    # Column mapping. The canonical channel for each module:
    #   M1 -> Vb_rms (clean DAC2 voltage channel)
    #   M3 -> Ph_VbIb (phase pair stabilised by GPIO27 sync)
    #   M4 -> Pb (per-phase active power on phase b)
    #   M5 -> THD_Vb (THD on canonical voltage channel)
    targets = [
        ('M1 V_rms (Vb_rms, V)',     'Vb_rms'),
        ('M3 Phase Vb/Ib (deg)',     'Ph_VbIb'),
        ('M4 Active power P (Pb, W)','Pb'),
        ('M5 THD (THD_Vb, %)',       'THD_Vb'),
    ]

    print()
    print(f"=== Module sigma statistics ===")
    print(f"  CSV: {os.path.basename(args.csv)}")
    print(f"  Window: replay_index = {args.i_start}..{args.i_end}  (n = {len(rows)})")
    print()
    print(f"  {'Module':<28}{'n':>5}{'mean':>13}{'sigma':>12}{'min':>11}{'max':>11}")
    print(f"  {'-'*80}")

    out = {'csv': os.path.basename(args.csv),
           'window': [args.i_start, args.i_end],
           'modules': {}}
    for label, col in targets:
        s = col_stats(rows, col)
        if s is None:
            print(f"  {label:<28}  -- column '{col}' not found or insufficient data --")
            continue
        print(f"  {label:<28}{s['n']:>5d}{s['mean']:>13.4f}{s['stdev']:>12.4f}"
              f"{s['min']:>11.4f}{s['max']:>11.4f}")
        out['modules'][label] = s

    out_path = os.path.splitext(args.csv)[0] + '_task_1_4.json'
    with open(out_path, 'w', encoding='utf-8') as fh:
        json.dump(out, fh, indent=2)
    print()
    print(f"  -> {out_path}")


if __name__ == '__main__':
    main()
