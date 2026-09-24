"""
replay_metrics.py
=================

Standalone post-processor for a captured replay CSV. Computes four
replay metrics (the maximum RoCoF tracking error of Section 3.1 of the paper
is metric 3):

    1. Max |f_instrument(t) - f_NESO(t)| across the full window     (mHz)
       Target: +- 10 mHz  (IEC 61400-21 Class I band)
    2. FREQ_alarm timestamp accuracy vs the 49.5 Hz NESO crossing    (ms)
       Target: +- 20 ms
    3. M2 RoCoF tracking error vs first-difference of NESO            (Hz/s)
       Target: 0.05 Hz/s resolution
    4. DNP3 frame count and CRC failure count over the window        (count)
       Target: zero CRC failures
       (CRC failures = gaps in frame_count; the M6 sketch increments
        frame_count only on successful frame emission.)

This script is intentionally separate from SCADA_HMI.py so that you can
re-run it as many times as you like against the same CSV without
re-running the bench session.

Usage:
    python replay_metrics.py replay_20260625_143052.csv
    python replay_metrics.py replay_*.csv --neso neso_trajectory.csv

Outputs a one-page text summary on stdout plus a JSON sidecar
'<csv>_metrics.json' alongside each input CSV.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
from typing import Optional


# ─────────────────────────────── Loaders ──────────────────────────────────
def load_neso(path: str) -> dict[int, float]:
    """Returns {replay_index: frequency_hz}."""
    out = {}
    with open(path, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            out[int(row['index'])] = float(row['frequency_hz'])
    return out


def load_capture(path: str) -> list[dict]:
    """Returns list of {i, f_Vb, freq_alarm, frame_count, iso_time, monotonic_s}
    for replay-tagged rows only (replay_index >= 0). Sorted by index."""
    rows = []
    with open(path, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            try:
                i = int(r['replay_index'])
                if i < 0:
                    continue
                rows.append({
                    'i':           i,
                    'f_Vb':        float(r['f_Vb']),
                    'freq_alarm':  int(r['FREQ_alarm']),
                    'frame_count': int(r['frame_count']),
                    'monotonic_s': float(r['monotonic_s']),
                    'iso_time':    r['iso_time'],
                })
            except (KeyError, ValueError):
                continue
    rows.sort(key=lambda r: r['i'])
    return rows


# ─────────────────────────────── Metrics ──────────────────────────────────
def metric1_max_abs_error(rows, neso) -> tuple[float, int]:
    """Returns (max_abs_error_mhz, worst_index)."""
    worst = 0.0
    worst_i = -1
    for r in rows:
        ref = neso.get(r['i'])
        if ref is None:
            continue
        e = abs(r['f_Vb'] - ref) * 1000.0
        if e > worst:
            worst, worst_i = e, r['i']
    return worst, worst_i


def metric2_alarm_timing(rows, neso, threshold_hz=49.5) -> Optional[float]:
    """Returns timestamp delta in ms between:
      - the first row where FREQ_alarm flips 0 -> 1
      - the first NESO index where ref drops below threshold_hz
    Negative = HMI alarm fired BEFORE NESO crossing (early/false trip).
    Positive = HMI alarm fired AFTER  NESO crossing (lag).
    Returns None if either event is absent.
    """
    # First NESO index where the reference goes sub-49.5
    neso_cross_i = None
    for i in sorted(neso):
        if neso[i] < threshold_hz:
            neso_cross_i = i
            break
    if neso_cross_i is None:
        return None

    # First HMI row where FREQ_alarm = 1 (and the previous frame was 0)
    alarm_row = None
    prev_alarm = 0
    for r in rows:
        if prev_alarm == 0 and r['freq_alarm'] == 1:
            alarm_row = r
            break
        prev_alarm = r['freq_alarm']
    if alarm_row is None:
        return None

    # NESO indices step at 1 s; HMI timestamps are in monotonic_s from
    # the start of the logger. Anchor to replay_index 0 = monotonic_s 0
    # by linear interpolation.
    base = next((r for r in rows if r['i'] == 0), rows[0])
    hmi_time_at_alarm = alarm_row['monotonic_s'] - base['monotonic_s']
    neso_time_at_cross = float(neso_cross_i)  # 1 s/index
    return (hmi_time_at_alarm - neso_time_at_cross) * 1000.0


def metric3_rocof_tracking(rows, neso) -> tuple[float, float]:
    """First-difference RoCoF tracking error.
    Returns (max_abs_error_hz_per_s, rms_error_hz_per_s).
    NESO sample spacing is 1 s, so first-difference == RoCoF in Hz/s.
    """
    # Build index-aligned arrays
    indices = sorted(set(r['i'] for r in rows) & set(neso.keys()))
    meas = {r['i']: r['f_Vb'] for r in rows}
    ref = neso
    errs = []
    for k in range(1, len(indices)):
        i, j = indices[k - 1], indices[k]
        if (j - i) != 1:
            continue  # skip gaps
        rocof_meas = meas[j] - meas[i]
        rocof_ref = ref[j] - ref[i]
        errs.append(rocof_meas - rocof_ref)
    if not errs:
        return 0.0, 0.0
    max_abs = max(abs(e) for e in errs)
    rms = math.sqrt(sum(e * e for e in errs) / len(errs))
    return max_abs, rms


def metric4_frame_count(rows) -> tuple[int, int, int]:
    """Returns (frames_received, expected_min, gaps).
    'gaps' = number of skipped frame_count values (proxy for CRC failures
    or other dropped frames since the M6 sketch increments frame_count
    only on emit).
    """
    if not rows:
        return 0, 0, 0
    counts = [r['frame_count'] for r in rows]
    gaps = 0
    for k in range(1, len(counts)):
        delta = counts[k] - counts[k - 1]
        if delta > 1:
            gaps += (delta - 1)
        elif delta < 0:
            # counter rolled (board reset mid-run) - flag as a gap
            gaps += 1
    # Expected: one frame per second of replay, i.e. (last_index - first_index + 1)
    first_i = rows[0]['i']
    last_i = rows[-1]['i']
    expected = last_i - first_i + 1
    return len(rows), expected, gaps


# ───────────────────────────── Main / CLI ─────────────────────────────────
def process(csv_path: str, neso_path: str) -> dict:
    rows = load_capture(csv_path)
    neso = load_neso(neso_path)
    if not rows:
        raise SystemExit(f"{csv_path}: no replay-tagged rows.")
    if not neso:
        raise SystemExit(f"{neso_path}: no NESO reference rows.")

    max_err_mhz, worst_i = metric1_max_abs_error(rows, neso)
    alarm_delta_ms = metric2_alarm_timing(rows, neso)
    rocof_max, rocof_rms = metric3_rocof_tracking(rows, neso)
    frames, expected, gaps = metric4_frame_count(rows)

    result = {
        'csv': os.path.basename(csv_path),
        'neso_reference': os.path.basename(neso_path),
        'metric_1_max_abs_freq_error_mhz': round(max_err_mhz, 3),
        'metric_1_worst_replay_index': worst_i,
        'metric_2_alarm_timing_error_ms': (round(alarm_delta_ms, 1)
                                            if alarm_delta_ms is not None else None),
        'metric_3_max_rocof_tracking_error_hz_per_s': round(rocof_max, 4),
        'metric_3_rms_rocof_tracking_error_hz_per_s': round(rocof_rms, 4),
        'metric_4_frames_received': frames,
        'metric_4_frames_expected_min': expected,
        'metric_4_frame_gaps': gaps,
    }

    # Print summary
    print()
    print(f"=== {os.path.basename(csv_path)} ===")
    print(f"  rows in window: {frames}  (expected >= {expected}, gaps {gaps})")
    print()
    print(f"  Metric 1 - max |f_inst - f_NESO|:        "
          f"{max_err_mhz:7.2f} mHz   "
          f"(at index {worst_i})            target +-10 mHz")
    if alarm_delta_ms is None:
        print(f"  Metric 2 - FREQ_alarm timing error:      "
              f"     n/a   (alarm did not fire or NESO never < 49.5)")
    else:
        print(f"  Metric 2 - FREQ_alarm timing error:      "
              f"{alarm_delta_ms:+7.1f} ms                              target +-20 ms")
    print(f"  Metric 3 - RoCoF tracking error (max):   "
          f"{rocof_max:7.4f} Hz/s                          target 0.05 Hz/s")
    print(f"  Metric 3 -                       (RMS):  "
          f"{rocof_rms:7.4f} Hz/s")
    print(f"  Metric 4 - DNP3 frame integrity:         "
          f"{frames} received, {gaps} gap(s)                target 0 gaps")
    print()

    # Write JSON sidecar
    json_path = os.path.splitext(csv_path)[0] + "_metrics.json"
    with open(json_path, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(result, fh, indent=2)
    print(f"  -> {json_path}")
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('csv', nargs='+',
                    help="Replay CSV(s) produced by SCADA_HMI.py. "
                         "Globs are expanded.")
    ap.add_argument('--neso', default='neso_trajectory.csv',
                    help="Path to the NESO reference CSV "
                         "(default: ./neso_trajectory.csv)")
    args = ap.parse_args()

    paths = []
    for pat in args.csv:
        matched = sorted(glob.glob(pat))
        paths.extend(matched if matched else [pat])

    for p in paths:
        if not os.path.isfile(p):
            print(f"WARN: not found: {p}")
            continue
        try:
            process(p, args.neso)
        except SystemExit as e:
            print(f"ERROR processing {p}: {e}")


if __name__ == '__main__':
    main()
