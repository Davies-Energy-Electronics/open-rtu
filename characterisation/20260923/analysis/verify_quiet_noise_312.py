#!/usr/bin/env python3
"""
verify_quiet_noise_312.py - the quiet-input and multiplexer figures of
Section 3.12, recomputed from the eight released quiet-input records.

RECORDS (thermal condition read from each header: "# warm 600 s" present =
warmed after the full 600 s settle; absent = taken from cold)
  single channel
    20260913/A1_quiet_single.txt           warmed   (13 Sep bundle)
    20260913/A2_quiet_single_run2.txt      cold     (13 Sep bundle)
    20260913/A3_quiet_single_run3.txt      cold     (13 Sep bundle)
    20260923/captures/A1_quiet_single_run2.txt   warmed
    20260923/captures/A1_quiet_single_run3.txt   warmed
  round robin
    20260923/captures/A2_quiet_roundrobin_run1.txt   warmed
    20260923/captures/A2_quiet_roundrobin_run2.txt   cold
    20260923/captures/A2_quiet_roundrobin_run3.txt   cold
  (20260923/captures/A1_quiet_single.txt, A2_quiet_roundrobin.txt and
   C1_sweep.txt are byte-identical copies of the three 13 Sep records and are
   not counted again.)

METHOD
  SD per record in codes, both raw (ddof=1) and after removing a least-squares
  line (analyse_noise.py's detrend); mV = codes * 3.3/4095.  Group figures are
  the quadrature (rms) mean of the records in the group - the arithmetic mean
  is printed alongside.  Channel-switching term = sqrt(RR^2 - single^2).
  Spread = half the range over the mean.  Effective bits = 12 - log2(sigma /
  sigma_q), sigma_q = LSB/sqrt12.  SNR = 20 log10((A/sqrt2)/sigma), A = 1.0928 V.

Usage:  python3 verify_quiet_noise_312.py
"""
from __future__ import annotations
import math, os
import numpy as np
import common_312 as c

SINGLE = c.SEPT_QUIET + c.SEPT23_SINGLE


def stats(paths):
    rows = []
    for p in paths:
        codes, meta = c.read_codes(p)
        raw, det = c.sd_raw_and_detrended(codes)
        rows.append({"file": os.path.relpath(p, os.path.join(c.HERE, "..", "..")),
                     "warmed": c.warmed(p), "sd_code_raw": raw, "sd_code_detrended": det,
                     "mv_raw": raw * c.LSB_V * 1e3, "mv_detrended": det * c.LSB_V * 1e3})
    return rows


def group(rows, key):
    xs = [r[key] for r in rows]
    return {"rms": c.rms(xs), "mean": float(np.mean(xs)),
            "half_range_over_mean_pct": 50 * (max(xs) - min(xs)) / np.mean(xs),
            "full_range_over_mean_pct": 100 * (max(xs) - min(xs)) / np.mean(xs)}


def main():
    single, rr = stats(SINGLE), stats(c.RR)
    print("=" * 84)
    print(f"  {'record':52s} {'cond':6s} {'raw mV':>8s} {'detr mV':>8s} {'raw code':>9s}")
    print("-" * 84)
    for r in single + rr:
        print(f"  {r['file']:52s} {'warm' if r['warmed'] else 'cold':6s} "
              f"{r['mv_raw']:8.4f} {r['mv_detrended']:8.4f} {r['sd_code_raw']:9.4f}")
    print("  [paper five records: 0.687, 0.694, 0.801, 0.634, 0.657 mV]")

    out = {"records": single + rr}
    sept = single[:3]
    for key in ("mv_raw", "mv_detrended"):
        g = group(sept, key)
        print(f"\n  13 Sep three records ({key}): arithmetic mean {g['mean']:.4f}, rms "
              f"{g['rms']:.4f} mV, spread half-range/mean {g['half_range_over_mean_pct']:.2f} %,"
              f" full range/mean {g['full_range_over_mean_pct']:.2f} %"
              f"   [paper 0.659 mV, 8.2 %]")
        out[f"sept13_{key}"] = g
    s = group(sept, "mv_detrended")["mean"] * 1e-3
    lsb = s / c.LSB_V
    print(f"  0.659 mV = {0.659e-3/c.LSB_V:.3f} LSB [0.82]; ideal sigma_q {c.SIGMA_Q_V*1e3:.4f} mV "
          f"[0.233]; ratio {0.659e-3/c.SIGMA_Q_V:.3f}x [2.8]; eff bits "
          f"{12-math.log2(0.659e-3/c.SIGMA_Q_V):.2f} [10.5]; SNR "
          f"{20*math.log10(c.A_BENCH/math.sqrt(2)/0.659e-3):.2f} dB [61.4]")
    print(f"  crossing-instant error 0.659 mV / 333 V/s = {0.659e-3/333*1e6:.3f} us [1.98]; "
          f"slope 2*pi*50*1.0928 = {2*math.pi*50*c.A_BENCH:.1f} V/s")

    cond = {}
    for name, rows in (("single", single), ("rr", rr)):
        for w in (True, False):
            sel = [r for r in rows if r["warmed"] == w]
            k = f"{name}_{'warm' if w else 'cold'}"
            cond[k] = {"n": len(sel), **{f"{kk}_{key}": v for key in ("mv_raw", "mv_detrended")
                                         for kk, v in group(sel, key).items()}}
            cond[k]["code_rms_raw"] = c.rms([r["sd_code_raw"] for r in sel])
    print()
    for k, g in cond.items():
        print(f"  {k:12s} n={g['n']}  raw rms {g['rms_mv_raw']:.4f} (mean {g['mean_mv_raw']:.4f}) "
              f"detr rms {g['rms_mv_detrended']:.4f} mV; codes rms {g['code_rms_raw']:.4f}; "
              f"spread (half-range/mean) raw {g['half_range_over_mean_pct_mv_raw']:.2f} %, "
              f"detr {g['half_range_over_mean_pct_mv_detrended']:.2f} %")
    print("  [paper: single warm 0.729, cold 0.646; RR warm 0.910, cold 0.736 mV;"
          " spreads 7.7 % warm, 1.8 % cold; 0.905 codes static quiet input]")
    for key in ("mv_raw", "mv_detrended"):
        sw, sc = cond["single_warm"][f"rms_{key}"], cond["single_cold"][f"rms_{key}"]
        rw, rc = cond["rr_warm"][f"rms_{key}"], cond["rr_cold"][f"rms_{key}"]
        print(f"  {key}: mux term warm {math.sqrt(rw**2-sw**2):.4f}, cold "
              f"{math.sqrt(rc**2-sc**2):.4f} mV [0.545, 0.354]; warm/cold single "
              f"{sw/sc:.3f} [1.13], RR {rw/rc:.3f} [1.24]")
    ws = [r["mv_raw"] for r in single if r["warmed"]] + [r["mv_raw"] for r in rr if r["warmed"]]
    print(f"  no-overlap check: RR warm {[round(r['mv_raw'],3) for r in rr if r['warmed']]} > "
          f"single warm max {max(r['mv_raw'] for r in single if r['warmed']):.3f}; "
          f"RR cold min {min(r['mv_raw'] for r in rr if not r['warmed']):.3f} > single cold max "
          f"{max(r['mv_raw'] for r in single if not r['warmed']):.3f}")
    out["conditions"] = cond
    c.save_json("quiet_noise_312.json", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
