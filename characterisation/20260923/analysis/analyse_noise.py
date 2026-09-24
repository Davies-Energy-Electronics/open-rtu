#!/usr/bin/env python3
"""
analyse_noise.py — input-referred noise of the open-RTU conditioning chain.

Reads a capture from characterise.ino taken with the conditioning-network
input SHORTED, and reports the additive noise floor of the whole chain
referred to the converter input.

This is the measurement that replaces the assumed 40 dB SNR in
crb_analysis.py with a measured figure, and that settles the effective-
resolution inference in Section 2.8.6.

Usage:
    python analyse_noise.py capture_A1.txt
    python analyse_noise.py capture_A1.txt --json noise_A1.json
    python analyse_noise.py A1.txt --compare A2.txt     # single vs round-robin

Writes a JSON sidecar so the figures can be chained into crb_analysis.py
rather than retyped.

Standard library plus numpy.
"""
from __future__ import annotations
import argparse, json, math, sys
import numpy as np

LSB_V      = 3.3 / 4095.0          # 0.806 mV, the paper's stated LSB
SIGMA_Q_V  = LSB_V / math.sqrt(12) # 0.233 mV, ideal uniform quantisation
AMPLITUDE  = 1.060                 # V, conditioned amplitude at the ADC (Sec 2.2)
FREQ_WIN   = 200                   # samples per M2 estimate
FS_NOMINAL = 2000.0


def read_capture(path):
    codes, meta = [], {}
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                body = line.lstrip("#").strip()
                if "," in body:
                    k, _, v = body.partition(",")
                    meta[k.strip()] = v.strip()
                continue
            try:
                codes.append(int(line))
            except ValueError:
                pass
    if not codes:
        sys.exit(f"{path}: no ADC codes found")
    return np.asarray(codes, dtype=np.int32), meta


def analyse(codes, meta, label):
    n = len(codes)
    fs = float(meta.get("fs_hz", FS_NOMINAL))

    mean_code = codes.mean()
    # Detrend: slow thermal/bias drift over an 8 s record is not noise.
    t = np.arange(n)
    slope, intercept = np.polyfit(t, codes.astype(float), 1)
    resid = codes - (slope * t + intercept)

    sd_raw_code  = codes.std(ddof=1)
    sd_code      = resid.std(ddof=1)
    sd_v         = sd_code * LSB_V
    drift_code   = slope * n

    # Effective resolution: how many bits the additive noise costs against
    # an ideal uniform quantiser. NOT an ENOB determination - see 2.8.6.
    bits_lost = math.log2(sd_v / SIGMA_Q_V) if sd_v > SIGMA_Q_V else 0.0
    eff_bits  = 12.0 - bits_lost

    snr_db = 10 * math.log10(AMPLITUDE ** 2 / (2 * sd_v ** 2))

    # Cramer-Rao bound at the estimator's own basis, with this measured noise.
    rho = AMPLITUDE ** 2 / (2 * sd_v ** 2)
    crb_hz = math.sqrt(6.0 * fs * fs /
                       ((2 * math.pi) ** 2 * rho * FREQ_WIN * (FREQ_WIN ** 2 - 1)))

    # Histogram flatness: a healthy noise floor exercises many codes smoothly.
    occupied = int(np.unique(codes).size)

    print(f"\n{'='*66}")
    print(f"  NOISE FLOOR — {label}")
    print(f"{'='*66}")
    print(f"  samples                {n}")
    print(f"  sample rate            {fs:.2f} Hz   (record {n/fs:.3f} s)")
    print(f"  mean code              {mean_code:.2f}  ({mean_code*LSB_V:.4f} V)")
    print(f"  codes occupied         {occupied}")
    print(f"  drift over record      {drift_code:+.2f} codes "
          f"({drift_code*LSB_V*1e3:+.2f} mV) — removed before sd")
    print(f"  {'-'*62}")
    print(f"  sd, raw                {sd_raw_code:.4f} codes")
    print(f"  sd, detrended          {sd_code:.4f} codes  = {sd_v*1e3:.4f} mV"
          f"   <-- sigma_v")
    print(f"  ideal quantisation     {SIGMA_Q_V/LSB_V:.4f} codes  "
          f"= {SIGMA_Q_V*1e3:.4f} mV")
    print(f"  excess over ideal      {sd_v/SIGMA_Q_V:.2f}x  "
          f"({bits_lost:.2f} bits)")
    print(f"  {'-'*62}")
    print(f"  inferred eff. resolution   {eff_bits:.2f} bits of a nominal 12")
    print(f"  implied SNR                {snr_db:.1f} dB")
    print(f"  CRB at N={FREQ_WIN}, this SNR   {crb_hz*1e3:.4f} mHz")
    print(f"{'='*66}")
    print("  Reference points:")
    print(f"    paper Sec 2.8.6 infers   0.45 mV  (0.56 LSB, ~11 bits)")
    print(f"    simulation of the V2 estimator against the 5.40 mHz")
    print(f"    realised floor predicts  2.24 mV  (2.78 LSB, ~8.7 bits)")
    print(f"    crb_analysis.py assumes  40 dB conservative / 74 dB ideal")
    print()
    if sd_v * 1e3 < 1.0:
        print("  >> READS LOW. If sigma_v is well under 1 mV the excess noise is")
        print("     NOT additive at the converter input, and the 5.40 mHz floor")
        print("     has another source. Do not proceed to the paper edit; find it.")
    elif 1.5 < sd_v * 1e3 < 3.2:
        print("  >> CONSISTENT with the simulation. Section 2.8.6 can be rewritten")
        print("     around this measured figure.")
    else:
        print("  >> OUTSIDE the predicted band. Record it and investigate before")
        print("     changing anything in the manuscript.")
    print()

    return {
        "label": label, "n": n, "fs_hz": fs,
        "mean_code": float(mean_code), "codes_occupied": occupied,
        "sd_code_detrended": float(sd_code),
        "sigma_v_volts": float(sd_v),
        "sigma_v_millivolts": float(sd_v * 1e3),
        "sigma_v_lsb": float(sd_v / LSB_V),
        "excess_over_quantisation": float(sd_v / SIGMA_Q_V),
        "bits_lost": bits_lost, "effective_bits": eff_bits,
        "snr_db": snr_db,
        "crb_mhz_at_freq_win": crb_hz * 1e3,
        "drift_codes_over_record": float(drift_code),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("capture")
    ap.add_argument("--compare", default=None,
                    help="second capture (e.g. round-robin) to difference against")
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)

    codes, meta = read_capture(a.capture)
    out = {"primary": analyse(codes, meta, meta.get("label", a.capture))}

    if a.compare:
        c2, m2 = read_capture(a.compare)
        out["compare"] = analyse(c2, m2, m2.get("label", a.compare))
        s1 = out["primary"]["sigma_v_volts"]
        s2 = out["compare"]["sigma_v_volts"]
        extra = math.sqrt(abs(s2 ** 2 - s1 ** 2))
        print(f"{'='*66}")
        print("  DIFFERENCE — what channel switching costs")
        print(f"{'='*66}")
        print(f"    {out['primary']['label']:<28} {s1*1e3:.4f} mV")
        print(f"    {out['compare']['label']:<28} {s2*1e3:.4f} mV")
        print(f"    in quadrature                {extra*1e3:.4f} mV"
              f"  ({extra/LSB_V:.2f} LSB)")
        if s2 > 1.3 * s1:
            print("\n  >> Multiplexer settling dominates. That is a FIRMWARE fix")
            print("     (settling delay or a discarded dummy read per channel),")
            print("     not a hardware limitation. Strong Gate 1 result.")
        else:
            print("\n  >> Channel switching is not the dominant term.")
        print()
        out["mux_contribution_millivolts"] = extra * 1e3

    if a.json:
        with open(a.json, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"  -> {a.json}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
