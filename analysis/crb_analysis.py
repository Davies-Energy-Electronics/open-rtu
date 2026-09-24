"""
crb_analysis.py
===============

Cramer-Rao Lower Bound (CRB) analysis for the M2 frequency estimator,
reproducing Table 17 of Section 5.3.12 of the open-RTU MDPI Sensors paper.

For an unbiased estimator of the frequency of a single real sinusoid in
additive white Gaussian noise, sampled at f_s over N samples at
signal-to-noise ratio rho (linear, not dB), the Cramer-Rao Lower Bound on
the estimator variance is the standard Rife-Boorstyn result:

    Var(f_hat) >= 12 * f_s^2 / [ (2*pi)^2 * eta * N * (N^2 - 1) ]       (Hz^2)

where eta = A^2 / (2*sigma^2) is the signal-to-noise ratio of a REAL sinusoid.
The CRB standard deviation is sqrt of the above, reported here in mHz.

CORRECTED 23 September 2026 (HQ task L46). Earlier releases used the
coefficient 6, which is the bound for a COMPLEX exponential. A real sampled
sinusoid carries half the Fisher information per sample about frequency, so
the coefficient is 12 and the bound is larger by exactly sqrt(2). Confirmed
three independent ways: the closed form; a Monte-Carlo single-tone
least-squares fit; and numerical inversion of the 3x3 Fisher information
matrix for (A, f, phi), which gives 3.9007 mHz against the closed form's
3.8985 at N = 200, f_s = 2 kHz and 40 dB. The released bound of 2.757 mHz was
therefore not a bound: the Kalman estimator's own Monte-Carlo dispersion sits
below it, which is impossible for a true lower bound.

SAMPLING RATE CORRECTED 23 September 2026 (HQ task L41). Earlier releases used
the 2000 Hz quoted in Sections 2.1.2 and 2.2. That figure is correct for the
deadline-scheduled characterisation loop, which measures 2000.0002 Hz with an
interval standard deviation of 0.0135 us. The canonical M2 build does not run
at that rate: it paces with a fixed 500 us delay AFTER its per-sample work, so
its measured per-channel interval is 550.50 us and its rate 1816.5 Hz, against
the 1855 Hz previously derived in Section 3.12. Twelve captures on the quiet
input and two with the production stimulus agree to 0.04 per cent. The bound is
linear in the sampling rate at fixed N, so the correction lowers it by 9 per
cent and raises every ratio to it, in the conservative direction.

This script:
  1. Evaluates the CRB standard deviation at two bracketing SNR scenarios:
       - conservative: SNR = 40 dB (practical measurement conditions)
       - ideal:        SNR = 74 dB (ADC SINAD ceiling, Section 4.1.4)
  2. Compares each of the four algorithms' achieved standard deviation
     (from Table 16) against the conservative-scenario CRB, reporting the
     ratio-to-CRB that quantifies headroom for further refinement.

The four achieved-std inputs default to the Table 16 values but can be
overridden from a hybrid/kalman JSON sidecar via --from-json so the wrapper
can chain real (not hardcoded) numbers.

Usage:
    python crb_analysis.py
    python crb_analysis.py --n 200 --fs 2000
    python crb_analysis.py --v1 130.85 --v2 208.93 --v3 116.75 --v4 202.52
    python crb_analysis.py --json crb_output.json

Writes crb_output.json. Standard library only. No external dependencies.

Reproducibility:
    Defaults reproduce Table 4 exactly:
      CRB std (40 dB) = 3.541 mHz ; CRB std (74 dB) = 0.071 mHz
      ratios: V1 37.0x, V2 59.0x, V3 33.0x, V4 54.9x
    (Superseded values, for anyone comparing against earlier deposits:
     2.757 mHz with 47.5x/75.8x/42.4x/70.5x, from the complex-exponential
     coefficient; and 3.899 mHz with 33.6x/53.6x/29.9x/49.8x, which used the
     corrected coefficient but the assumed 2 kHz sampling rate.)

Author: Jack Davies
"""

from __future__ import annotations
import argparse, json, math, os


# ============================================================================
# Constants (paper defaults)
# ============================================================================

N_DEFAULT       = 200          # samples per M2 window
FS_DEFAULT      = 1816.5       # Hz, M2 per-channel rate, MEASURED 23 Sep 2026
SNR_CONSERV_DB  = 40.0         # practical measurement conditions
SNR_IDEAL_DB    = 74.0         # ADC SINAD ceiling (Section 4.1.4)

# Table 16 achieved standard deviations (mHz) - defaults
STD_DEFAULTS_MHZ = {
    "V1 baseline":  130.85,
    "V2 canonical": 208.93,
    "V3 hybrid":    116.75,
    "V4 Kalman":    194.25,   # TRUE Kalman output (was estimated 202.52 in V8)
}


# ============================================================================
# CRB
# ============================================================================

def crb_std_hz(n: int, fs: float, snr_db: float) -> float:
    """Cramer-Rao Lower Bound standard deviation (Hz) for single-tone
    frequency estimation in AWGN, REAL sinusoid (Rife & Boorstyn 1974;
    Kay 1993 eq. 3.41). Coefficient 12, not 6 - see the module docstring."""
    rho = 10.0 ** (snr_db / 10.0)
    var = 12.0 * fs * fs / ((2.0 * math.pi) ** 2 * rho * n * (n * n - 1))
    return math.sqrt(var)


# ============================================================================
# Main
# ============================================================================

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=N_DEFAULT,
                    help=f"samples per window (default {N_DEFAULT})")
    ap.add_argument("--fs", type=float, default=FS_DEFAULT,
                    help=f"sample rate Hz (default {FS_DEFAULT})")
    ap.add_argument("--snr-conservative", type=float, default=SNR_CONSERV_DB)
    ap.add_argument("--snr-ideal", type=float, default=SNR_IDEAL_DB)
    ap.add_argument("--v1", type=float, default=STD_DEFAULTS_MHZ["V1 baseline"])
    ap.add_argument("--v2", type=float, default=STD_DEFAULTS_MHZ["V2 canonical"])
    ap.add_argument("--v3", type=float, default=STD_DEFAULTS_MHZ["V3 hybrid"])
    ap.add_argument("--v4", type=float, default=STD_DEFAULTS_MHZ["V4 Kalman"])
    ap.add_argument("--from-json", type=str, default=None,
                    help="Optional hybrid/kalman JSON to pull achieved std values from")
    ap.add_argument("--json", type=str, default="crb_output.json")
    args = ap.parse_args(argv)

    achieved = {
        "V1 baseline":  args.v1,
        "V2 canonical": args.v2,
        "V3 hybrid":    args.v3,
        "V4 Kalman":    args.v4,
    }

    # Optional override from a chained JSON (e.g. kalman_output.json for V4)
    if args.from_json and os.path.isfile(args.from_json):
        with open(args.from_json, encoding="utf-8") as fh:
            j = json.load(fh)
        # Accept a few plausible shapes without being fragile
        for key, std in (j.get("achieved_std_mhz", {}) or {}).items():
            if key in achieved:
                achieved[key] = std

    crb_conserv = crb_std_hz(args.n, args.fs, args.snr_conservative) * 1000.0  # mHz
    crb_ideal   = crb_std_hz(args.n, args.fs, args.snr_ideal) * 1000.0         # mHz

    print()
    print("=" * 70)
    print("  CRAMER-RAO LOWER BOUND ANALYSIS (Table 17)")
    print("=" * 70)
    print(f"  N = {args.n} samples,  f_s = {args.fs:.0f} Hz")
    print(f"  CRB std (conservative, {args.snr_conservative:.0f} dB SNR): "
          f"{crb_conserv:.3f} mHz")
    print(f"  CRB std (ideal,        {args.snr_ideal:.0f} dB SNR): "
          f"{crb_ideal:.3f} mHz")
    print()
    print(f"  {'Algorithm':<16}{'Achieved std [mHz]':>20}{'CRB std [mHz]':>16}"
          f"{'Ratio to CRB':>16}")
    print(f"  {'-'*16}{'-'*20}{'-'*16}{'-'*16}")

    rows = []
    for name in ["V1 baseline", "V2 canonical", "V3 hybrid", "V4 Kalman"]:
        std = achieved[name]
        ratio = std / crb_conserv
        print(f"  {name:<16}{std:>20.2f}{crb_conserv:>16.3f}{ratio:>15.1f}x")
        rows.append({"algorithm": name,
                     "achieved_std_mhz": std,
                     "crb_std_mhz": crb_conserv,
                     "ratio_to_crb": ratio})
    print("=" * 70)

    out = {
        "schema_version": "1.0",
        "tool": "crb_analysis.py",
        "config": {"n": args.n, "fs_hz": args.fs,
                   "snr_conservative_db": args.snr_conservative,
                   "snr_ideal_db": args.snr_ideal},
        "crb_std_conservative_mhz": crb_conserv,
        "crb_std_ideal_mhz": crb_ideal,
        "table_17_rows": rows,
    }
    # newline="\n" so the SHA-256 is identical on every platform.
    with open(args.json, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(out, fh, indent=2)
    print(f"  -> {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
