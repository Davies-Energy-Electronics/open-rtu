#!/usr/bin/env python3
"""
crb_confirm.py
==============

Three independent evaluations of the Cramer-Rao bound (CRB) on the frequency
of a single REAL sinusoid in white Gaussian noise, as stated in Section 3.5
of the open-RTU MDPI Sensors paper:

  (a) CLOSED FORM. crb_std_hz() imported from crb_analysis.py, i.e.

          Var(f_hat) >= 12 f_s^2 / [ (2 pi)^2 eta N (N^2 - 1) ]      (Hz^2)

      with eta = A^2 / (2 sigma^2). The script checks that the imported
      function equals this expression to rounding.

  (b) FISHER MATRIX. The 3 x 3 Fisher information matrix for the parameters
      (A, f, phi) of s[n] = A cos(2 pi f n / f_s + phi), n = 0 .. N-1, is
      formed from the analytic derivatives,

          I = (1 / sigma^2) J^T J,   J[n, :] = ds[n]/d(A, f, phi),

      and inverted numerically. sqrt([I^-1]_ff) is the bound on the standard
      deviation of f_hat when A and phi are also unknown. It depends slightly
      on the initial phase phi, so it is evaluated on a grid of 360 phases
      (0, 1, ..., 359 degrees) and reported as mean, minimum and maximum.

  (c) MONTE CARLO. For each trial the initial phase is drawn uniformly on
      [0, 2 pi), white Gaussian noise of variance sigma^2 is added, and
      (A, f, phi) are estimated by a nonlinear least-squares single-tone fit
      (scipy.optimize.least_squares, Levenberg-Marquardt, analytic Jacobian,
      started from the true values). Under Gaussian noise that fit is the
      maximum-likelihood estimator, which attains the bound at these SNRs.
      The sample standard deviation of f_hat over the trials is reported with
      its standard error, SE = sd / sqrt(2 (M - 1)) for M trials.

Operating points (N = 200, f = 50 Hz, A = 1):

    f_s = 2000   Hz, 40.0 dB   characterisation loop, conservative SNR
    f_s = 1816.5 Hz, 40.0 dB   canonical build, conservative SNR (Table 4)
    f_s = 1816.5 Hz, 61.4 dB   canonical build, measured SNR (Section 3.12)

The random stream is NumPy's PCG64 seeded with 20260925, so the run is
deterministic. Output values are rounded to six significant figures and
written with sorted keys and LF line endings to crb_confirm_output.json at the
repository root (override with --json). This script is not called by
reproduce_all.py and the paper's numbers do not depend on its output; it
exists to make the Section 3.5 confirmations reproducible.

Usage (from the repository root):
    python analysis/crb_confirm.py
    python analysis/crb_confirm.py --trials 5000 --json crb_confirm_output.json

Requires NumPy and SciPy. Runtime is of the order of a minute.
"""

from __future__ import annotations
import argparse, json, math, os, sys, time

import numpy as np
from scipy.optimize import least_squares

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = (os.path.dirname(_HERE)
         if os.path.isdir(os.path.join(os.path.dirname(_HERE), "analysis"))
         else _HERE)
sys.path.insert(0, _HERE)
from crb_analysis import crb_std_hz   # noqa: E402  (the released closed form)


N_SAMPLES   = 200
F0_HZ       = 50.0
AMPLITUDE   = 1.0
SEED        = 20260925
TRIALS      = 5000
N_PHASES    = 360
OPERATING_POINTS = [          # (label, f_s in Hz, SNR in dB)
    ("2 kHz, 40 dB",        2000.0, 40.0),
    ("1816.5 Hz, 40 dB",    1816.5, 40.0),
    ("1816.5 Hz, 61.4 dB",  1816.5, 61.4),
]


def sig6(x: float) -> float:
    """Round to six significant figures, for a stable JSON."""
    return float(f"{x:.6g}")


def closed_form_reproduced_hz(n: int, fs: float, snr_db: float) -> float:
    """The crb_analysis.py expression written out again, as a check."""
    eta = 10.0 ** (snr_db / 10.0)
    return math.sqrt(12.0 * fs * fs / ((2.0 * math.pi) ** 2 * eta * n * (n * n - 1)))


def model(p, t):
    a, f, ph = p
    return a * np.cos(2.0 * np.pi * f * t + ph)


def jacobian(p, t):
    a, f, ph = p
    th = 2.0 * np.pi * f * t + ph
    s = np.sin(th)
    return np.column_stack((np.cos(th), -a * s * 2.0 * np.pi * t, -a * s))


def fisher_std_hz(n: int, fs: float, snr_db: float, phase: float) -> float:
    """sqrt of the (f, f) element of the inverse 3x3 Fisher matrix, in Hz."""
    sigma2 = AMPLITUDE ** 2 / (2.0 * 10.0 ** (snr_db / 10.0))
    t = np.arange(n) / fs
    J = jacobian((AMPLITUDE, F0_HZ, phase), t)
    fim = (J.T @ J) / sigma2
    return math.sqrt(np.linalg.inv(fim)[1, 1])


def monte_carlo_std_hz(n: int, fs: float, snr_db: float, trials: int,
                       rng: np.random.Generator):
    """Sample SD, mean bias and SE of the SD of f_hat from an LS single-tone fit."""
    sigma = math.sqrt(AMPLITUDE ** 2 / (2.0 * 10.0 ** (snr_db / 10.0)))
    t = np.arange(n) / fs
    f_hat = np.empty(trials)
    for k in range(trials):
        phase = rng.uniform(0.0, 2.0 * np.pi)
        truth = np.array([AMPLITUDE, F0_HZ, phase])
        y = model(truth, t) + sigma * rng.standard_normal(n)
        res = least_squares(lambda p: model(p, t) - y, truth,
                            jac=lambda p: jacobian(p, t), method="lm",
                            xtol=1e-14, ftol=1e-14, gtol=1e-14)
        f_hat[k] = res.x[1]
    sd = float(np.std(f_hat, ddof=1))
    se = sd / math.sqrt(2.0 * (trials - 1))
    bias = float(np.mean(f_hat) - F0_HZ)
    return sd, se, bias


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trials", type=int, default=TRIALS,
                    help=f"Monte-Carlo trials per operating point (default {TRIALS})")
    ap.add_argument("--seed", type=int, default=SEED,
                    help=f"PCG64 seed (default {SEED})")
    ap.add_argument("--json", default=os.path.join(_ROOT, "crb_confirm_output.json"),
                    help="output path (default: crb_confirm_output.json at the repository root)")
    args = ap.parse_args(argv)

    rng = np.random.default_rng(args.seed)
    phases = np.arange(N_PHASES) * (2.0 * np.pi / N_PHASES)
    t0 = time.time()
    all_ok = True
    points = []

    print("=" * 96)
    print("  CRAMER-RAO BOUND, REAL SINUSOID: THREE INDEPENDENT EVALUATIONS (Section 3.5)")
    print(f"  N = {N_SAMPLES}, f = {F0_HZ:.0f} Hz, eta = A^2/(2 sigma^2), "
          f"{args.trials} Monte-Carlo trials per point, seed {args.seed}")
    print("=" * 96)
    print(f"  {'Operating point':<20}{'closed form':>13}{'Fisher mean':>13}"
          f"{'Fisher min-max':>19}{'Monte Carlo':>13}{'MC SE':>9}{'MC/CF':>8}")
    print(f"  {'':<20}{'(mHz)':>13}{'(mHz)':>13}{'(mHz)':>19}{'(mHz)':>13}{'(mHz)':>9}{'':>8}")
    print(f"  {'-'*20}{'-'*13}{'-'*13}{'-'*19}{'-'*13}{'-'*9}{'-'*8}")

    for label, fs, snr_db in OPERATING_POINTS:
        cf = crb_std_hz(N_SAMPLES, fs, snr_db)
        cf_again = closed_form_reproduced_hz(N_SAMPLES, fs, snr_db)
        if abs(cf - cf_again) > 1e-12 * cf:
            print(f"  [FAIL] imported closed form differs from the expression: {cf} vs {cf_again}")
            all_ok = False

        fis = np.array([fisher_std_hz(N_SAMPLES, fs, snr_db, ph) for ph in phases])
        mc_sd, mc_se, mc_bias = monte_carlo_std_hz(N_SAMPLES, fs, snr_db,
                                                   args.trials, rng)
        z = (mc_sd - cf) / mc_se
        fis_mean, fis_min, fis_max = float(fis.mean()), float(fis.min()), float(fis.max())
        agree_fisher = abs(fis_mean - cf) / cf < 0.01
        agree_mc = abs(z) < 3.0
        all_ok &= agree_fisher and agree_mc

        print(f"  {label:<20}{cf*1e3:>13.4f}{fis_mean*1e3:>13.4f}"
              f"{fis_min*1e3:>9.4f} - {fis_max*1e3:<7.4f}{mc_sd*1e3:>13.4f}"
              f"{mc_se*1e3:>9.4f}{mc_sd/cf:>8.3f}")

        points.append({
            "label": label,
            "fs_hz": fs,
            "snr_db": snr_db,
            "closed_form_mhz": sig6(cf * 1e3),
            "closed_form_coefficient_6_mhz": sig6(cf * 1e3 / math.sqrt(2.0)),
            "fisher_mean_over_phase_mhz": sig6(fis_mean * 1e3),
            "fisher_min_over_phase_mhz": sig6(fis_min * 1e3),
            "fisher_max_over_phase_mhz": sig6(fis_max * 1e3),
            "monte_carlo_sd_mhz": sig6(mc_sd * 1e3),
            "monte_carlo_se_mhz": sig6(mc_se * 1e3),
            "monte_carlo_mean_bias_mhz": sig6(mc_bias * 1e3),
            "monte_carlo_z_vs_closed_form": sig6(z),
        })

    print("=" * 96)
    print("  Fisher: sqrt([I^-1]_ff) for (A, f, phi), over 360 initial phases.")
    print("  Monte Carlo: SD of f from a Levenberg-Marquardt single-tone LS fit; SE = sd/sqrt(2(M-1)).")
    print("  The coefficient-6 (complex-exponential) value is closed form / sqrt(2).")
    verdict = ("AGREE (Fisher within 1 %, Monte Carlo within 3 SE of the closed form)"
               if all_ok else "DISAGREE - see above")
    print(f"  Verdict: {verdict}")

    out = {
        "schema_version": "1.0",
        "tool": "crb_confirm.py",
        "config": {
            "n_samples": N_SAMPLES,
            "tone_hz": F0_HZ,
            "amplitude": AMPLITUDE,
            "snr_definition": "eta = A^2 / (2 sigma^2)",
            "fisher_parameters": "A, f, phi",
            "fisher_phase_grid_deg": "0, 1, ..., 359",
            "monte_carlo_trials": args.trials,
            "monte_carlo_seed": args.seed,
            "monte_carlo_rng": "numpy PCG64 (default_rng)",
            "monte_carlo_fit": "scipy.optimize.least_squares, method lm, "
                               "analytic Jacobian, started from the true values",
            "monte_carlo_phase": "uniform on [0, 2 pi) per trial",
            "rounding": "six significant figures",
        },
        "points": points,
        "all_agree": bool(all_ok),
    }
    with open(args.json, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"  -> {args.json}   ({time.time() - t0:.1f} s)")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
