"""
tost_metrics.py
===============

Two One-Sided Tests (TOST) for measurement equivalence of the open RTU's
M2 frequency output against the NESO reference, with a non-parametric
bootstrap interval for the mean error in each window.

Implements:
    1. Schuirmann (1987) TOST procedure with the equivalence margin delta and
       significance level alpha set explicitly. Reports p-values and a verdict
       for each of three windows of the canonical replay capture (pre-event
       nominal, post-event settled, full replay).
    2. Non-parametric (i.i.d.) percentile bootstrap 95% confidence interval
       for the mean error in each window, computed via 10,000 resamples with
       replacement, compared with the asymptotic Wald interval
       (mean +/- 1.96 SE). Agreement between the two 95% intervals shows that
       the sampling distribution of the mean is close to normal near the
       centre. It does NOT confirm the p-values in the far tail, and the
       bootstrap resamples frames independently, so it does not address
       serial correlation.

NORMAL APPROXIMATION. The Student-t CDF in the parametric TOST step is
replaced by the standard-normal CDF (math.erfc). This is mildly
anti-conservative: the normal tail is thinner than the t tail, so the p-values
printed here are slightly smaller than the exact t values near the verdict
threshold, and many orders of magnitude smaller at the 100 mHz margin, where
|t| is large. In addition, the upper-tail probability is formed as
1 - CDF(t1), which cancels to exactly 0 in double precision for large t1, so
the 100 mHz p-values printed here are not even accurate normal-tail values.
Use this script for its verdicts only; they agree with the exact test on all
three windows. The exact Student-t p-values of Table 2 of the paper are
computed by paper_tables.py (function tost_exact), not by this script.

TIER NAMES. The console labels the margins Tier 1 (±10 mHz) and Tier 2
(±100 mHz), the declared accuracy tiers of the paper (Section 2.8); they are
not IEC 61400-21 classes. The JSON sidecar keeps its historical window keys
"Class I (±10 mHz)" and "Class II (±100 mHz)" so that the released sidecar is
reproduced byte for byte and reproduce_all.py can read it; those keys hold the
Tier 1 and Tier 2 results respectively.

Usage:
    python tost_metrics.py replay_20260627_102507V2.csv
    python tost_metrics.py replay_20260627_102507V2.csv --no-bootstrap
    python tost_metrics.py replay_20260627_102507V2.csv --bootstrap-seed 42

Writes a JSON sidecar <csv>_tost.json alongside the input CSV, overwriting
any existing file of that name. For the canonical capture that is the
manifested captures/replay_20260627_102507V2_tost.json. Any run with a
non-default --alpha, --bootstrap-seed, --bootstrap-n or --no-bootstrap should
therefore be given --json with a path outside captures/, for example
    python tost_metrics.py captures/replay_20260627_102507V2.csv \
        --bootstrap-seed 42 --json /tmp/tost_seed42.json
so that the released sidecar is not replaced.

Reproducibility:
    Canonical input:  replay_20260627_102507V2.csv
    Expected SHA-256: dadbcfe247dd4e538b71a891b2b26f070d02c42e3dd4622b5d2125aff169560f
    Bootstrap seed:   20260628 (project-wide reproducibility convention)
    Bootstrap N:      10000 resamples
    Bootstrap CI:     percentile method, alpha = 0.05 (95% CI)

Standard library only. No external dependencies.

References:
    Schuirmann, D.J. (1987). J. Pharmacokinet. Biopharm. 15:657-680.
    Bland, J.M. and Altman, D.G. (1986). Lancet 327:307-310.
    Efron, B. and Tibshirani, R.J. (1993). An Introduction to the Bootstrap.
"""

from __future__ import annotations
import argparse, csv, json, math, os, random, statistics


# ============================================================================
# Constants
# ============================================================================

BOOTSTRAP_N_RESAMPLES = 10000
BOOTSTRAP_SEED_DEFAULT = 20260628
SCHEMA_VERSION = "1.1"


# ============================================================================
# Data loading (unchanged from v1.0)
# ============================================================================

def load_errors(path: str):
    """Returns three dicts {window_label: [errors_mhz]} for pre, settled, full."""
    rows = []
    with open(path, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            try:
                i = int(r['replay_index'])
                if i < 0:
                    continue
                err = (float(r['f_Vb']) - float(r['f_target_hz'])) * 1000.0
                rows.append((i, err))
            except (KeyError, ValueError):
                continue
    pre  = [e for (i, e) in rows if   0 <= i <= 45]
    sett = [e for (i, e) in rows if 301 <= i <= 359]
    full = [e for (i, e) in rows]
    return {'pre_event_nominal': pre,
            'post_event_settled': sett,
            'full_replay_window': full}


# ============================================================================
# Parametric TOST (unchanged from v1.0; Schuirmann 1987)
# ============================================================================

def t_cdf_normal_approx(t: float, df: int) -> float:
    """Standard-normal approximation of the t-distribution CDF.
    The df argument is ignored. The approximation is mildly
    anti-conservative (p-values slightly too small) and poor in the far tail;
    see the module docstring. The bootstrap below compares 95% intervals only
    and does not validate far-tail p-values."""
    return 0.5 * math.erfc(-t / math.sqrt(2))


def tost(errs, delta_mhz: float, alpha: float = 0.05):
    """Returns dict {n, mean, sd, se, t1, p1, t2, p2, equivalent}."""
    n = len(errs)
    if n < 2:
        return None
    mean = statistics.mean(errs)
    sd = statistics.stdev(errs)
    se = sd / math.sqrt(n)
    df = n - 1
    t1 = (mean - (-delta_mhz)) / se
    p1 = 1.0 - t_cdf_normal_approx(t1, df)
    t2 = (mean - delta_mhz) / se
    p2 = t_cdf_normal_approx(t2, df)
    return {'n': n, 'mean': mean, 'sd': sd, 'se': se, 'df': df,
            't1': t1, 'p1': p1, 't2': t2, 'p2': p2,
            'equivalent': (p1 < alpha) and (p2 < alpha)}


# ============================================================================
# Non-parametric bootstrap interval (NEW in v1.1)
# ============================================================================

def bootstrap_ci_mean(errs, n_resamples: int = BOOTSTRAP_N_RESAMPLES,
                      alpha: float = 0.05,
                      seed: int = BOOTSTRAP_SEED_DEFAULT):
    """
    Non-parametric percentile bootstrap for the 100*(1-alpha)% confidence
    interval of the mean of `errs` (a list of frequency-error values in mHz).

    Resamples `errs` with replacement n_resamples times, computes the sample
    mean of each resample, and returns the alpha/2 and 1-alpha/2 quantiles of
    the resulting bootstrap distribution. The percentile method makes no
    normality or symmetry assumption, but resampling is i.i.d., so it assumes
    independent frames. Agreement with the Wald interval shows that the
    sampling distribution of the mean is close to normal at the 95% level; it
    says nothing about p-values in the far tail.

    Returns a dict with all summary statistics plus an explicit comparison
    against the asymptotic Wald interval (mean +/- z_{1-alpha/2} * SE), so
    the agreement between the two methods is reported on the same JSON entry.
    """
    n = len(errs)
    if n < 2:
        return None

    rng = random.Random(seed)
    boot_means = []
    for _ in range(n_resamples):
        sample = rng.choices(errs, k=n)
        boot_means.append(sum(sample) / n)
    boot_means.sort()

    lo_idx = int(math.floor(alpha / 2.0 * n_resamples))
    hi_idx = min(int(math.ceil((1.0 - alpha / 2.0) * n_resamples)) - 1,
                 n_resamples - 1)
    ci_lower = boot_means[lo_idx]
    ci_upper = boot_means[hi_idx]
    ci_halfwidth = (ci_upper - ci_lower) / 2.0

    observed_mean = statistics.mean(errs)
    sd = statistics.stdev(errs)
    se = sd / math.sqrt(n)
    z = 1.959963984540054
    asym_lower = observed_mean - z * se
    asym_upper = observed_mean + z * se
    asym_halfwidth = z * se

    dev_lower = abs(ci_lower - asym_lower)
    dev_upper = abs(ci_upper - asym_upper)
    max_dev = max(dev_lower, dev_upper)
    max_dev_pct_of_asym_halfwidth = (100.0 * max_dev / asym_halfwidth
                                     if asym_halfwidth > 0 else 0.0)

    # Agreement verdict: pass on either absolute or relative criterion.
    # The absolute criterion captures the physical-relevance of the deviation
    # (sub-millihertz is far below M2's resolution); the relative criterion
    # captures the statistical-shape agreement. Small-n windows may show
    # higher relative deviation despite small absolute deviation, so passing
    # either criterion is sufficient evidence of agreement.
    if max_dev < 1.0:
        verdict = 'agreement (sub-mHz absolute)'
    elif max_dev_pct_of_asym_halfwidth < 5.0:
        verdict = 'agreement (<5% relative)'
    else:
        verdict = 'check'

    return {
        'n': n,
        'n_resamples': n_resamples,
        'seed': seed,
        'alpha': alpha,
        'observed_mean_mhz': observed_mean,
        'bootstrap_mean_of_means_mhz': sum(boot_means) / n_resamples,
        'bootstrap_ci_lower_mhz': ci_lower,
        'bootstrap_ci_upper_mhz': ci_upper,
        'bootstrap_ci_halfwidth_mhz': ci_halfwidth,
        'asymptotic_ci_lower_mhz': asym_lower,
        'asymptotic_ci_upper_mhz': asym_upper,
        'asymptotic_ci_halfwidth_mhz': asym_halfwidth,
        'max_deviation_mhz': max_dev,
        'max_deviation_pct_of_asym_halfwidth': max_dev_pct_of_asym_halfwidth,
        'agreement_verdict': verdict,
    }


# ============================================================================
# Main
# ============================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('csv', help="Replay CSV from SCADA_HMI.py")
    ap.add_argument('--alpha', type=float, default=0.05,
                    help='Significance level (default 0.05)')
    ap.add_argument('--no-bootstrap', action='store_true',
                    help='Skip bootstrap (v1.0-compatible output)')
    ap.add_argument('--bootstrap-seed', type=int, default=BOOTSTRAP_SEED_DEFAULT,
                    help=f'PRNG seed (default {BOOTSTRAP_SEED_DEFAULT})')
    ap.add_argument('--bootstrap-n', type=int, default=BOOTSTRAP_N_RESAMPLES,
                    help=f'Resample count (default {BOOTSTRAP_N_RESAMPLES})')
    ap.add_argument('--json', default=None,
                    help='Output path for the JSON sidecar (default: '
                         '<csv>_tost.json next to the input, which overwrites '
                         'the manifested sidecar for the canonical capture)')
    args = ap.parse_args()

    windows = load_errors(args.csv)
    results = {
        'csv': os.path.basename(args.csv),
        'schema_version': SCHEMA_VERSION,
        'alpha': args.alpha,
        'windows': {},
    }

    print()
    print("=== TOST equivalence testing + bootstrap interval ===")
    print(f"  CSV:                 {os.path.basename(args.csv)}")
    print(f"  alpha (TOST and CI): {args.alpha}")
    print("  p-values:            standard-normal approximation (verdicts only;"
          " exact t: paper_tables.py)")
    if not args.no_bootstrap:
        print(f"  bootstrap N:         {args.bootstrap_n} resamples")
        print(f"  bootstrap seed:      {args.bootstrap_seed}")
    print()

    for label, errs in windows.items():
        print(f"  {label}  (n = {len(errs)}):")
        win_result = {}

        # band_name is the historical JSON key and is kept unchanged so that
        # the released sidecar is reproduced; tier_label is what is printed.
        for delta, band_name, tier_label in [
                (10, 'Class I (\u00b110 mHz)', 'Tier 1 (\u00b110 mHz)'),
                (100, 'Class II (\u00b1100 mHz)', 'Tier 2 (\u00b1100 mHz)')]:
            r = tost(errs, delta, args.alpha)
            if r is None:
                continue
            verdict = "EQUIVALENT" if r['equivalent'] else "NOT EQUIVALENT"
            print(f"    {tier_label}:  mean = {r['mean']:+.2f} mHz, "
                  f"sd = {r['sd']:.2f} mHz,  "
                  f"p1 = {r['p1']:.4g}, p2 = {r['p2']:.4g}  --> {verdict}")
            win_result[band_name] = r

        if not args.no_bootstrap:
            b = bootstrap_ci_mean(errs,
                                  n_resamples=args.bootstrap_n,
                                  alpha=args.alpha,
                                  seed=args.bootstrap_seed)
            if b is not None:
                print(f"    Bootstrap 95% CI for mean:  "
                      f"[{b['bootstrap_ci_lower_mhz']:+.2f}, "
                      f"{b['bootstrap_ci_upper_mhz']:+.2f}] mHz")
                print(f"    Asymptotic 95% CI for mean: "
                      f"[{b['asymptotic_ci_lower_mhz']:+.2f}, "
                      f"{b['asymptotic_ci_upper_mhz']:+.2f}] mHz")
                print(f"    Max deviation:              "
                      f"{b['max_deviation_mhz']:.3f} mHz "
                      f"({b['max_deviation_pct_of_asym_halfwidth']:.2f}% of "
                      f"asymptotic half-width)  --> {b['agreement_verdict'].upper()}")
                win_result['bootstrap_ci_mean'] = b

        results['windows'][label] = win_result
        print()

    out_path = args.json or (os.path.splitext(args.csv)[0] + '_tost.json')
    with open(out_path, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(results, fh, indent=2)
    print(f"  -> {out_path}")


if __name__ == '__main__':
    main()
