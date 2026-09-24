#!/usr/bin/env python3
"""
common_312.py - shared paths, constants and loaders for the Section 3.12
verification scripts (verify_*_312.py, make_figure_19_noise_residual.py, make_figure_20_error_budget.py).

Nothing here is new analysis. The loaders and the estimator are imported
from the released scripts so that every verification runs the same code the
released analysis ran:

  read_codes      analyse_noise.read_capture      (code-only captures)
  read_tcodes     estimate_from_codes.read_log    (t_us,code captures)
  fit_sine, refine_frequency, measure_freq_median
                  analyse_ac.fit_sine / refine_frequency / estimate_freq
                  (estimate_freq is the "exact port of measureFreq()")
  estimate_windows
                  estimate_from_codes.estimate    (mean- or median-pooled
                  offline estimator, float crossing instants, MAX_EDGES=12)

Constants are those of the released scripts and the canonical firmware
(firmware/Feather1codeMimoDemoV2/Feather1codeMimoDemoV2.ino: FREQ_WIN 200,
MAX_EDGES_PER_TYPE 12, analogRead()*3.3/4095).
"""
from __future__ import annotations
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import analyse_noise as _an          # noqa: E402
import analyse_ac as _aac            # noqa: E402
import estimate_from_codes as _efc   # noqa: E402

CAP_0923 = os.path.normpath(os.path.join(HERE, "..", "captures"))
DIR_0913 = os.path.normpath(os.path.join(HERE, "..", "..", "20260913"))
FIG_DIR = os.path.normpath(os.path.join(HERE, "..", "figures"))
OUT_DIR = os.path.join(HERE, "out_312")
os.makedirs(OUT_DIR, exist_ok=True)

LSB_V = 3.3 / 4095.0
SIGMA_Q_V = LSB_V / math.sqrt(12.0)
FREQ_WIN = 200
F0 = 50.0
A_BENCH = 1.0928          # V, conditioned amplitude (README; ac.json 1.09278)

# The three quiet-input records of the 13 September bundle (Fig. 19a, 0.659 mV)
SEPT_QUIET = [os.path.join(DIR_0913, f) for f in
              ("A1_quiet_single.txt", "A2_quiet_single_run2.txt",
               "A3_quiet_single_run3.txt")]
# The two additional single-channel quiet-input records of 23 September
SEPT23_SINGLE = [os.path.join(CAP_0923, f) for f in
                 ("A1_quiet_single_run2.txt", "A1_quiet_single_run3.txt")]
# Round-robin quiet-input records of 23 September (run1 warmed, run2/3 cold)
RR = [os.path.join(CAP_0923, f) for f in
      ("A2_quiet_roundrobin_run1.txt", "A2_quiet_roundrobin_run2.txt",
       "A2_quiet_roundrobin_run3.txt")]
B1 = os.path.join(DIR_0913, "B1_ac_running.txt")

read_codes = _an.read_capture
read_tcodes = _efc.read_log
fit_sine = _aac.fit_sine
refine_frequency = _aac.refine_frequency
measure_freq_median = _aac.estimate_freq
estimate_windows = _efc.estimate


def warmed(path):
    """True if the capture header logs the full 600 s thermal settle."""
    with open(path, errors="replace") as fh:
        return any(l.startswith("# warm 600 s") for l in fh)


def header(path):
    meta = {}
    with open(path, errors="replace") as fh:
        for l in fh:
            if l.startswith("#") and "," in l:
                k, _, v = l.lstrip("#").strip().partition(",")
                meta.setdefault(k.strip(), v.strip())
    return meta


def sd_raw_and_detrended(codes):
    """Raw SD and linear-detrended SD in codes (analyse_noise.py's method)."""
    import numpy as np
    c = np.asarray(codes, dtype=float)
    t = np.arange(len(c))
    p = np.polyfit(t, c, 1)
    return float(c.std(ddof=1)), float((c - np.polyval(p, t)).std(ddof=1))


def rms(xs):
    return math.sqrt(sum(x * x for x in xs) / len(xs))


def mhz_from_sigma_v(sigma_v, amplitude=A_BENCH, f0=F0):
    """characterise_loop.py's propagation: sigma_f = f0^2*sqrt2*sigma_v/S."""
    S = 2 * math.pi * f0 * amplitude
    return f0 * f0 * math.sqrt(2.0) * (sigma_v / S) * 1e3


def save_json(name, obj):
    p = os.path.join(OUT_DIR, name)
    with open(p, "w") as fh:
        json.dump(obj, fh, indent=2)
    print(f"  -> {os.path.relpath(p, HERE)}")
    return p
