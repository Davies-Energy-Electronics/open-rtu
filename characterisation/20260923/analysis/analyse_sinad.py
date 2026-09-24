#!/usr/bin/env python3
"""
analyse_sinad.py — IEEE Std 1241 sine-wave test on the open-RTU acquisition chain.

Reads a capture from characterise.ino taken with a full-scale sinusoid applied,
and reports SINAD, ENOB, THD and — the figure the Cramer-Rao bound actually
needs — SNR with the harmonic bins excluded.

Why both:
  IEEE 1241 derives effective bits from SINAD, which folds distortion in with
  noise. A benchtop generator's own distortion is not negligible against a
  12-bit converter (ideal 12-bit SINAD is 74.0 dB, equivalent to 0.02 % THD+N),
  so a SINAD measured this way is a LOWER BOUND on the converter's SINAD and
  the ENOB from it is a lower bound too. Say so when reporting it.

  SNR excluding harmonics is different: the source's distortion lands in known
  harmonic bins, which are excluded from the noise sum, so the figure is
  largely immune to source quality. That is the number to feed to
  crb_analysis.py.

Usage:
    python analyse_sinad.py capture_B1.txt
    python analyse_sinad.py capture_B1.txt --harmonics 12 --json sinad_B1.json

CRB FIGURE — NOT THE BOUND THE PAPER USES. The "CRB at N=200, this SNR" line
printed by this script, and the sidecar field crb_mhz_at_freq_win, use the
coefficient 6 in the Rife-Boorstyn expression. That is the bound for a
COMPLEX exponential. For the real sampled sinusoid the paper analyses the
coefficient is 12 (analysis/crb_analysis.py, confirmed by
analysis/crb_confirm.py), so the figure printed here is LOW by exactly
sqrt(2). It is left unchanged because the archived sidecars were generated
with it; it is not used anywhere in the paper. The Section 3.12 bounds
(0.301 and 0.332 mHz at 61.4 dB) come from analysis/crb_analysis.py.

Standard library plus numpy.
"""
from __future__ import annotations
import argparse, json, math, sys
import numpy as np

LSB_V      = 3.3 / 4095.0
SIGMA_Q_V  = LSB_V / math.sqrt(12)
FREQ_WIN   = 200
FS_NOMINAL = 2000.0


def read_capture(path):
    codes, meta = [], {}
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                b = s.lstrip("#").strip()
                if "," in b:
                    k, _, v = b.partition(",")
                    meta[k.strip()] = v.strip()
                continue
            try:
                codes.append(int(s))
            except ValueError:
                pass
    if not codes:
        sys.exit(f"{path}: no ADC codes found")
    return np.asarray(codes, dtype=float), meta


def blackman_harris_7(n):
    """7-term Blackman-Harris. Sidelobes below -180 dB, so leakage from a
    non-coherent record cannot masquerade as noise."""
    a = [0.27105140069342, -0.43329793923448, 0.21812299954311,
         -0.06592544638803, 0.01081174209837, -0.00077658482522,
          0.00001388721735]
    k = np.arange(n)
    w = np.zeros(n)
    for i, ai in enumerate(a):
        w += ai * np.cos(2 * math.pi * i * k / n)
    return w


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("capture")
    ap.add_argument("--harmonics", type=int, default=10,
                    help="highest harmonic included in THD (default 10)")
    ap.add_argument("--leak", type=int, default=8,
                    help="bins either side of a tone treated as its leakage (default 8)")
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)

    codes, meta = read_capture(a.capture)
    n  = len(codes)
    fs = float(meta.get("fs_hz", FS_NOMINAL))

    v = codes * LSB_V
    v = v - v.mean()

    w = blackman_harris_7(n)
    # Power normalisation, not amplitude normalisation: scaling by N*sum(w^2)
    # makes the sum of psd across a tone's mainlobe equal that tone's
    # mean-square power, so summing over a band is correct. Dividing by the
    # coherent gain instead would overstate a windowed tone by sqrt(N*sum(w^2))
    # / sum(w) — 1.622x for this window. Ratios are unaffected either way;
    # absolute amplitude and absolute noise voltage are not.
    spec = np.fft.rfft(v * w)
    psd  = np.abs(spec) ** 2 / (n * np.sum(w ** 2))
    psd[1:] *= 2.0                         # one-sided
    freqs = np.fft.rfftfreq(n, d=1.0 / fs)

    # ---- fundamental -------------------------------------------------------
    search = psd.copy()
    search[:a.leak + 1] = 0.0              # ignore DC and its skirt
    k_fund = int(np.argmax(search))
    f_fund = freqs[k_fund]

    def band(k):
        lo = max(0, k - a.leak)
        hi = min(len(psd), k + a.leak + 1)
        return slice(lo, hi)

    used = np.zeros(len(psd), dtype=bool)
    used[band(0)] = True                   # DC
    p_fund = psd[band(k_fund)].sum()
    used[band(k_fund)] = True

    # ---- harmonics, folded back across Nyquist ----------------------------
    harm = []
    p_harm = 0.0
    nyq = fs / 2.0
    for h in range(2, a.harmonics + 1):
        fh = h * f_fund
        # alias into the first Nyquist zone
        fa = fh % fs
        if fa > nyq:
            fa = fs - fa
        kh = int(round(fa / (fs / n)))
        if kh <= a.leak or kh >= len(psd) - a.leak:
            continue
        if used[kh]:
            continue
        p = psd[band(kh)].sum()
        used[band(kh)] = True
        p_harm += p
        harm.append({"h": h, "f_hz": float(freqs[kh]), "power": float(p)})

    p_noise = psd[~used].sum()             # noise, harmonics excluded
    p_nad   = p_noise + p_harm             # noise and distortion

    sinad_db = 10 * math.log10(p_fund / p_nad)
    snr_db   = 10 * math.log10(p_fund / p_noise)
    thd_db   = 10 * math.log10(p_harm / p_fund) if p_harm > 0 else float("-inf")
    enob     = (sinad_db - 1.76) / 6.02
    enob_snr = (snr_db - 1.76) / 6.02

    amp   = math.sqrt(2 * p_fund)          # sine amplitude, volts
    sigma_noise = math.sqrt(p_noise)       # rms noise, volts, harmonics excluded

    # CRB at the estimator basis using the harmonics-excluded noise
    rho = amp ** 2 / (2 * sigma_noise ** 2)
    # NOTE: coefficient 6 is the COMPLEX-exponential bound. The paper uses the
    # real-sinusoid bound (coefficient 12, analysis/crb_analysis.py), which is
    # larger by sqrt(2). This figure is kept as archived and is not used by
    # the paper; see the module docstring.
    crb_hz = math.sqrt(6.0 * fs * fs /
                       ((2 * math.pi) ** 2 * rho * FREQ_WIN * (FREQ_WIN ** 2 - 1)))

    print(f"\n{'='*66}")
    print(f"  IEEE 1241 SINE-WAVE TEST — {meta.get('label', a.capture)}")
    print(f"{'='*66}")
    print(f"  samples                {n}   ({n/fs:.3f} s)")
    print(f"  sample rate            {fs:.2f} Hz     bin {fs/n:.5f} Hz")
    print(f"  fundamental            {f_fund:.4f} Hz  (bin {k_fund})")
    print(f"  amplitude              {amp:.4f} V   "
          f"({amp/LSB_V:.0f} codes, {200*amp/3.3:.1f} % of span)")
    print(f"  {'-'*62}")
    print(f"  SINAD                  {sinad_db:6.2f} dB   -> ENOB {enob:5.2f} bits")
    print(f"  SNR, harmonics excl.   {snr_db:6.2f} dB   -> {enob_snr:5.2f} bits")
    print(f"  THD                    {thd_db:6.2f} dB   "
          f"({100*10**(thd_db/20):.3f} %)")
    print(f"  rms noise (excl. harm) {sigma_noise*1e3:.4f} mV  "
          f"= {sigma_noise/LSB_V:.2f} LSB")
    print(f"  {'-'*62}")
    print(f"  CRB at N={FREQ_WIN}, this SNR   {crb_hz*1e3:.4f} mHz")
    print(f"{'='*66}")
    if harm:
        print("  harmonics:")
        for h in harm[:8]:
            print(f"    H{h['h']:<3} {h['f_hz']:9.3f} Hz   "
                  f"{10*math.log10(h['power']/p_fund):7.2f} dBc")
    print()
    print("  Report SINAD/ENOB as a LOWER BOUND unless the source's own")
    print("  distortion has been separately characterised. Report the")
    print("  harmonics-excluded SNR as the Cramer-Rao input.")
    print(f"  Ideal 12-bit: SINAD 74.0 dB, ENOB 12.00, noise "
          f"{SIGMA_Q_V*1e3:.3f} mV.")
    print()

    out = {
        "label": meta.get("label", a.capture), "n": n, "fs_hz": fs,
        "f_fundamental_hz": f_fund, "amplitude_volts": amp,
        "sinad_db": sinad_db, "enob_bits": enob,
        "snr_excl_harmonics_db": snr_db, "enob_from_snr_bits": enob_snr,
        "thd_db": thd_db, "thd_percent": 100 * 10 ** (thd_db / 20),
        "noise_rms_volts_excl_harmonics": sigma_noise,
        "noise_rms_lsb": sigma_noise / LSB_V,
        "crb_mhz_at_freq_win": crb_hz * 1e3,
        "harmonics": harm,
    }
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"  -> {a.json}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
