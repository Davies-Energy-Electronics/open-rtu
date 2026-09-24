"""
regenerate_figure_9_overlays.py
================================

Regenerates the HMI's replay-overlay PNGs from archived and live capture
CSVs, using `generate_overlay_plot()` from SCADA_HMI.py (|error| as
percentage of nominal 50 Hz, the +-10 mHz and +-100 mHz lines, and five
region-shaded performance bands). The file name is historical: the overlay
was the replay figure of an earlier draft (numbered Figure 9 there). The
manuscript's own replay figures (Figures 11 and 15-17) are drawn by
paper_figures.py, not by this script.

This exists because SCADA_HMI.py itself only produces the overlay at the end
of a live replay session; regenerating from archived CSVs would otherwise
require a bench run. This script imports `generate_overlay_plot()` directly
and calls it against each provided CSV, so the overlays can be
refreshed without touching the boards.

Usage:
    python regenerate_figure_9_overlays.py
    (regenerates against the four default CSVs listed below)

    python regenerate_figure_9_overlays.py file1.csv file2.csv ...
    (regenerates against the CSVs you specify)

Requires SCADA_HMI.py, neso_trajectory.csv, and the CSV files themselves to
sit in the same directory. Matplotlib must be installed (pip install
matplotlib). Tkinter is NOT required - this script does not open the GUI.

Author: Jack Davies
"""

from __future__ import annotations
import os
import sys


# Default set: the two archived canonical captures (V2 and the V1 baseline),
# plus the two 7 July 2026 live captures that back Table 6.
DEFAULT_CSVS = [
    "replay_20260627_102507V2.csv",   # archived V2 canonical (cf. Figure 11)
    "replay_20260625_180520V1.csv",   # archived V1 baseline  (Table 4, V1 row)
    "replay_20260707_132321.csv",     # live V2 (7 Jul, Table 6)
    "replay_20260707_133804V1.csv",   # live V1 (7 Jul, Table 6)
]


def _import_overlay_function():
    """Import generate_overlay_plot() from SCADA_HMI.py WITHOUT running its
    main() or importing tkinter. Tkinter is not available in headless Python
    environments, and we do not need the GUI here."""
    # Guard against Tk not being available on headless machines
    try:
        import tkinter  # noqa: F401
        # If tkinter imports fine, we can just import SCADA_HMI normally
        from SCADA_HMI import generate_overlay_plot
        return generate_overlay_plot
    except ImportError:
        # Fall back: exec just the two functions we need
        import csv as _csv, math as _math
        src = open("SCADA_HMI.py", encoding="utf-8").read()
        ns = {"csv": _csv, "math": _math, "os": os}
        # load_neso_reference
        i = src.find("def load_neso_reference")
        j = src.find("\ndef ", i + 1)
        exec(src[i:j], ns)
        # generate_overlay_plot
        i = src.find("def generate_overlay_plot")
        j = src.find("\ndef ", i + 1)
        exec(src[i:j], ns)
        return ns["generate_overlay_plot"]


def main(argv=None):
    argv = argv or sys.argv[1:]
    csvs = argv if argv else DEFAULT_CSVS

    if not os.path.isfile("neso_trajectory.csv"):
        print("ERROR: neso_trajectory.csv not found in the current directory.")
        print("       This file provides the ground-truth reference for the")
        print("       overlay; move it into this folder and try again.")
        return 2

    overlay = _import_overlay_function()

    print()
    print("=" * 68)
    print("  Regenerating replay overlays with percentage error scale")
    print("  and region-shaded performance bands")
    print("=" * 68)

    missing = []
    written = []
    for path in csvs:
        if not os.path.isfile(path):
            missing.append(path)
            print(f"  [skip] not found: {path}")
            continue
        print(f"\n  Rendering: {path}")
        png_path = overlay(path, "neso_trajectory.csv")
        if png_path:
            written.append(png_path)

    print()
    print("=" * 68)
    print(f"  {len(written)} PNG(s) written, {len(missing)} skipped.")
    for p in written:
        print(f"    -> {p}")
    if missing:
        print(f"  Missing input CSVs (not fatal - skipped):")
        for p in missing:
            print(f"    -> {p}")
    print("=" * 68)
    return 0 if written else 1


if __name__ == "__main__":
    raise SystemExit(main())
