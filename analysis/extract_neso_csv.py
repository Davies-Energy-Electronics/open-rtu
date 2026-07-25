"""
extract_neso_csv.py
===================

Pulls the 360-element NESO_TRAJECTORY[] array out of neso_trajectory.h and
writes it as neso_trajectory.csv (index,frequency_hz) for SCADA_HMI.py's
sanity-check overlay to load.

Run once. The CSV does not change unless the .h file is regenerated.

Usage:
    python extract_neso_csv.py                       # uses neso_trajectory.h in cwd
    python extract_neso_csv.py path/to/header.h      # explicit path
"""

import re
import sys
import csv
from pathlib import Path


def extract(h_path: Path, csv_path: Path) -> None:
    text = h_path.read_text(encoding='utf-8')

    # Find the curly-brace body containing the array initialiser. The header
    # comment in neso_trajectory.h declares
    #     const float NESO_TRAJECTORY[NESO_TRAJECTORY_LENGTH] = {
    #         50.045f,  50.064f,  ...
    #     };
    m = re.search(r'NESO_TRAJECTORY\s*\[[^\]]*\]\s*=\s*\{([^}]*)\}',
                  text, flags=re.DOTALL)
    if not m:
        raise SystemExit(f"Could not find NESO_TRAJECTORY array body in {h_path}")
    body = m.group(1)

    # Pull every float literal (with or without trailing 'f' / 'F')
    floats = re.findall(r'[+-]?\d+\.\d+(?:[eE][+-]?\d+)?[fF]?', body)
    if not floats:
        raise SystemExit("No float literals found inside NESO_TRAJECTORY body")
    values = [float(v.rstrip('fF')) for v in floats]

    # Sanity: expected length is 360 (window −60 s … +299 s around nadir)
    if len(values) != 360:
        print(f"WARN: expected 360 elements, got {len(values)}. "
              f"Writing CSV anyway - check the header file.")

    with csv_path.open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['index', 'frequency_hz'])
        for i, f in enumerate(values):
            w.writerow([i, f'{f:.6f}'])

    nadir_i = min(range(len(values)), key=lambda i: values[i])
    print(f"Wrote {csv_path} - {len(values)} samples, "
          f"nadir {values[nadir_i]:.4f} Hz at index {nadir_i}, "
          f"range {min(values):.3f} - {max(values):.3f} Hz")


if __name__ == '__main__':
    h_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('neso_trajectory.h')
    csv_path = h_path.with_name('neso_trajectory.csv')
    if not h_path.is_file():
        raise SystemExit(f"Header file not found: {h_path}")
    extract(h_path, csv_path)
