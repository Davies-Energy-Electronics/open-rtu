#!/usr/bin/env python3
"""
make_hashes.py — write or verify hashes.txt, the release manifest.

hashes.txt is what lets a reader confirm that the bundle they downloaded is the
bundle the paper describes. It went stale the moment the analysis layer moved on
from the v1.0.0 release, which is one half of HQ task L42; this script exists so
that it cannot go stale silently again.

    python analysis/make_hashes.py            # rewrite hashes.txt
    python analysis/make_hashes.py --check    # verify, change nothing, exit 1 on drift

Run from anywhere: the repository root is found from this file's own location.

Every hash is taken over the file's bytes exactly as stored. Combined with the
.gitattributes rule that disables line-ending conversion, that makes the manifest
reproducible on Windows, macOS and Linux from the same clone.

Jack Davies, September 2026
"""
from __future__ import annotations
import argparse, hashlib, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = (os.path.dirname(HERE)
        if os.path.isdir(os.path.join(os.path.dirname(HERE), "analysis"))
        else HERE)

MANIFEST = os.path.join(ROOT, "hashes.txt")

# Directories that are part of the release. Anything outside them is working
# material and is deliberately not manifested.
INCLUDE_DIRS = ("analysis", "captures", "trajectories", "firmware", "hmi", "figures")
INCLUDE_ROOT_FILES = ("README.md", "LICENSE", "CITATION.cff", ".gitattributes")

SKIP_NAMES = {".DS_Store", "Thumbs.db", "hashes.txt"}
SKIP_DIRS = {".git", "__pycache__", ".ipynb_checkpoints", ".vscode", ".idea"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def walk():
    out = []
    for name in INCLUDE_ROOT_FILES:
        p = os.path.join(ROOT, name)
        if os.path.isfile(p):
            out.append(name)
    for d in INCLUDE_DIRS:
        base = os.path.join(ROOT, d)
        if not os.path.isdir(base):
            continue
        for dirpath, dirnames, filenames in os.walk(base):
            dirnames[:] = [x for x in dirnames if x not in SKIP_DIRS]
            for fn in sorted(filenames):
                if fn in SKIP_NAMES:
                    continue
                full = os.path.join(dirpath, fn)
                out.append(os.path.relpath(full, ROOT).replace(os.sep, "/"))
    return sorted(set(out))


def read_manifest():
    if not os.path.isfile(MANIFEST):
        return {}
    d = {}
    with open(MANIFEST, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.replace("  ", " ").split(" ", 1)
            if len(parts) != 2:
                continue
            a, b = parts[0].strip(), parts[1].strip().lstrip("*")
            # tolerate both "<hash> <path>" and "<path> <hash>"
            if len(a) == 64 and all(c in "0123456789abcdefABCDEF" for c in a):
                d[b] = a.lower()
            elif len(b) == 64:
                d[a] = b.lower()
    return d


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="verify without writing; exit 1 if anything has drifted")
    a = ap.parse_args(argv)

    files = walk()
    fresh = {f: sha256(os.path.join(ROOT, f)) for f in files}

    if not a.check:
        with open(MANIFEST, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("# SHA-256 manifest of the open-rtu release.\n")
            fh.write("# Regenerate with:  python analysis/make_hashes.py\n")
            fh.write("# Verify with:      python analysis/make_hashes.py --check\n")
            for f in files:
                fh.write(f"{fresh[f]}  {f}\n")
        print(f"  wrote {MANIFEST}  ({len(files)} files)")
        return 0

    old = read_manifest()
    changed = [f for f in files if f in old and old[f] != fresh[f]]
    missing = [f for f in old if f not in fresh]
    unlisted = [f for f in files if f not in old]

    for f in changed:
        print(f"  CHANGED   {f}")
    for f in sorted(missing):
        print(f"  MISSING   {f}   (in the manifest, not on disk)")
    for f in unlisted:
        print(f"  UNLISTED  {f}   (on disk, not in the manifest)")

    if changed or missing or unlisted:
        print(f"\n  {len(changed)} changed, {len(missing)} missing, "
              f"{len(unlisted)} unlisted. Run without --check to rewrite.")
        return 1
    print(f"  hashes.txt matches all {len(files)} files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
