#!/usr/bin/env python3
"""Apply a delivered archive to the working tree, and optionally commit it.

The problem this solves is bookkeeping. Changes have been arriving as a series of patches, each assuming the one before it was applied, and a patch skipped in the middle leaves the tree in a state no later patch matches. This script removes the question: it copies the archive over the working tree, so the result depends on the archive alone and not on which patches were applied along the way.

What it copies is deliberately narrow. Source, scripts, documentation, input decks, tests and sample figures come from the archive. Everything produced by running or building is left alone: the run directories, the exported bundles under results/, the generated figures, the compiled binaries and object files. Some of those represent days of cluster time, and no delivery should overwrite them.

Nothing is written until --apply is given, and nothing reaches GitHub until --push is given as well. The default lists what would change and stops.

    python3 sync_from_archive.py ~/LSA_Rayleigh_complete.zip
    python3 sync_from_archive.py ~/LSA_Rayleigh_complete.zip --apply
    python3 sync_from_archive.py ~/LSA_Rayleigh_complete.zip --apply --push

Only the standard library is used, so this runs on a login node with no environment loaded.
"""

from __future__ import annotations

import argparse
import filecmp
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile

#: Directories that hold generated output. Never copied from the archive, never removed from the working tree.
SKIP_DIRS = {"runs", "results", "figures", "build", "bin", "__pycache__",
             "data", ".pytest_cache", ".git"}
SKIP_SUFFIX = (".pyc", ".o", ".mod", ".a", ".so")

TREES = ("LSA_Rayleigh_Fortran", "LSA_Rayleigh_Python")
TOP_FILES = ("analysis_report.typ", "references.bib")


def digest(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(65536), b""):
            h.update(block)
    return h.hexdigest()


def walk(root):
    """Relative paths of every file under root, minus the generated ones."""
    out = []
    for base, dirs, files in os.walk(root):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        for name in files:
            if name.endswith(SKIP_SUFFIX):
                continue
            full = os.path.join(base, name)
            out.append(os.path.relpath(full, root))
    return sorted(out)


def find_root():
    """The repository root, from AII or by walking up from this file."""
    candidate = os.environ.get("AII", "")
    if candidate and os.path.isdir(os.path.join(candidate, "NS_fluids_lib")):
        return candidate
    here = os.path.dirname(os.path.abspath(__file__))
    while here != "/":
        if os.path.isdir(os.path.join(here, "NS_fluids_lib")):
            return here
        here = os.path.dirname(here)
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("archive")
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--push", action="store_true")
    ap.add_argument("--root", default=None,
                    help="repository root; found automatically if omitted")
    args = ap.parse_args(argv)

    if not os.path.isfile(args.archive):
        ap.error(f"no such archive: {args.archive}")
    root = args.root or find_root()
    if not root:
        ap.error("could not locate NS_fluids_lib; set AII or pass --root")
    dest = os.path.join(root, "NS_fluids_lib")

    stage = tempfile.mkdtemp(prefix="lsa_sync_")
    try:
        with zipfile.ZipFile(args.archive) as zf:
            zf.extractall(stage)

        src = stage
        if not os.path.isdir(os.path.join(src, TREES[0])):
            for base, dirs, _ in os.walk(stage):
                if TREES[0] in dirs:
                    src = base
                    break
        if not os.path.isdir(os.path.join(src, TREES[0])):
            print(f"error: archive does not contain {TREES[0]}", file=sys.stderr)
            return 1

        print(f"archive : {args.archive}")
        print(f"target  : {dest}")
        print(f"mode    : {'APPLY' if args.apply else 'dry run, nothing written'}\n")

        added, changed = [], []
        jobs = []
        for tree in TREES:
            s_root, d_root = os.path.join(src, tree), os.path.join(dest, tree)
            if not os.path.isdir(s_root):
                continue
            for rel in walk(s_root):
                s, d = os.path.join(s_root, rel), os.path.join(d_root, rel)
                shown = f"{tree}/{rel}"
                if not os.path.exists(d):
                    added.append(shown); jobs.append((s, d))
                elif not filecmp.cmp(s, d, shallow=False):
                    changed.append(shown); jobs.append((s, d))
        for name in TOP_FILES:
            s, d = os.path.join(src, name), os.path.join(dest, name)
            if not os.path.isfile(s):
                continue
            if not os.path.exists(d):
                added.append(name); jobs.append((s, d))
            elif not filecmp.cmp(s, d, shallow=False):
                changed.append(name); jobs.append((s, d))

        for label, items in (("ADDED", added), ("CHANGED", changed)):
            print(f"--- {label} ({len(items)}) ---")
            for item in items:
                print(f"    {item}")
            if not items:
                print("    none")
            print()

        if not jobs:
            print("Working tree already matches the archive. Nothing to do.")
            return 0
        if not args.apply:
            print("Nothing was written. Re-run with --apply to make these changes.")
            return 0

        for s, d in jobs:
            os.makedirs(os.path.dirname(d), exist_ok=True)
            shutil.copy2(s, d)
        print(f"wrote {len(jobs)} file(s)")

        # A zip round trip does not preserve the executable bit reliably, and a script without one fails in a way that looks like a missing file.
        for tree in TREES:
            sdir = os.path.join(dest, tree, "scripts")
            if os.path.isdir(sdir):
                for name in os.listdir(sdir):
                    if name.endswith(".sh"):
                        p = os.path.join(sdir, name)
                        os.chmod(p, os.stat(p).st_mode | 0o111)

        print("\n--- verifying ---")
        fort = os.path.join(dest, "LSA_Rayleigh_Fortran")
        if os.path.isfile(os.path.join(fort, "Makefile")):
            subprocess.run(["make", "-s"], cwd=fort, capture_output=True)
            out = subprocess.run(["make", "-s", "test"], cwd=fort,
                                 capture_output=True, text=True)
            line = [l for l in out.stdout.splitlines()
                    if "ALL TEST" in l or "SOME" in l]
            print("   Fortran:", line[0] if line else "check by hand")
        pyt = os.path.join(dest, "LSA_Rayleigh_Python")
        out = subprocess.run([sys.executable, "-m", "pytest", "-q"],
                             cwd=pyt, capture_output=True, text=True)
        print("   Python :", (out.stdout.strip().splitlines() or ["check by hand"])[-1])

        if args.push:
            print("\n--- committing ---")
            subprocess.run(["git", "add", "-A"], cwd=root, check=True)
            staged = subprocess.run(["git", "diff", "--cached", "--quiet"],
                                    cwd=root)
            if staged.returncode == 0:
                print("nothing to commit")
                return 0
            msg = (f"Sync analysis package from delivered archive\n\n"
                   f"Applied {os.path.basename(args.archive)}: "
                   f"{len(added)} file(s) added, {len(changed)} changed. "
                   f"Run output, exported results, generated figures and "
                   f"build artefacts were left untouched.")
            subprocess.run(["git", "commit", "-m", msg], cwd=root, check=True)
            subprocess.run(["git", "push"], cwd=root, check=True)
            print("pushed")
        return 0
    finally:
        shutil.rmtree(stage, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
