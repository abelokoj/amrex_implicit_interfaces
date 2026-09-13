"""Export the figure data of a calculation as plain text, once.

`dump_amplitude` writes the amplitude record, which carries the growth curve and the local growth rate. Two further figures need more: the interface figure needs the radius profile r(z) at several instants, and the field figure needs a plotfile variable on the (r,z) plane with the interface laid over it. This program extracts both and writes them in the portable format defined in `bundle.py`, so that every figure of a case can afterwards be drawn from a few hundred kilobytes of text rather than from the plotfiles.

The separation matters because the plotfiles usually stay where the calculation ran. A sweep leaves tens of gigabytes on a cluster that carries no matplotlib, so the figures are drawn on a workstation; moving the plotfiles there is impractical, while moving a bundle is a single copy. Running this program once, on the machine that holds the plotfiles, is what makes that possible, and it is invoked automatically for every case by `scripts/run_all.sh`.

Profiles are exported for every snapshot, since the file remains small. Field slices are exported only at the instants named by `--field-times`, or at three instants spread across the record when none are named, because a field slice is three orders of magnitude larger than a profile and only a few panels are ever drawn.

Usage, from amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Python:

    python3 python/dump_profiles.py <rundir> <wavelength> <outdir>
    python3 python/dump_profiles.py runs/k0.7_eps0.02 8.976 runs/bundles/k0.7 \\
        --field y_velocity --field-times 1.0 4.0 6.7
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import bundle  # noqa: E402
import interface_mode  # noqa: E402
import nsmfp  # noqa: E402


def _axes(pf):
    """Cell-centred radial and axial coordinates of a plotfile."""
    n_r, n_z = pf.shape[0], pf.shape[1]
    r = (np.arange(n_r) + 0.5) * float(pf.domain_right[0]) / n_r
    z = (np.arange(n_z) + 0.5) * float(pf.domain_right[1]) / n_z
    return r, z


def export_profiles(rundir, wavelength, outdir, ls_field="L0101",
                    prefix="nddataPLT"):
    """Write `profiles.dat` for every snapshot of a run."""
    paths = nsmfp.find_plotfiles(rundir, prefix=prefix)
    if not paths:
        raise FileNotFoundError(f"no {prefix}* plotfiles under {rundir}")

    entries = []
    for path in paths:
        t, z, radius = interface_mode.interface_radius(
            path, level_set_field=ls_field)
        entries.append((t, z, radius))
    entries.sort(key=lambda item: item[0])

    z = entries[0][1]
    times = [item[0] for item in entries]
    radii = np.vstack([item[2] for item in entries])
    out = os.path.join(outdir, "profiles.dat")
    bundle.write_profiles(out, wavelength, z, times, radii)
    return out, len(times)


def export_field(rundir, field, instants, outdir, prefix="nddataPLT"):
    """Write `field_<name>.dat` at the snapshots nearest `instants`."""
    paths = nsmfp.find_plotfiles(rundir, prefix=prefix)
    if not paths:
        raise FileNotFoundError(f"no {prefix}* plotfiles under {rundir}")

    # The times are read first from a one-field load so that the selection costs one cheap pass rather than holding every snapshot in memory.
    stamps = []
    for path in paths:
        pf = nsmfp.read_plotfile(path, fields=(field,))
        stamps.append((pf.time, path))
    stamps.sort()
    all_times = np.array([item[0] for item in stamps])

    chosen, seen = [], set()
    for target in instants:
        index = int(np.argmin(np.abs(all_times - target)))
        if index not in seen:
            seen.add(index)
            chosen.append(index)
    chosen.sort()

    planes, times = [], []
    r = z = None
    for index in chosen:
        t, path = stamps[index]
        pf = nsmfp.read_plotfile(path, fields=(field,))
        if r is None:
            r, z = _axes(pf)
        planes.append(pf.data.reshape(pf.shape)[:, :, 0])
        times.append(t)

    out = os.path.join(outdir, f"field_{field}.dat")
    bundle.write_field(out, field, r, z, times, np.stack(planes))
    return out, len(times)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="export interface profiles and field slices as text")
    ap.add_argument("rundir",
                    help="run directory holding the nddataPLT* plotfiles")
    ap.add_argument("wavelength", type=float,
                    help="axial wavelength, equal to ns.yblob of the deck")
    ap.add_argument("outdir", help="destination directory for the bundle")
    ap.add_argument("--field", action="append", default=None,
                    help="plotfile variable to slice; repeatable")
    ap.add_argument("--field-times", type=float, nargs="+",
                    dest="field_times", default=None,
                    help="instants to slice; default three across the record")
    ap.add_argument("--ls-field", dest="ls_field", default="L0101",
                    help="level-set field carrying the interface")
    ap.add_argument("--prefix", default="nddataPLT",
                    help="plotfile directory prefix written by the solver")
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)

    try:
        path, count = export_profiles(args.rundir, args.wavelength,
                                      args.outdir, ls_field=args.ls_field,
                                      prefix=args.prefix)
    except FileNotFoundError as exc:
        print(f"ERROR reading the run: {exc}", file=sys.stderr)
        return 1
    print(f"wrote {count} interface profiles to {path}")

    for field in args.field or []:
        instants = args.field_times
        if instants is None:
            record = bundle.read_profiles(os.path.join(args.outdir,
                                                       "profiles.dat"))
            times = record["times"]
            instants = list(np.linspace(times.min(), times.max(), 3))
        try:
            path, count = export_field(args.rundir, field, instants,
                                       args.outdir, prefix=args.prefix)
        except Exception as exc:                      # noqa: BLE001
            # A missing variable is a naming question rather than a failure of the run, so the profiles already written are kept and the remaining fields are still attempted.
            print(f"ERROR exporting field {field}: {exc}", file=sys.stderr)
            continue
        print(f"wrote {count} slices of {field} to {path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
