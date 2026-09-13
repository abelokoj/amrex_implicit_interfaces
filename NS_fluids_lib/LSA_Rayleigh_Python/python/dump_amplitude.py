"""Write the interface-mode amplitude record of a calculation as plain text.

The record is the quantity the dispersion relation predicts, namely the amplitude of the k-Fourier component of the interface radius; see `interface_mode.py` for why a field norm is unsuitable for an interfacial instability. Output is three columns: time, mode amplitude, and mean interface radius. The last is included because a drifting mean signals a mass-conservation error, which would invalidate the growth rate.

The format is byte-for-byte the format written by the Fortran tool `bin/dump_amplitude`, so that the records of the two implementations are interchangeable and the plotting scripts need not distinguish them.

Usage, from amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Python:

    python3 python/dump_amplitude.py <rundir> <wavelength> <outfile>
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import interface_mode  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="interface-mode amplitude record as plain text")
    ap.add_argument("rundir",
                    help="run directory holding the nddataPLT* plotfiles")
    ap.add_argument("wavelength", type=float,
                    help="axial wavelength, equal to ns.yblob of the deck")
    ap.add_argument("outfile", help="destination text file")
    ap.add_argument("--field", default="L0101",
                    help="level-set field carrying the interface")
    ap.add_argument("--harmonic", type=int, default=1,
                    help="Fourier harmonic imposed by the deck")
    ap.add_argument("--prefix", default="nddataPLT",
                    help="plotfile directory prefix written by the solver")
    args = ap.parse_args(argv)

    try:
        times, amps, means = interface_mode.amplitude_series(
            args.rundir, args.wavelength, level_set_field=args.field,
            harmonic=args.harmonic, prefix=args.prefix)
    except FileNotFoundError as exc:
        print(f"ERROR reading the run: {exc}", file=sys.stderr)
        return 1

    with open(args.outfile, "w") as fh:
        fh.write("# time  amplitude  mean_radius\n")
        for t, a, r in zip(times, amps, means):
            fh.write(f"{t:24.16E}{a:24.16E}{r:24.16E}\n")

    print(f"wrote {len(times)} rows to {args.outfile}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
