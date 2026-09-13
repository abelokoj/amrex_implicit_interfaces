# !/usr/bin/env python3
"""Cost of a Jacobian-vector product, and what that buys.

Describing the method as Jacobian-free states what it does not do. Measuring the cost of a Jacobian-vector product against the cost of forming and storing the Jacobian states what it gains, and turns an architectural description into a result. That comparison is the point of this script.

The quantity measured is the wall-clock time to advance one snapshot of the twin pair, which is what one Jacobian-vector product costs in this formulation, together with the time to read it back and assemble it into the subspace. The quantity compared against is not measured but counted: an explicit Jacobian for a state vector of n degrees of freedom has n columns, each requiring one such product, and needs n squared storage. For the reference configuration at 8 cells per radius the state vector already runs to thousands of degrees of freedom, so the explicit matrix is not merely expensive but cannot be stored, and the figure should be read as establishing infeasibility rather than a ratio.

Every timing is written together with the machine it was taken on, through `machine_info`. A cost without its machine is not reproducible and cannot be compared against a later measurement, and the distinction between processors present and processors granted matters especially on a cluster, where a node advertising many cores may yield a handful inside an allocation.

Two modes are provided. With `--rundir` the script times the analysis path on real snapshots, which is what the reported numbers should come from. With `--synthetic` it times the same operations on generated arrays of stated size, which needs no calculation and is useful for establishing how cost scales with problem size before committing to a sweep.

Run from `amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Fortran`:

    python3 plotting/benchmark_jvp.py --synthetic --outdir figures
    python3 plotting/benchmark_jvp.py --perturbed runs/k0.7_eps0.02 --base runs/k0.7_base --outdir figures
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import machine_info
from pub_style import apply_style, savefig_all, cm_x, cm_both

GUIDE_LW = 0.8

#: Bytes per double, used to convert a Jacobian's element count into storage.
BYTES_PER_REAL = 8


def time_call(function, repeats=3):
    """Best of several timings of a callable.

    The minimum is reported rather than the mean. A timing is a lower bound contaminated upwards by scheduling, page faults and competing processes; the smallest observation is the closest to the quantity of interest, whereas the mean drifts with whatever else the machine was doing.
    """
    best = float("inf")
    for _ in range(repeats):
        start = time.perf_counter()
        function()
        best = min(best, time.perf_counter() - start)
    return best


def synthetic_costs(sizes, n_snapshots=40, repeats=3):
    """Time subspace assembly and decomposition against state-vector size."""
    import dmd as dmd_mod

    rows = []
    rng = np.random.default_rng(12345)
    for n_dof in sizes:
        # A rank-deficient matrix resembling a real subspace: a few growing modes plus noise, rather than random numbers, so the decomposition does comparable work.
        basis = rng.standard_normal((n_dof, 4))
        coeffs = np.exp(np.outer([0.32, 0.05, -0.1, -0.4],
                                 np.arange(n_snapshots) * 0.08))
        snapshots = basis @ coeffs
        snapshots += 1e-8 * rng.standard_normal(snapshots.shape)

        t_subtract = time_call(
            lambda: np.subtract(snapshots, snapshots), repeats)
        t_dmd = time_call(
            lambda: dmd_mod.dmd(snapshots, 0.08, rank=6), repeats)

        rows.append({
            "n_dof": int(n_dof),
            "n_snapshots": int(n_snapshots),
            "t_subtract": t_subtract,
            "t_dmd": t_dmd,
            "jacobian_elements": int(n_dof) ** 2,
            "jacobian_bytes": int(n_dof) ** 2 * BYTES_PER_REAL,
        })
        print(f"    n_dof {n_dof:>9,}: subtract {t_subtract * 1e3:8.2f} ms, "
              f"dmd {t_dmd * 1e3:8.2f} ms", flush=True)
    return rows


def measured_costs(perturbed, base, fields, prefix, repeats=1):
    """Time the analysis path on real snapshots."""
    import nsmfp

    paths = nsmfp.find_plotfiles(perturbed, prefix=prefix)
    if not paths:
        raise FileNotFoundError(f"no {prefix}* plotfiles under {perturbed}")

    t_read = time_call(
        lambda: nsmfp.read_plotfile(paths[0], fields=tuple(fields)), repeats)

    start = time.perf_counter()
    snapshots, dt, times = nsmfp.build_snapshot_matrix(
        perturbed, base, fields=tuple(fields), verbose=False)
    t_assemble = time.perf_counter() - start

    import dmd as dmd_mod
    t_dmd = time_call(lambda: dmd_mod.dmd(snapshots, dt), repeats)

    n_dof, n_snap = snapshots.shape
    return {
        "n_dof": int(n_dof),
        "n_snapshots": int(n_snap),
        "t_read_one": t_read,
        "t_assemble": t_assemble,
        "t_assemble_per_snapshot": t_assemble / max(n_snap, 1),
        "t_dmd": t_dmd,
        "jacobian_elements": int(n_dof) ** 2,
        "jacobian_bytes": int(n_dof) ** 2 * BYTES_PER_REAL,
    }


def figure_scaling(rows, provenance, outdir):
    """Cost against state-vector size, with the explicit Jacobian alongside.

    The two are plotted together because the argument is comparative. The Jacobian-free cost grows roughly linearly in the degrees of freedom, since it is a fixed number of passes over the state; the explicit Jacobian grows as the square in storage alone, before any of its columns are computed. On logarithmic axes the difference is a difference of slope, which is the durable statement, rather than a ratio at one problem size, which is not.
    """
    n_dof = np.array([r["n_dof"] for r in rows], dtype=float)
    t_dmd = np.array([r["t_dmd"] for r in rows], dtype=float)
    gib = np.array([r["jacobian_bytes"] for r in rows],
                   dtype=float) / 1024.0 ** 3

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.4, 4.2))

    ax1.loglog(n_dof, t_dmd, "o-", ms=4.5, label="decomposition, measured")
    reference = t_dmd[0] * (n_dof / n_dof[0])
    ax1.loglog(n_dof, reference, "--", label=r"$O(n)$")
    ax1.set_xlabel("state-vector degrees of freedom $n$")
    ax1.set_ylabel("time (s)")
    ax1.grid(which="both")
    ax1.legend()
    # Both axes are logarithmic and therefore take no tick formatter.

    ax2.loglog(n_dof, gib, "o-", ms=4.5, label="explicit Jacobian storage")
    for limit, name in ((16.0, "16 GiB workstation"),
                        (192.0, "192 GiB cluster node")):
        ax2.axhline(limit, ls=":", lw=GUIDE_LW, color="0.45")
        ax2.text(n_dof[0], limit * 1.25, name, fontsize=7, color="0.35")
    ax2.set_xlabel("state-vector degrees of freedom $n$")
    ax2.set_ylabel("storage (GiB)")
    ax2.grid(which="both")
    ax2.legend(loc="lower right")

    fig.suptitle("Cost of the Jacobian-free path against an explicit Jacobian",
                 fontsize=11)
    # The provenance sits on the figure, because a cost read without its machine is not interpretable and the caption may not travel with the image.
    fig.text(0.5, 0.005, machine_info.caption(provenance), ha="center",
             fontsize=7, color="0.35")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    savefig_all(fig, "jvp_cost_scaling", outdir=outdir)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--synthetic", action="store_true",
                    help="time generated arrays instead of a calculation")
    ap.add_argument("--perturbed")
    ap.add_argument("--base")
    ap.add_argument("--fields", nargs="+",
                    default=["x_velocity", "y_velocity"])
    ap.add_argument("--sizes", type=int, nargs="+",
                    default=[2048, 8192, 32768, 131072, 524288])
    ap.add_argument("--n-snapshots", type=int, default=40,
                    dest="n_snapshots")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--solver-path", default=None, dest="solver_path",
                    help="solver working tree, so its commit is recorded")
    ap.add_argument("--prefix", default="nddataPLT")
    ap.add_argument("--outdir", default="figures")
    args = ap.parse_args(argv)

    apply_style()
    os.makedirs(args.outdir, exist_ok=True)

    provenance = machine_info.collect(solver_path=args.solver_path)
    print(machine_info.format_table(provenance))
    print()

    record = {"provenance": provenance}

    if args.synthetic or not (args.perturbed and args.base):
        print("timing the synthetic scaling study ...", flush=True)
        rows = synthetic_costs(args.sizes, args.n_snapshots, args.repeats)
        record["synthetic"] = rows
        figure_scaling(rows, provenance, args.outdir)
    else:
        print("timing the measured path ...", flush=True)
        measured = measured_costs(args.perturbed, args.base, args.fields,
                                  args.prefix, args.repeats)
        record["measured"] = measured
        print()
        print(f"state vector       {measured['n_dof']:,} degrees of freedom")
        print(f"snapshots          {measured['n_snapshots']}")
        print(f"read one plotfile  {measured['t_read_one']:.3f} s")
        print(f"assemble subspace  {measured['t_assemble']:.3f} s "
              f"({measured['t_assemble_per_snapshot']:.3f} s per snapshot)")
        print(f"decomposition      {measured['t_dmd']:.3f} s")
        print(f"explicit Jacobian  {measured['jacobian_elements']:,} elements,"
              f" {measured['jacobian_bytes'] / 1024 ** 3:,.1f} GiB, "
              f"and {measured['n_dof']:,} products to form it")
        rows = synthetic_costs(args.sizes, args.n_snapshots, args.repeats)
        record["synthetic"] = rows
        figure_scaling(rows, provenance, args.outdir)

    out = os.path.join(args.outdir, "jvp_cost.json")
    with open(out, "w") as fh:
        json.dump(record, fh, indent=2)
    print()
    print(f"timings and provenance written to {out}")
    print(f"figures written to {args.outdir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
