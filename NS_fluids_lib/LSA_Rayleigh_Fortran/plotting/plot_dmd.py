# !/usr/bin/env python3
"""Figures for the modal extraction: spectrum, amplitudes and convergence.

The growth rate reported by this study is obtained from a decomposition of a snapshot subspace, yet the decomposition itself appears nowhere in the figures. A reader is asked to accept a number without seeing the spectrum it came from, and cannot tell whether the physical mode was cleanly separated from the numerical ones or merely happened to be the largest. These three figures supply that evidence.

The spectrum places the discrete-time eigenvalues in the complex plane against the unit circle, which is where a growth rate becomes legible geometrically: an eigenvalue outside the circle grows, one inside decays, and one on it is neutral. The capillary mode should sit on the positive real axis just outside the circle, since it grows without oscillating, and the numerical modes should sit at or within the circle. That separation is the claim the figure makes, and it is the one a reader will check.

The amplitude ranking shows how much of the first snapshot each mode carries. It is what justifies attending to one mode and discarding the rest, and it is more honest than asserting a spectral gap, since the gap is visible or it is not.

The convergence figure varies the size and the identity of the subspace, following the procedure of Ranjan and coauthors, and plots the recovered rate against the number of snapshots retained. A rate that is a property of the flow is flat across that variation; a rate that is an artefact of the window is not. This is what distinguishes modal extraction from curve fitting, and for a study whose contribution is methodological it is arguably the single most important figure of the set.

Snapshots are read from a twin pair, so this script needs the plotfiles and needs yt. It is therefore run on the machine that holds them, and it writes the numbers behind each figure to a `.npz` so that the figures can afterwards be redrawn without repeating the decomposition.

Run from `amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Fortran`:

    python3 plotting/plot_dmd.py --perturbed runs/k0.7_eps0.02 --base runs/k0.7_base --outdir figures
    python3 plotting/plot_dmd.py --replot --outdir figures
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import dispersion
from pub_style import (apply_style, savefig_all, cm_both, cm_x,
                       save_data, load_data, has_data)

GUIDE_LW = 0.8

#: Stem under which the numbers behind these figures are stored.
DATA = "dmd_figures"


def figure_spectrum(ritz, amplitude, dt, sigma_expected, outdir):
    """Discrete-time eigenvalues against the unit circle.

    Marker area follows amplitude, so the mode carrying the flow is visually dominant without being singled out by hand. The unit circle is drawn because it is the neutral locus: the radial distance of an eigenvalue from it is the growth over one sampling interval, and a reader can therefore read stability off the figure directly.
    """
    ritz = np.asarray(ritz)
    amplitude = np.asarray(amplitude, dtype=float)

    fig, ax = plt.subplots(figsize=(5.4, 5.0))
    theta = np.linspace(0.0, 2.0 * np.pi, 512)
    ax.plot(np.cos(theta), np.sin(theta), "-", lw=GUIDE_LW, color="0.45",
            label="unit circle (neutral)")

    scale = amplitude / amplitude.max() if amplitude.max() > 0 else amplitude
    ax.scatter(ritz.real, ritz.imag, s=12.0 + 90.0 * scale,
               facecolors="none", edgecolors="k", linewidths=0.8,
               label="Ritz values")

    # The mode of largest amplitude is the one the growth rate is taken from, so it is named rather than left for the reader to infer from marker size alone.
    lead = int(np.argmax(amplitude))
    ax.plot(ritz[lead].real, ritz[lead].imag, "o", ms=7, mfc="none",
            label="leading mode")

    if sigma_expected is not None:
        ax.plot([np.exp(sigma_expected * dt)], [0.0], "x", ms=7,
                label=rf"predicted, $\sigma={sigma_expected:.4f}$")

    ax.set_xlabel(r"$\mathrm{Re}\,\lambda$")
    ax.set_ylabel(r"$\mathrm{Im}\,\lambda$")
    ax.set_aspect("equal")
    span = max(1.15, float(np.abs(ritz).max()) * 1.08)
    ax.set_xlim(-span, span)
    ax.set_ylim(-span, span)
    cm_both(ax)
    ax.legend(loc="lower left")
    fig.tight_layout()
    savefig_all(fig, "dmd_spectrum", outdir=outdir)
    plt.close(fig)
    return lead


def figure_amplitudes(growth, amplitude, sigma_expected, outdir):
    """Mode amplitude against growth rate.

    Plotting amplitude against growth rate rather than against mode index puts the two quantities that matter on the same axes: a mode is worth attending to if it both grows and carries weight, and this figure shows at a glance whether any mode other than the physical one does both.
    """
    growth = np.asarray(growth, dtype=float)
    amplitude = np.asarray(amplitude, dtype=float)
    positive = amplitude > 0

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.semilogy(growth[positive], amplitude[positive], "o", ms=4.5,
                mfc="none", color="k")
    ax.axvline(0.0, lw=GUIDE_LW, color="0.45")
    if sigma_expected is not None:
        ax.axvline(sigma_expected, ls="--", lw=GUIDE_LW, color="0.35",
                   label=rf"predicted $\sigma={sigma_expected:.4f}$")
        ax.legend(loc="lower left")
    ax.set_xlabel(r"$\sigma = \mathrm{Re}\,(\log\lambda)/\Delta t$")
    ax.set_ylabel("mode amplitude")
    # The ordinate is logarithmic and therefore takes no tick formatter.
    cm_x(ax)
    fig.tight_layout()
    savefig_all(fig, "dmd_amplitudes", outdir=outdir)
    plt.close(fig)


def figure_convergence(sizes, rates, starts, sigma_expected, outdir):
    """Recovered growth rate against subspace size, for several start points.

    Varying the identity of the subspace as well as its size is what separates a converged rate from one that merely stopped changing. Curves from different starting snapshots collapsing onto the same value is the evidence; a single curve flattening is not, since it is consistent with a window-dependent artefact that happens to be stable.
    """
    sizes = np.asarray(sizes)
    rates = np.asarray(rates, dtype=float)

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for index, start in enumerate(starts):
        finite = np.isfinite(rates[index])
        if not finite.any():
            continue
        ax.plot(sizes[finite], rates[index][finite], "o-", ms=4.0,
                label=rf"from snapshot ${start}$")

    if sigma_expected is not None:
        ax.axhline(sigma_expected, ls="--", lw=GUIDE_LW, color="0.35",
                   label=rf"predicted $\sigma={sigma_expected:.4f}$")

    ax.set_xlabel("snapshots retained in the subspace")
    ax.set_ylabel(r"recovered $\sigma$")
    ax.set_xlim(float(sizes.min()), float(sizes.max()))
    cm_both(ax)
    ax.legend()
    fig.tight_layout()
    savefig_all(fig, "dmd_convergence", outdir=outdir)
    plt.close(fig)


def figure_singular_values(s, rank, outdir):
    """Singular values of the snapshot matrix, with the truncation marked.

    The truncation rank is a choice, and a figure that shows where it fell relative to the decay of the spectrum is what turns that choice from arbitrary into justified.
    """
    s = np.asarray(s, dtype=float)
    normalised = s / s[0] if s[0] > 0 else s

    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    ax.semilogy(np.arange(1, s.size + 1), normalised, "o-", ms=4.0)
    ax.axvline(rank, ls="--", lw=GUIDE_LW, color="0.35",
               label=f"truncation, rank {rank}")
    ax.set_xlabel("index")
    ax.set_ylabel(r"$s_i / s_1$")
    ax.set_xlim(1, s.size)
    cm_x(ax)
    ax.legend()
    fig.tight_layout()
    savefig_all(fig, "dmd_singular_values", outdir=outdir)
    plt.close(fig)


def compute(args):
    """Assemble the subspace, decompose it, and store what the figures need."""
    import dmd as dmd_mod
    import nsmfp

    print("assembling the snapshot matrix ...", flush=True)
    # build_snapshot_matrix returns the sampling interval alongside the matrix, having already verified that the snapshots are uniformly spaced; taking a difference of the first two times here would silently accept a non-uniform record that it would have rejected.
    snapshots, dt, times = nsmfp.build_snapshot_matrix(
        args.perturbed, args.base, fields=tuple(args.fields),
        skip=args.skip, stride=args.stride)
    dt = float(dt)
    print(f"    {snapshots.shape[1]} snapshots, dt = {dt:.6g}", flush=True)

    result = dmd_mod.dmd(snapshots, dt, rank=args.rank)
    print(f"    rank {result.rank}, leading sigma "
          f"{result.leading_growth_rate():.6f}", flush=True)

    sizes = [n for n in np.unique(np.linspace(
        8, snapshots.shape[1], 10).astype(int)) if n > 3]
    starts = [int(v) for v in args.starts]
    rates = np.full((len(starts), len(sizes)), np.nan)
    for i, start in enumerate(starts):
        for j, n in enumerate(sizes):
            if start + n > snapshots.shape[1]:
                continue
            sub = dmd_mod.dmd(snapshots[:, start:start + n], dt,
                              rank=args.rank)
            rates[i, j] = sub.leading_growth_rate()

    save_data(DATA,
              ritz_real=result.ritz.real, ritz_imag=result.ritz.imag,
              amplitude=result.amplitude, growth=result.growth_rate,
              singular_values=result.singular_values,
              sizes=np.asarray(sizes), rates=rates,
              starts=np.asarray(starts),
              meta={"dt": dt, "rank": int(result.rank),
                    "n_snapshots": int(snapshots.shape[1]),
                    "kr0": args.kr0})
    return DATA


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--perturbed", help="run directory of the perturbed case")
    ap.add_argument("--base", help="run directory of the unperturbed twin")
    ap.add_argument("--replot", action="store_true",
                    help="redraw from stored data, without decomposing again")
    ap.add_argument("--fields", nargs="+",
                    default=["x_velocity", "y_velocity"])
    ap.add_argument("--skip", type=int, default=60,
                    help="snapshots discarded as startup transient")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--rank", type=int, default=None)
    ap.add_argument("--starts", type=int, nargs="+", default=[0, 5, 10])
    ap.add_argument("--kr0", type=float, default=0.7)
    ap.add_argument("--r0", type=float, default=1.0)
    ap.add_argument("--mu", type=float, default=0.02)
    ap.add_argument("--outdir", default="figures")
    args = ap.parse_args(argv)

    apply_style()
    os.makedirs(args.outdir, exist_ok=True)

    if not args.replot:
        if not (args.perturbed and args.base):
            ap.error("--perturbed and --base are required unless --replot")
        compute(args)
    elif not has_data(DATA):
        ap.error(f"no stored data for {DATA}; run once without --replot")

    d = load_data(DATA)
    dt = float(d["meta"]["dt"])
    kr0 = float(d["meta"].get("kr0", args.kr0))
    sigma_expected = float(dispersion.growth_rate_viscous(
        kr0 / args.r0, r0=args.r0, mu=args.mu))

    ritz = d["ritz_real"] + 1j * d["ritz_imag"]
    lead = figure_spectrum(ritz, d["amplitude"], dt, sigma_expected,
                           args.outdir)
    figure_amplitudes(d["growth"], d["amplitude"], sigma_expected,
                      args.outdir)
    figure_convergence(d["sizes"], d["rates"], list(d["starts"]),
                       sigma_expected, args.outdir)
    figure_singular_values(d["singular_values"], int(d["meta"]["rank"]),
                           args.outdir)

    print(f"leading mode     lambda = {ritz[lead]:.6f}")
    print(f"                 sigma  = {d['growth'][lead]:.6f}")
    print(f"predicted        sigma  = {sigma_expected:.6f}")
    print(f"figures written to {args.outdir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
