# !/usr/bin/env python3
"""Growth rate over the wavenumber and Ohnesorge plane.

The dispersion figure reports the growth rate along one line of the parameter space, at the single Ohnesorge number of the reference calculation. That is enough to validate a measurement but it says nothing about how the instability behaves as viscosity is varied, and a reader is entitled to ask whether the agreement found at Oh = 0.02 is characteristic or particular. This script extends the curve to a map.

Two panels are drawn. The surface shows the shape of the growth rate over the plane, and is the panel a reader looks at first because the ridge and its descent are immediately legible. The contour map is the panel a reader measures from, since a value can be read off it and the locus of the most unstable wavenumber can be drawn on it as a line. Neither replaces the other, which is why both are produced.

Two features carry the argument. The growth rate falls monotonically with Ohnesorge number at every wavenumber, since viscosity can only retard a capillary instability, and any measured point lying above the inviscid curve is therefore an error rather than a discovery. And the most unstable wavenumber drifts to longer waves as viscosity rises, departing from the inviscid value of 0.697, which is the quantitative statement that the reference measurement at $kr_0 = 0.7$ sits near but not exactly at the peak.

Everything here is evaluated from the dispersion relation and no calculation is read, so the figure is a prediction against which measurements are placed rather than a result. Measured points may be overlaid from a sweep with `--records`, in which case they are drawn on the contour panel at their own Ohnesorge number.

Run from `amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Python`:

    python3 python/plot_ohnesorge.py --outdir figures
    python3 python/plot_ohnesorge.py --oh-max 0.5 --n-k 240 --n-oh 180
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
                       surface, zaxis_left, subcaption, SURFACE_CMAP)

GUIDE_LW = 0.8


def sweep(kr0, oh, r0=1.0, sigma=1.0, rho=1.0):
    """Growth rate on the (kr0, Oh) grid.

    The Ohnesorge number fixes the viscosity through Oh = mu / sqrt(rho sigma r0), so each row of the grid is one viscosity and the existing viscous dispersion relation is evaluated along it. Building the map from the same function the one-dimensional figure uses is deliberate: a separate implementation for the surface could disagree with the curve and the disagreement would be invisible.
    """
    growth = np.empty((oh.size, kr0.size))
    for row, oh_value in enumerate(oh):
        mu = oh_value * np.sqrt(rho * sigma * r0)
        for col, k in enumerate(kr0):
            growth[row, col] = dispersion.growth_rate_viscous(
                k / r0, r0=r0, sigma=sigma, rho=rho, mu=mu)
    # Only the unstable band is of interest; above the cutoff the mode oscillates rather than growing, and a negative or zero rate plotted as height would be read as decay at a rate the figure does not represent.
    return np.where(growth > 0.0, growth, np.nan)


def peak_locus(kr0, oh, growth):
    """Most unstable wavenumber at each Ohnesorge number.

    The location is refined below the grid spacing by fitting a parabola through the largest sample and its two neighbours. Taking the grid maximum alone quantises the locus onto the abscissa, and the resulting staircase reads as noise in a quantity that is in fact smooth, inviting a question about numerical scatter where none exists.
    """
    peaks = np.full(oh.size, np.nan)
    for row in range(oh.size):
        line = growth[row]
        if np.all(np.isnan(line)):
            continue
        best = int(np.nanargmax(line))
        if 0 < best < kr0.size - 1 and np.all(np.isfinite(
                line[best - 1:best + 2])):
            left, mid, right = line[best - 1:best + 2]
            denom = left - 2.0 * mid + right
            # A vanishing denominator means the three samples are collinear, so there is no interior maximum to refine towards and the grid point stands.
            shift = 0.5 * (left - right) / denom if denom != 0.0 else 0.0
            shift = float(np.clip(shift, -0.5, 0.5))
            peaks[row] = kr0[best] + shift * (kr0[1] - kr0[0])
        else:
            peaks[row] = kr0[best]
    return peaks


def figure_surface(kr0, oh, growth, outdir):
    """The growth rate as a surface over the plane."""
    fig = plt.figure(figsize=(6.4, 4.8))
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    zaxis_left(ax)
    K, O = np.meshgrid(kr0, oh)
    # NaN above the cutoff would tear the surface, so the stable region is flattened to zero height and the caption records that the surface is clipped there rather than falling to zero physically.
    Z = np.nan_to_num(growth, nan=0.0)
    surface(fig, ax, K, O, Z, zlabel=r"$\sigma$", cmap=SURFACE_CMAP)
    ax.set_xlabel(r"$k r_0$")
    ax.set_ylabel(r"$\mathrm{Oh}$")
    subcaption(ax, "(a) growth rate over the parameter plane")
    fig.tight_layout()
    savefig_all(fig, "ohnesorge_surface", outdir=outdir)
    plt.close(fig)


def figure_contour(kr0, oh, growth, measured, outdir):
    """The growth rate as a contour map, with the locus of the peak."""
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    levels = np.linspace(0.0, float(np.nanmax(growth)), 21)
    filled = ax.contourf(kr0, oh, growth, levels=levels, cmap=SURFACE_CMAP)
    lines = ax.contour(kr0, oh, growth, levels=levels[::4], colors="k",
                       linewidths=GUIDE_LW)
    ax.clabel(lines, fmt="%.2f", fontsize=7)

    peaks = peak_locus(kr0, oh, growth)
    ax.plot(peaks, oh, "--", color="w", label="most unstable wavenumber")
    ax.axvline(dispersion.most_unstable_wavenumber()[0], ls=":",
               lw=GUIDE_LW, color="w")

    if measured is not None and len(measured):
        m = np.asarray(measured, dtype=float)
        ax.plot(m[:, 0], m[:, 1], "o", ms=4.5, mfc="none", color="k",
                label="measured")

    ax.set_xlabel(r"$k r_0$")
    ax.set_ylabel(r"$\mathrm{Oh}$")
    ax.set_xlim(float(kr0.min()), float(kr0.max()))
    ax.set_ylim(float(oh.min()), float(oh.max()))
    # The filled contours would otherwise carry the global dashed grid across them.
    ax.grid(False)
    cm_both(ax)
    ax.legend(loc="upper left", framealpha=0.85)
    fig.colorbar(filled, ax=ax, pad=0.02, label=r"$\sigma$")
    subcaption(ax, "(b) contours of the growth rate")
    fig.tight_layout()
    savefig_all(fig, "ohnesorge_contour", outdir=outdir)
    plt.close(fig)


def figure_peak_drift(oh, peaks, outdir):
    """Drift of the most unstable wavenumber with viscosity.

    Read from the map this is a curve on a surface; drawn alone it is a quantitative statement, and it is what licenses the claim that the reference wavenumber sits near the peak rather than at it.
    """
    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    ax.plot(oh, peaks, "-")
    ax.axhline(dispersion.most_unstable_wavenumber()[0], ls=":",
               lw=GUIDE_LW, color="0.4", label="inviscid peak, $0.697$")
    ax.set_xlabel(r"$\mathrm{Oh}$")
    ax.set_ylabel(r"$(k r_0)_{\mathrm{peak}}$")
    ax.set_xlim(float(oh.min()), float(oh.max()))
    cm_both(ax)
    ax.legend()
    fig.tight_layout()
    savefig_all(fig, "ohnesorge_peak_drift", outdir=outdir)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kr0-min", type=float, default=0.02)
    ap.add_argument("--kr0-max", type=float, default=1.0)
    ap.add_argument("--oh-min", type=float, default=0.0)
    ap.add_argument("--oh-max", type=float, default=0.4)
    ap.add_argument("--n-k", type=int, default=161, dest="n_k")
    ap.add_argument("--n-oh", type=int, default=121, dest="n_oh")
    ap.add_argument("--r0", type=float, default=1.0)
    ap.add_argument("--sigma", type=float, default=1.0)
    ap.add_argument("--rho", type=float, default=1.0)
    ap.add_argument("--measured", type=float, nargs="+", default=None,
                    help="pairs of kr0 and Oh to overlay on the contour map")
    ap.add_argument("--outdir", default="figures")
    args = ap.parse_args(argv)

    apply_style()
    os.makedirs(args.outdir, exist_ok=True)

    kr0 = np.linspace(args.kr0_min, args.kr0_max, args.n_k)
    oh = np.linspace(args.oh_min, args.oh_max, args.n_oh)
    growth = sweep(kr0, oh, r0=args.r0, sigma=args.sigma, rho=args.rho)

    measured = None
    if args.measured:
        if len(args.measured) % 2:
            ap.error("--measured takes pairs of kr0 and Oh")
        measured = np.asarray(args.measured, dtype=float).reshape(-1, 2)

    figure_surface(kr0, oh, growth, args.outdir)
    figure_contour(kr0, oh, growth, measured, args.outdir)
    peaks = peak_locus(kr0, oh, growth)
    figure_peak_drift(oh, peaks, args.outdir)

    finite = peaks[np.isfinite(peaks)]
    print(f"grid             {args.n_k} wavenumbers by {args.n_oh} "
          f"Ohnesorge numbers")
    print(f"largest sigma    {np.nanmax(growth):.5f} at Oh = {oh.min():.3f}")
    print(f"peak wavenumber  {finite[0]:.4f} at Oh = {oh[0]:.3f}"
          f"  ->  {finite[-1]:.4f} at Oh = {oh[-1]:.3f}")
    print(f"inviscid peak    "
          f"{dispersion.most_unstable_wavenumber()[0]:.4f}")
    print(f"figures written to {args.outdir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
