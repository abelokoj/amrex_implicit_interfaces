# !/usr/bin/env python3
"""Schematic of the twin-run architecture.

The central methodological decision of this work is described in prose and drawn nowhere. NS-MFP as originally formulated adds a body force B_f = dQ/dt - F(Q) to hold the base state stationary, which requires modifying the solver. Here that body force is replaced by a second calculation: an unperturbed twin is advanced from the same initial condition under an identical discretisation, and the two solutions are differenced snapshot by snapshot, so the drift cancels to machine precision and the solver is used exactly as distributed.

That substitution is easier to see than to read. One diagram carries it: two parallel calculations, a difference taken at matched times, and a decomposition applied to the resulting sequence. For a defence in particular it does more work than any three quantitative figures, because it lets an audience follow the rest of the argument without holding the architecture in memory.

The figure is drawn rather than data-derived, so it carries no measurement and makes no claim that could be checked against a calculation. It uses the house palette and typography so that it does not read as an intruder among the measured figures.

Run from `amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Fortran`:

    python3 plotting/plot_schematic.py --outdir figures
"""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pub_style import apply_style, savefig_all, PALETTE

#: Muted fills, so that the boxes read as background and the text as content. Saturated palette entries are reserved for the arrows, which carry the flow of the diagram.
FILL_BASE = "#eaf4fd"
FILL_PERT = "#fdeaea"
FILL_WORK = "#f0f0f0"
FILL_OUT = "#eaf6ec"


def box(ax, x, y, w, h, text, fill, fontsize=8.5, weight="normal"):
    """A rounded box with centred text."""
    patch = FancyBboxPatch((x, y), w, h,
                           boxstyle="round,pad=0.012,rounding_size=0.02",
                           linewidth=0.8, edgecolor="0.35", facecolor=fill,
                           zorder=2)
    ax.add_patch(patch)
    ax.text(x + w / 2.0, y + h / 2.0, text, ha="center", va="center",
            fontsize=fontsize, zorder=3, weight=weight, linespacing=1.5)
    return (x + w / 2.0, y + h / 2.0)


def arrow(ax, start, end, colour="0.3", style="-|>", lw=1.0, text=None,
          offset=(0.0, 0.02), fontsize=7.5, dashed=False):
    """A connector, optionally labelled."""
    patch = FancyArrowPatch(start, end, arrowstyle=style,
                            mutation_scale=9, linewidth=lw, color=colour,
                            zorder=1, shrinkA=2, shrinkB=2,
                            linestyle="--" if dashed else "-")
    ax.add_patch(patch)
    if text:
        mx = 0.5 * (start[0] + end[0]) + offset[0]
        my = 0.5 * (start[1] + end[1]) + offset[1]
        ax.text(mx, my, text, ha="center", va="bottom", fontsize=fontsize,
                color=colour, zorder=3)


def draw(ax):
    """Lay out the diagram on a unit canvas."""
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis("off")

    blue, red = PALETTE[0], PALETTE[1]

    # --- initial condition -------------------------------------------
    ic = box(ax, 0.02, 0.60, 0.155, 0.16,
             "quiescent\ncolumn\n$r = r_0$", FILL_WORK)

    # --- the two calculations ----------------------------------------
    base = box(ax, 0.245, 0.755, 0.30, 0.155,
               "unperturbed twin\n"
               r"$\bar{Q}(t)$, solver unmodified", FILL_BASE)
    pert = box(ax, 0.245, 0.455, 0.30, 0.155,
               "perturbed calculation\n"
               r"$Q(t)$, $r = r_0 + \varepsilon\cos kz$", FILL_PERT)

    arrow(ax, (ic[0] + 0.08, ic[1] + 0.035), (0.245, base[1]), colour=blue)
    arrow(ax, (ic[0] + 0.08, ic[1] - 0.035), (0.245, pert[1]), colour=red)

    # --- the difference ----------------------------------------------
    diff = box(ax, 0.615, 0.605, 0.20, 0.155,
               "difference at\nmatched times\n"
               r"$Q' = Q - \bar{Q}$", FILL_WORK, weight="bold")
    arrow(ax, (0.545, base[1]), (0.615, diff[1] + 0.045), colour=blue)
    arrow(ax, (0.545, pert[1]), (0.615, diff[1] - 0.045), colour=red)

    # --- the outputs --------------------------------------------------
    snap = box(ax, 0.615, 0.335, 0.20, 0.13,
               "snapshot\nsubspace\n" r"$\{Q'(t_n)\}$", FILL_WORK)
    arrow(ax, (diff[0], 0.605), (snap[0], 0.465))

    dmd = box(ax, 0.615, 0.10, 0.20, 0.14,
              "modal extraction\n" r"$\lambda,\ \ \sigma = \mathrm{Re}\,"
              r"\log\lambda\,/\,\Delta t$", FILL_OUT)
    arrow(ax, (snap[0], 0.335), (dmd[0], 0.24))

    inter = box(ax, 0.245, 0.10, 0.28, 0.14,
                "interface mode\n" r"$a(t)$ from the zero contour",
                FILL_OUT)
    arrow(ax, (pert[0] - 0.04, 0.455), (inter[0], 0.24), colour=red,
          dashed=True)

    # --- what is replaced ---------------------------------------------
    # Stated as a struck-through alternative rather than omitted, because the contribution is the substitution and a diagram showing only the adopted path would not convey what was avoided.
    ax.text(0.02, 0.30,
            "replaces the constraining body force\n"
            r"$B_f = \partial\bar{Q}/\partial t - F(\bar{Q})$,"
            "\nwhich would require modifying\nthe solver",
            fontsize=7.5, va="top", ha="left", color="0.35",
            linespacing=1.6)
    # The dashed path is drawn because it is a result in itself: any study confined to growth rates needs only the perturbed calculation, which halves its cost.
    ax.text(0.385, 0.055,
            "the interface observable needs no twin,\n"
            "which halves the cost of a growth-rate study",
            fontsize=7, va="top", ha="center", color=red, alpha=0.9,
            linespacing=1.5)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="figures")
    args = ap.parse_args(argv)

    apply_style()
    os.makedirs(args.outdir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    # The global dashed grid would run across a diagram that has no axes.
    ax.grid(False)
    draw(ax)
    fig.tight_layout()
    savefig_all(fig, "schematic_twin_run", outdir=args.outdir)
    plt.close(fig)

    print(f"schematic written to {args.outdir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
