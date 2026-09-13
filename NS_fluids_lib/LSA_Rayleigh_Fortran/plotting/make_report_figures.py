# !/usr/bin/env python3
"""Generate the figures reported in ``analysis_report.typ``.

The script reads the plain-text amplitude records written by the Fortran tool ``bin/dump_amplitude`` together with the archived record from the earlier solver build, so it requires neither yt nor a live plotfile tree. All typographic settings come from ``pub_style``; every figure is written as PNG, PDF and SVG at 600 dpi.

Run from ``amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Python``:

    python3 python/make_report_figures.py --datadir runs/records \\
        --archive validation_data/validation_kr0_0.7.json --outdir figures
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import dispersion
from pub_style import apply_style, savefig_all, cm_both, cm_x, cm_y

#: Width of a reference or guide line, matched to the axis frame so that a guide never competes with a data curve. Data curves take lines.linewidth from pub_style, which is 1.0, and therefore carry no explicit width anywhere in this file.
GUIDE_LW = 0.8


# Growth rates predicted for the configuration used throughout: r0 = sigma = rho = 1 and mu = 0.02, giving an Ohnesorge number of 0.02.
KR0 = 0.7
SIGMA_INVISCID = float(dispersion.growth_rate(KR0, 1.0, 1.0, 1.0))
SIGMA_VISCOUS = float(dispersion.growth_rate_viscous(KR0, 1.0, 1.0, 1.0, 0.02))


def local_sigma(t, a, width=0.8, n_window=8):
    """Growth rate over sliding windows of the given width.

    Returns empty arrays when the record is shorter than one window, which keeps the caller free of special cases for very short runs.
    """
    lo, hi = t.min() + 0.5 * width, t.max() - 0.5 * width
    if hi <= lo:
        return np.array([]), np.array([])
    centres, rates = [], []
    for centre in np.linspace(lo, hi, n_window):
        mask = (t >= centre - 0.5 * width) & (t <= centre + 0.5 * width) & (a > 0)
        if mask.sum() >= 4:
            centres.append(centre)
            rates.append(np.polyfit(t[mask], np.log(a[mask]), 1)[0])
    return np.asarray(centres), np.asarray(rates)


def fit(t, a, t_min, t_max=None):
    """Least-squares slope of log a(t) with the coefficient of determination."""
    mask = (t >= t_min) & (a > 0)
    if t_max is not None:
        mask &= t <= t_max
    slope, intercept = np.polyfit(t[mask], np.log(a[mask]), 1)
    resid = np.log(a[mask]) - (slope * t[mask] + intercept)
    ss_tot = np.sum((np.log(a[mask]) - np.log(a[mask]).mean()) ** 2)
    return float(slope), 1.0 - float(np.sum(resid**2)) / float(ss_tot)


def figure_validation(t, a, outdir):
    """Amplitude and local growth rate from the long reference calculation."""
    t_fit = 4.5
    slope, r2 = fit(t, a, t_fit)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.0))

    ax1.semilogy(t, a, "o", ms=3.0, mfc="none", color="k", label="measured")
    mask = t >= t_fit
    tw = t[mask]
    anchor = a[mask][0]
    ax1.semilogy(tw, anchor * np.exp(slope * (tw - tw[0])), "-",
                 label=rf"fit, $\sigma={slope:.4f}$")
    ax1.semilogy(tw, anchor * np.exp(SIGMA_VISCOUS * (tw - tw[0])), "--",
                 label=rf"viscous, $\sigma={SIGMA_VISCOUS:.4f}$")
    ax1.semilogy(tw, anchor * np.exp(SIGMA_INVISCID * (tw - tw[0])), ":",
                 label=rf"inviscid, $\sigma={SIGMA_INVISCID:.4f}$")
    ax1.axvspan(t.min(), t_fit, color="0.88", zorder=0, lw=0)
    ax1.set_xlabel("$t$")
    ax1.set_ylabel("$a(t)$")
    ax1.set_xlim(t.min(), t.max())
    # Logarithmic ordinate: the Computer Modern tick formatter is applied to the abscissa only, since it would replace the powers of ten.
    cm_x(ax1)
    ax1.legend(loc="lower right", fontsize=8)

    centres, rates = local_sigma(t, a)
    ax2.plot(centres, rates, "o-", ms=4.0, label=r"local $\sigma(t)$")
    ax2.axhline(SIGMA_VISCOUS, ls="--", color="0.35",
                label="viscous estimate")
    ax2.axhline(SIGMA_INVISCID, ls=":", color="0.55",
                label="inviscid")
    ax2.set_xlabel("$t$")
    ax2.set_ylabel(r"$\mathrm{d}(\ln a)/\mathrm{d}t$")
    ax2.set_xlim(float(centres.min()), float(centres.max()))
    ax2.set_ylim(0.0, SIGMA_INVISCID * 1.12)
    cm_both(ax2)
    ax2.legend(loc="lower right", fontsize=8)

    fig.tight_layout()
    savefig_all(fig, "fig1_validation", outdir=outdir)
    plt.close(fig)
    return slope, r2


def figure_replication(t_new, a_new, t_old, a_old, outdir):
    """Comparison of the two solver builds on the same configuration."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.0))

    ax1.semilogy(t_old, a_old, "-", color="0.6",
                 label="cbcee7c, locally corrected")
    ax1.semilogy(t_new, a_new, "o", ms=3.6, mfc="none", color="k",
                 label="0b8df42, corrected upstream")
    ax1.set_xlabel("$t$")
    ax1.set_ylabel("$a(t)$")
    ax1.set_xlim(0.0, t_new.max())
    cm_x(ax1)
    ax1.legend(loc="lower right", fontsize=8)

    # Relative difference at the times sampled by the newer calculation.
    a_interp = np.interp(t_new, t_old, a_old)
    rel = 100.0 * np.abs(a_new - a_interp) / a_interp
    finite = rel[1:]
    positive = finite[finite > 0.0]

    # Drawn as a distribution rather than as a time series. Plotted against time and joined by a line, a quantity that is pure round-off appears to oscillate, and the appearance is an artefact of two things: the excursions are differences between adjacent representable doubles, not a physical signal, and the samples that happen to agree exactly fall to zero, which a logarithmic axis can only render by clipping them at its floor. A reader is then invited to interpret a sawtooth that means nothing. What the panel is actually evidence for is that every difference lies far below the measurement, and a histogram says that without implying a time dependence.
    if positive.size:
        edges = np.logspace(np.log10(positive.min()),
                            np.log10(positive.max()), 18)
        ax2.hist(positive, bins=edges, color="0.6", edgecolor="k")
        ax2.set_xscale("log")
        ax2.axvline(float(positive.max()), ls="--", lw=GUIDE_LW, color="k")
        ax2.annotate(f"largest {positive.max():.2g} per cent",
                     xy=(float(positive.max()), 0.0),
                     xytext=(0.97, 0.90), textcoords="axes fraction",
                     ha="right", fontsize=8,
                     arrowprops=dict(arrowstyle="->", lw=0.8))
    n_exact = int(finite.size - positive.size)
    ax2.set_xlabel("relative difference (per cent)")
    ax2.set_ylabel("number of snapshots")
    ax2.set_title(f"{finite.size} snapshots, of which {n_exact} agree exactly",
                  fontsize=9)
    # The abscissa is logarithmic and takes no tick formatter; the ordinate is a linear count and takes one.
    cm_y(ax2)

    fig.tight_layout()
    savefig_all(fig, "fig2_replication", outdir=outdir)
    plt.close(fig)
    return float(rel[1:].max())


def figure_proportionality(records, outdir):
    """Collapse of the scaled amplitude and dependence of sigma on epsilon."""
    eps_values = sorted(records)
    t_common = min(records[e][0].max() for e in eps_values)
    grid = np.linspace(0.3, t_common, 60)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.0))

    for e in eps_values:
        t, a = records[e]
        ax1.plot(grid, np.interp(grid, t, a) / e,
                 label=rf"$\varepsilon={e:g}$")
    ax1.set_xlabel("$t$")
    ax1.set_ylabel(r"$a(t)/\varepsilon$")
    ax1.set_xlim(grid.min(), grid.max())
    cm_both(ax1)
    ax1.legend(loc="upper left", fontsize=8)

    scaled = np.vstack([np.interp(grid, records[e][0], records[e][1]) / e
                        for e in eps_values])
    reference = scaled.mean(axis=0)
    for e, row in zip(eps_values, scaled):
        ax2.plot(grid, 100.0 * (row - reference) / reference,
                 label=rf"$\varepsilon={e:g}$")
    ax2.axhline(0.0, ls=":", lw=GUIDE_LW, color="0.5")
    ax2.set_xlabel("$t$")
    ax2.set_ylabel("departure from the mean (per cent)")
    ax2.set_xlim(grid.min(), grid.max())
    cm_both(ax2)
    ax2.legend(loc="upper left", fontsize=8)

    fig.tight_layout()
    savefig_all(fig, "fig3_proportionality", outdir=outdir)
    plt.close(fig)

    spread = np.abs(scaled - reference) / reference
    return float(spread.max()), float(spread.mean())


def figure_dispersion(sigma_measured, outdir):
    """Measured point against the analytic dispersion relation."""
    x = np.linspace(1e-3, 1.0, 400)
    inviscid = np.array([float(dispersion.growth_rate(k, 1.0, 1.0, 1.0))
                         for k in x])
    viscous = np.array([float(dispersion.growth_rate_viscous(
        k, 1.0, 1.0, 1.0, 0.02)) for k in x])

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.plot(x, inviscid, "-", label="inviscid, Rayleigh (1878)")
    ax.plot(x, viscous, "--", label="viscous estimate, Oh = 0.02")
    ax.plot([KR0], [sigma_measured], "o", ms=7, mfc="none", mew=1.8,
            color="k", label="measured, NS-MFP")
    k_peak, s_peak = dispersion.most_unstable_wavenumber(1.0, 1.0, 1.0)
    ax.axvline(k_peak, ls=":", lw=GUIDE_LW, color="0.5")
    ax.annotate(rf"$kr_0={k_peak:.3f}$", xy=(k_peak, 0.02),
                xytext=(k_peak - 0.30, 0.045), fontsize=9,
                arrowprops=dict(arrowstyle="->", lw=0.8))
    ax.set_xlabel("$k r_0$")
    ax.set_ylabel(r"$\sigma$")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, SIGMA_INVISCID * 1.12)
    cm_both(ax)
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    savefig_all(fig, "fig4_dispersion", outdir=outdir)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datadir", default="runs/records")
    parser.add_argument("--archive",
                        default="validation_data/validation_kr0_0.7.json")
    parser.add_argument("--outdir", default="figures")
    args = parser.parse_args(argv)

    apply_style()
    os.makedirs(args.outdir, exist_ok=True)

    archive = json.load(open(args.archive))
    t_old = np.asarray(archive["times"])
    a_old = np.asarray(archive["amplitude"])

    slope, r2 = figure_validation(t_old, a_old, args.outdir)
    print(f"fig1  sigma = {slope:.4f}, r^2 = {r2:.6f}")

    new = np.loadtxt(os.path.join(args.datadir, "amp_eps0.02.dat"))
    worst = figure_replication(new[:, 0], new[:, 1], t_old, a_old,
                               args.outdir)
    print(f"fig2  worst relative difference = {worst:.4f} per cent")

    records = {}
    for eps, name in [(0.01, "amp_eps0.01.dat"),
                      (0.02, "amp_eps0.02.dat"),
                      (0.04, "amp_eps0.04.dat")]:
        d = np.loadtxt(os.path.join(args.datadir, name))
        records[eps] = (d[:, 0], d[:, 1])
    smax, smean = figure_proportionality(records, args.outdir)
    print(f"fig3  scaled spread: max {100*smax:.3f}, "
          f"mean {100*smean:.3f} per cent")

    figure_dispersion(slope, args.outdir)
    print("fig4  dispersion relation with the measured point")
    print(f"figures written to {args.outdir}/")


if __name__ == "__main__":
    main()
