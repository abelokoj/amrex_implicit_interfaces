# !/usr/bin/env python3
"""Break-up of the column, and the linear rate extrapolated to it.

The linear analysis returns a growth rate. What a reader wants to know is whether that rate predicts anything, and the most direct test available is the time at which the column breaks. Extrapolating exponential growth from the linear regime to the moment the minimum radius reaches zero gives a predicted break-up time; running the calculation on into the nonlinear regime gives an observed one. The comparison closes the loop between a linear eigenvalue and a physical outcome, and it is the most legible result in the study for an audience that does not work in stability theory.

The prediction is expected to be an underestimate of the observed time rather than an exact match, and the reason is physical rather than numerical. Exponential growth is the linear solution; as the neck thins, the local curvature and hence the driving capillary pressure change, and the approach to pinch-off follows a self-similar law rather than an exponential. Reporting the discrepancy and naming its cause is a stronger position than reporting agreement, and a figure that concealed the departure by fitting over the nonlinear range would be worse than useless.

Three quantities are drawn. The minimum and maximum interface radius against time show the neck thinning while the bulge grows, and their separation is the visible signature of nonlinearity. The mode amplitude against time carries the linear fit extrapolated to the axis, with the predicted and observed break-up times marked. And a sequence of interface profiles through break-up shows the shape the column takes, which no scalar conveys.

Everything is read from a profile bundle, so this script needs neither yt nor the plotfiles. The run itself is a different matter, and the reference configuration will not do. Two settings must change.

The calculation must be carried past break-up. Since `make_inputs.py` measures its stop time in e-folding times of the viscous rate, and the reference amplitude of 2e-2 must grow by a factor of fifty before the neck closes, break-up falls near t = 14 at the measured rate. Allowing for the nonlinear slowdown, `--n-periods 6` gives a stop time near 18 and a margin of roughly a quarter, which is the recommended default.

The mesh must resolve the neck. At 8 cells per radius a cell is an eighth of the radius, so a neck of any interesting thinness falls below one cell and the reconstruction cannot represent it; a break-up time measured there is a statement about the mesh. Thirty-two cells per radius is the recommended default, at which two cell widths are 6 per cent of the radius, and the case should be repeated at 16 so that the sensitivity of the break-up time to the mesh is reported rather than assumed. That sensitivity is the first thing a reader should ask about, since break-up is a singularity and no mesh resolves it to the end.

Run from `amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Python`:

    python3 python/plot_pinchoff.py --bundle runs/bundles/k0.7_long --kr0 0.7 --outdir figures
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import bundle as bundle_io
from pub_style import apply_style, savefig_all, cm_both, cm_x

GUIDE_LW = 0.8


def extrema(times, radius):
    """Minimum and maximum interface radius at each instant."""
    r_min = np.nanmin(radius, axis=1)
    r_max = np.nanmax(radius, axis=1)
    return r_min, r_max


def mode_amplitude(z, radius, wavelength):
    """Amplitude of the first axial harmonic at each instant."""
    k = 2.0 * np.pi / wavelength
    kernel = np.exp(-1j * k * z)
    amps = np.empty(radius.shape[0])
    for i, row in enumerate(radius):
        good = np.isfinite(row)
        if good.sum() < 4:
            amps[i] = np.nan
            continue
        centred = row[good] - row[good].mean()
        amps[i] = 2.0 * np.abs((centred * kernel[good]).mean())
    return amps


def linear_fit(times, amps, t_lo, t_hi):
    """Slope and intercept of log a over the stated window."""
    mask = (np.isfinite(amps) & (amps > 0) & (times >= t_lo)
            & (times <= t_hi))
    if mask.sum() < 3:
        raise ValueError("fewer than three usable samples in the fit window")
    slope, intercept = np.polyfit(times[mask], np.log(amps[mask]), 1)
    return float(slope), float(intercept)


def predicted_breakup(slope, intercept, r0):
    """Time at which extrapolated linear growth would close the neck.

    The neck closes when the mode amplitude reaches the unperturbed radius, since the minimum radius of r0 + a cos(k z) is r0 - a. That criterion is exact within the linear description and is the natural extrapolation; it is not a claim that the column actually breaks that way.
    """
    return (np.log(r0) - intercept) / slope


def observed_breakup(times, r_min, threshold):
    """First instant at which the minimum radius falls below a threshold.

    A threshold rather than zero, because the interface is reconstructed on a mesh and the last representable neck is of order one cell. The threshold is reported with the time so that the two are read together.
    """
    below = np.where(r_min <= threshold)[0]
    if below.size == 0:
        return None
    index = int(below[0])
    if index == 0:
        return float(times[0])
    # Linear interpolation between the bracketing samples, so that the answer is not quantised onto the snapshot interval.
    t0, t1 = times[index - 1], times[index]
    r_a, r_b = r_min[index - 1], r_min[index]
    if r_a == r_b:
        return float(t1)
    return float(t0 + (r_a - threshold) * (t1 - t0) / (r_a - r_b))


def figure_extrema(times, r_min, r_max, r0, t_obs, threshold, outdir):
    """Minimum and maximum interface radius against time."""
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(times, r_min, "-", label=r"minimum radius (neck)")
    ax.plot(times, r_max, "-", label=r"maximum radius (bulge)")
    ax.axhline(r0, ls=":", lw=GUIDE_LW, color="0.45",
               label="unperturbed radius")
    if t_obs is not None:
        ax.axvline(t_obs, ls="--", lw=GUIDE_LW, color="0.35",
                   label=rf"break-up, $t={t_obs:.3f}$")
    ax.set_xlabel("$t$")
    ax.set_ylabel(r"$r_{\mathrm{interface}}$")
    ax.set_title("Neck thinning and bulge growth")
    ax.set_xlim(float(times.min()), float(times.max()))
    ax.set_ylim(0.0, float(np.nanmax(r_max)) * 1.05)
    cm_both(ax)
    ax.legend(loc="upper left")
    fig.tight_layout()
    savefig_all(fig, "pinchoff_extrema", outdir=outdir)
    plt.close(fig)


def figure_extrapolation(times, amps, slope, intercept, t_lo, t_hi,
                         t_pred, t_obs, r0, outdir):
    """Mode amplitude with the linear fit extrapolated to break-up."""
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.semilogy(times, amps, "o", ms=3.2, mfc="none", color="k",
                label="measured")

    span = np.linspace(t_lo, max(t_pred, times.max()), 200)
    ax.semilogy(span, np.exp(slope * span + intercept), "-",
                label=rf"linear fit, $\sigma={slope:.4f}$")
    ax.axhline(r0, ls=":", lw=GUIDE_LW, color="0.45",
               label=r"$a = r_0$, neck closed")
    ax.axvspan(t_lo, t_hi, color="0.88", zorder=0, lw=0)

    ax.axvline(t_pred, ls="--", lw=GUIDE_LW, color="0.35",
               label=rf"predicted, $t={t_pred:.3f}$")
    if t_obs is not None:
        ax.axvline(t_obs, ls="-.", lw=GUIDE_LW, color="0.2",
                   label=rf"observed, $t={t_obs:.3f}$")

    ax.set_xlabel("$t$")
    ax.set_ylabel("interface mode amplitude $a(t)$")
    ax.set_title("Linear growth extrapolated to break-up")
    ax.set_xlim(float(times.min()), max(t_pred, float(times.max())) * 1.02)
    cm_x(ax)
    ax.legend(loc="lower right")
    fig.tight_layout()
    savefig_all(fig, "pinchoff_extrapolation", outdir=outdir)
    plt.close(fig)


def figure_sequence(z, times, radius, instants, wavelength, outdir):
    """Interface profiles through break-up."""
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    shown = []
    for target in instants:
        index = int(np.argmin(np.abs(times - target)))
        ax.plot(z, radius[index], label=rf"$t={times[index]:.2f}$")
        shown.append(radius[index])
    ax.axhline(0.0, lw=GUIDE_LW, color="0.45")
    ax.set_xlabel("$z$")
    ax.set_ylabel(r"$r_{\mathrm{interface}}$")
    ax.set_title("Interface through break-up")
    ax.set_xlim(0.0, float(wavelength))
    # Scaled to the profiles actually drawn, not to the whole record. The bulge continues to grow after break-up, so the record maximum would leave a third of the panel empty and shrink the very curves the figure exists to show.
    ax.set_ylim(0.0, float(np.nanmax(np.vstack(shown))) * 1.06)
    cm_both(ax)
    # Placed outside the frame for the same reason as the interface figure: near break-up the profiles span the whole ordinate and leave no clear region inside the axes.
    ncol = len(instants) if len(instants) <= 6 else 4
    ax.legend(ncol=ncol, loc="upper center", bbox_to_anchor=(0.5, -0.16),
              frameon=False, columnspacing=1.6, handlelength=1.8)
    fig.tight_layout()
    savefig_all(fig, "pinchoff_sequence", outdir=outdir,
                bbox_extra_artists=(ax.get_legend(),))
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", required=True,
                    help="bundle directory holding profiles.dat")
    ap.add_argument("--kr0", type=float, default=0.7)
    ap.add_argument("--r0", type=float, default=1.0)
    ap.add_argument("--t-lo", type=float, default=None, dest="t_lo",
                    help="start of the linear fit window")
    ap.add_argument("--t-hi", type=float, default=None, dest="t_hi",
                    help="end of the linear fit window; default where the "
                         "amplitude first reaches a tenth of r0")
    ap.add_argument("--threshold", type=float, default=None,
                    help="neck radius counted as break-up; default two cell "
                         "widths, from --cells-per-r0")
    ap.add_argument("--cells-per-r0", type=int, default=32,
                    dest="cells_per_r0",
                    help="radial resolution of the run, used to set the "
                         "break-up threshold at two cell widths")
    ap.add_argument("--n-profiles", type=int, default=6, dest="n_profiles")
    ap.add_argument("--outdir", default="figures")
    args = ap.parse_args(argv)

    apply_style()
    os.makedirs(args.outdir, exist_ok=True)

    profiles = bundle_io.read_profiles(
        os.path.join(args.bundle, "profiles.dat"))
    z = profiles["z"]
    times = profiles["times"]
    radius = profiles["radius"]
    wavelength = profiles["wavelength"]

    r_min, r_max = extrema(times, radius)
    amps = mode_amplitude(z, radius, wavelength)

    # The fit window ends where the amplitude reaches a tenth of the radius, which is the bound beyond which the linear description is not claimed to hold. Fitting past it would raise the apparent rate and flatter the prediction.
    t_hi = args.t_hi
    if t_hi is None:
        beyond = np.where(amps >= 0.1 * args.r0)[0]
        t_hi = float(times[beyond[0]]) if beyond.size else float(times.max())
    t_lo = args.t_lo if args.t_lo is not None else 0.65 * t_hi

    slope, intercept = linear_fit(times, amps, t_lo, t_hi)
    t_pred = predicted_breakup(slope, intercept, args.r0)

    # The break-up threshold is referred to the mesh, not to the radius. A neck thinner than a cell is not represented by the interface reconstruction, so a fixed fraction of r0 asks a question the calculation cannot answer: at 8 cells per radius two cell widths are a quarter of the radius, and no threshold below that is meaningful. Two cell widths is the smallest neck the reconstruction resolves with a cell either side of it.
    cell = args.r0 / max(args.cells_per_r0, 1)
    threshold = (args.threshold if args.threshold is not None
                 else 2.0 * cell)
    t_obs = observed_breakup(times, r_min, threshold)

    figure_extrema(times, r_min, r_max, args.r0, t_obs, threshold, args.outdir)
    figure_extrapolation(times, amps, slope, intercept, t_lo, t_hi,
                         t_pred, t_obs, args.r0, args.outdir)

    # The profile sequence is weighted towards late times, where the shape changes fastest and where the figure earns its place.
    last = t_obs if t_obs is not None else float(times.max())
    instants = list(np.linspace(t_lo, last, args.n_profiles))
    figure_sequence(z, times, radius, instants, wavelength, args.outdir)

    print(f"resolution       {args.cells_per_r0} cells per radius, "
          f"cell width {cell:.4f}")
    print(f"break threshold  {threshold:.4f} (two cell widths)")
    print(f"fit window       {t_lo:.3f} < t < {t_hi:.3f}")
    print(f"linear rate      {slope:.5f}")
    print(f"predicted break  t = {t_pred:.4f}")
    if t_obs is None:
        print(f"observed break   not reached; minimum radius fell only to "
              f"{np.nanmin(r_min):.4f} against a threshold of "
              f"{threshold:.4f}. Advance the calculation further.")
    else:
        print(f"observed break   t = {t_obs:.4f}  "
              f"(neck below {threshold:.4f})")
        print(f"departure        {100 * (t_pred - t_obs) / t_obs:+.2f} "
              f"per cent")
    print(f"figures written to {args.outdir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
