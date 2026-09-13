# !/usr/bin/env python3
"""Publication-quality figures for the Rayleigh-Plateau NS-MFP analysis.

The script produces the four figures used to report a linear stability measurement:

    growth_kr0_XXX          interface mode amplitude a(t) with the fitted
                            exponential and the two theoretical rates
    local_sigma_kr0_XXX     local growth rate d(ln a)/dt against time
    interface_kr0_XXX       interface radius r(z) at several instants
    field_<name>_kr0_XXX    a plotfile variable in the (r,z) plane with the
                            reconstructed interface overlaid

Every figure is written as PNG, PDF and SVG through `savefig_all`, and all typographic settings, including line width, come from `pub_style`; nothing is restated here.

Three input paths are supported, and they differ only in where the numbers come from.

`--rundir` reads the plotfiles directly. This requires yt and requires the plotfiles, so it is used on the machine that ran the calculation, and it is the only path that draws all four figures without a prior export step.

`--bundle` reads a directory written by `dump_amplitude` and `dump_profiles`, holding `amplitude.dat`, `profiles.dat` and any `field_<name>.dat`. This is the path intended for routine use. It draws the same four figures from a few hundred kilobytes of plain text, needs neither yt nor the plotfiles, and therefore works on a workstation while the plotfiles remain on the cluster that produced them.

`--json` reads a single amplitude record and draws the two amplitude figures only, since a record carries no interface profile. It is retained because the archived validation data is in that form.

Run from `amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Fortran`:

    python3 plotting/plot_results.py --bundle runs/bundles/k0.7 --kr0 0.7
    python3 plotting/plot_results.py --json validation_data/validation_kr0_0.7.json
    python3 plotting/plot_results.py --rundir runs/k0.7_eps0.02 \\
        --wavelength 8.976 --kr0 0.7 --field y_velocity
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import bundle as bundle_io
import dispersion
import transient
from pub_style import (apply_style, savefig_all, cm_both, cm_x,
                       CMTickFormatter, PALETTE)

#: Width of a reference or guide line, matched to the axis frame so that a guide never competes with a data curve. Data curves take lines.linewidth from pub_style, which is 1.0, and are therefore drawn without an explicit width anywhere in this file.
GUIDE_LW = 0.8


def _load_json(path):
    """Read an amplitude record written by a previous analysis run."""
    with open(path) as handle:
        record = json.load(handle)
    times = np.asarray(record["times"], dtype=float)
    amps = np.asarray(record["amplitude"], dtype=float)
    means = np.asarray(record.get("mean_radius", []), dtype=float)
    return times, amps, means


def _load_rundir(rundir, wavelength, ls_field, prefix):
    """Read the plotfiles directly.

    The import is deferred so that the bundle and JSON paths remain usable on a machine without yt installed.
    """
    import interface_mode as imode

    return imode.amplitude_series(
        rundir, wavelength, level_set_field=ls_field, prefix=prefix
    )


def _fit(times, amps, t_min):
    """Least-squares slope of log a(t) over t >= t_min.

    Duplicated here rather than imported from `interface_mode` so that the bundle and JSON paths carry no dependency on yt.
    """
    mask = np.isfinite(amps) & (amps > 0) & (times >= t_min)
    if mask.sum() < 3:
        raise ValueError("fewer than three usable samples in the fit window")
    t, y = times[mask], np.log(amps[mask])
    slope, intercept = np.polyfit(t, y, 1)
    residual = y - (slope * t + intercept)
    ss_res = float(np.sum(residual**2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(slope), float(intercept), r2


def _local_sigma(times, amps, width):
    """Growth rate estimated over sliding windows of the given width.

    The curve, rather than any single fitted number, indicates whether the asymptotic regime has been reached.
    """
    lo = times.min() + 0.5 * width
    hi = times.max() - 0.5 * width
    if hi <= lo:
        return np.array([]), np.array([])
    centres, sigmas = [], []
    for centre in np.linspace(lo, hi, 8):
        mask = (times >= centre - 0.5 * width) & (times <= centre + 0.5 * width)
        mask &= np.isfinite(amps) & (amps > 0)
        if mask.sum() >= 4:
            centres.append(centre)
            sigmas.append(np.polyfit(times[mask], np.log(amps[mask]), 1)[0])
    return np.asarray(centres), np.asarray(sigmas)


def figure_growth(times, amps, t_fit, sigma_inv, sigma_vis, kr0, tag, outdir):
    """Amplitude against time on a logarithmic ordinate.

    The transient is shaded and labelled rather than discarded, so that the extent of the excluded record remains visible to a reader who would otherwise wonder why the fit begins where it does. Both theoretical slopes are anchored to the fit at its left edge, which makes the comparison one of gradient rather than of offset.
    """
    slope, intercept, r2 = _fit(times, amps, t_fit)

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.semilogy(times, amps, "o", ms=3.2, mfc="none", color="k",
                label="measured")

    mask = times >= t_fit
    t_win = times[mask]
    ax.semilogy(t_win, np.exp(slope * t_win + intercept), "-",
                label=rf"fit, $\sigma={slope:.4f}$ ($r^2={r2:.5f}$)")

    anchor = float(np.exp(slope * t_win[0] + intercept))
    ax.semilogy(t_win, anchor * np.exp(sigma_vis * (t_win - t_win[0])),
                "--", label=rf"viscous exact, $\sigma={sigma_vis:.4f}$")
    ax.semilogy(t_win, anchor * np.exp(sigma_inv * (t_win - t_win[0])),
                ":", label=rf"inviscid, $\sigma={sigma_inv:.4f}$")

    ax.axvspan(times.min(), t_fit, color="0.88", zorder=0, lw=0)
    finite = amps[np.isfinite(amps) & (amps > 0)]
    if finite.size:
        ax.text(0.5 * (times.min() + t_fit),
                float(np.sqrt(finite.min() * finite.max())),
                "transient\n(excluded)", ha="center", va="center",
                color="0.35")
    ax.set_xlabel("$t$")
    ax.set_ylabel("interface mode amplitude $a(t)$")
    ax.set_title(rf"Rayleigh-Plateau growth, $k r_0 = {kr0:g}$")
    ax.set_xlim(times.min(), times.max())
    # The ordinate is logarithmic and therefore takes no tick formatter; the scalar formatter would replace the powers of ten with plain numbers.
    cm_x(ax)
    ax.legend(loc="lower right")
    fig.tight_layout()
    savefig_all(fig, f"growth_{tag}", outdir=outdir)
    plt.close(fig)
    return slope, r2


def figure_local_sigma(times, amps, width, sigma_inv, sigma_vis, kr0, tag,
                       outdir, settled=True):
    """Local growth rate against time.

    This figure justifies the choice of fit window: the estimate rises out of the startup transient and flattens once modal growth is established. A curve still rising at the right-hand edge indicates that the calculation was stopped before the asymptotic regime, and no rate from it should be quoted.
    """
    centres, sigmas = _local_sigma(times, amps, width)
    if centres.size == 0:
        print("record too short for a local growth rate; figure skipped")
        return centres, sigmas

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.plot(centres, sigmas, "o-", ms=4.0, label=r"local $\sigma(t)$")
    ax.axhline(sigma_vis, ls="--", color="0.35", label="viscous (exact)")
    ax.axhline(sigma_inv, ls=":", color="0.55", label="inviscid")
    ax.set_xlabel("$t$")
    ax.set_ylabel(r"$\mathrm{d}(\ln a)/\mathrm{d}t$")
    verdict = ("plateau reached" if settled else
               "STILL RISING: the run is too short")
    ax.set_title(rf"Local growth rate, $k r_0 = {kr0:g}$"
                 f"\n({verdict})")
    ax.set_xlim(float(centres.min()), float(centres.max()))
    ax.set_ylim(0.0, max(sigma_inv, float(sigmas.max())) * 1.12)
    cm_both(ax)
    ax.legend(loc="lower right")
    fig.tight_layout()
    savefig_all(fig, f"local_sigma_{tag}", outdir=outdir)
    plt.close(fig)
    return centres, sigmas


def figure_interface(z, times, radii, instants, wavelength, tag, outdir,
                     r_limit=None, r_window=(0.825, 1.175)):
    """Interface radius against axial coordinate at selected instants.

    The ordinate is fixed rather than scaled to the data, for two reasons. A figure whose axis moves between cases cannot be compared across them; and auto-scaling to a record that includes break-up compresses the linear growth this figure exists to show into a band a few per cent tall. The default window spans the amplitude range the linear description covers, so the growth of the mode fills the frame.

    Each instant takes its own colour from the house palette, assigned explicitly by position rather than left to the property cycle, since a cycle advanced by earlier artists on the same axes does not begin at the first colour.
    """
    times = np.asarray(times, dtype=float)
    radii = np.asarray(radii, dtype=float)

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for order, target in enumerate(instants):
        index = int(np.argmin(np.abs(times - target)))
        ax.plot(z, radii[index], color=PALETTE[order % len(PALETTE)],
                label=rf"$t={times[index]:.2f}$")
    ax.axhline(1.0, ls=":", lw=GUIDE_LW, color="0.45")

    ax.set_xlabel("$z$")
    ax.set_ylabel(r"$r_{\mathrm{interface}}$")
    ax.set_title("Interface shape: the capillary mode growing")
    ax.set_xlim(0.0, float(wavelength))
    if r_limit is not None:
        ax.set_ylim(1.0 - r_limit, 1.0 + r_limit)
    else:
        ax.set_ylim(*r_window)
    cm_both(ax)
    # The curves cross at the nodes and reach their extremes at the antinodes, so at six instants no region inside the axes is free. The legend therefore sits outside the frame, which also keeps the zero margins of the house style intact.
    ncol = len(instants) if len(instants) <= 6 else 4
    ax.legend(ncol=ncol, loc="upper center", bbox_to_anchor=(0.5, -0.16),
              frameon=False, columnspacing=1.6, handlelength=1.8)
    fig.tight_layout()
    savefig_all(fig, f"interface_{tag}", outdir=outdir,
                bbox_extra_artists=(ax.get_legend(),))
    plt.close(fig)


def figure_field(field, r, z, times, planes, radii_z, radii, tag, outdir,
                 amplitudes=None):
    """Filled contours of a plotfile variable with the interface overlaid.

    The panels share one colour scale so that the reader compares magnitudes rather than the widths of individual colour bands, and the scale is symmetric about zero so that the sign of the variable is read from the hue.

    Parameters
    ----------
    planes : array of shape (n_panel, nr, nz)
        The variable on the (r,z) plane at each instant.
    radii_z, radii : arrays
        Axial stations, and the interface radius at each panel's instant, drawn as the black contour.
    """
    planes = np.asarray(planes, dtype=float)
    n_panel = planes.shape[0]
    vmax = float(np.abs(planes).max())
    if vmax <= 0.0:
        vmax = 1.0

    # Panel width shrinks a little as panels are added, so that four still fit across a page at the margins used for the report rather than being scaled down as a block.
    per_panel = 2.7 if n_panel <= 3 else 2.3
    fig, axes = plt.subplots(1, n_panel, figsize=(per_panel * n_panel, 4.6),
                             sharey=True)
    axes = np.atleast_1d(axes)

    levels = np.linspace(-vmax, vmax, 21)
    r_limit = min(2.0, float(r.max()))
    for index, ax in enumerate(axes):
        contour = ax.contourf(r, z, planes[index].T, levels=levels,
                              cmap="RdBu_r", extend="both")
        ax.plot(radii[index], radii_z, "-", color="k")
        ax.set_xlabel("$r$")
        # The mode amplitude is stated on each panel. In the linear regime the interface displacement is a few per cent of the radius, which on a panel spanning two radii is a line width or two, so a reader who is not told the amplitude sees an apparently straight interface and reasonably wonders whether anything is happening. The number resolves that: the deformation is small because the measurement requires it to be small.
        if amplitudes is not None:
            ax.set_title(rf"$t={times[index]:.2f}$"
                         "\n" rf"$a/r_0={amplitudes[index]:.4f}$",
                         fontsize=9)
        else:
            ax.set_title(rf"$t={times[index]:.2f}$")
        ax.set_xlim(0.0, r_limit)
        ax.set_ylim(float(z.min()), float(z.max()))
        # The global style draws a dashed grid, which would otherwise be overlaid on the filled contours.
        ax.grid(False)
        cm_both(ax)
    axes[0].set_ylabel("$z$")
    # Computer Modern carries no underscore glyph, so a raw plotfile variable name such as y_velocity renders with a dot in place of the underscore. The label therefore uses a space while the file name keeps the variable name intact.
    label = field.replace("_", " ")
    fig.suptitle(f"{label} (colour) with the interface (black)")
    bar = fig.colorbar(contour, ax=axes, shrink=0.75, pad=0.03, label=label)
    bar.locator = MaxNLocator(nbins=7, symmetric=True)
    bar.formatter = CMTickFormatter()
    bar.update_ticks()
    savefig_all(fig, f"field_{field}_{tag}", outdir=outdir)
    plt.close(fig)


def _profiles_from_rundir(rundir, prefix, ls_field):
    """Interface profiles read from the plotfiles, in the bundle's shape."""
    import interface_mode as imode
    import nsmfp

    entries = []
    for path in nsmfp.find_plotfiles(rundir, prefix=prefix):
        t, z, radius = imode.interface_radius(path, level_set_field=ls_field)
        entries.append((t, z, radius))
    entries.sort(key=lambda item: item[0])
    times = np.array([item[0] for item in entries])
    radii = np.vstack([item[2] for item in entries])
    return entries[0][1], times, radii


def _field_from_rundir(rundir, field, instants, prefix):
    """Field slices read from the plotfiles, in the bundle's shape."""
    import nsmfp

    stamps = []
    for path in nsmfp.find_plotfiles(rundir, prefix=prefix):
        pf = nsmfp.read_plotfile(path, fields=(field,))
        stamps.append((pf.time, path))
    stamps.sort()
    all_times = np.array([item[0] for item in stamps])

    planes, times, r, z = [], [], None, None
    for target in instants:
        index = int(np.argmin(np.abs(all_times - target)))
        t, path = stamps[index]
        pf = nsmfp.read_plotfile(path, fields=(field,))
        if r is None:
            n_r, n_z = pf.shape[0], pf.shape[1]
            r = (np.arange(n_r) + 0.5) * float(pf.domain_right[0]) / n_r
            z = (np.arange(n_z) + 0.5) * float(pf.domain_right[1]) / n_z
        planes.append(pf.data.reshape(pf.shape)[:, :, 0])
        times.append(t)
    return r, z, np.array(times), np.stack(planes)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--rundir", help="directory holding the plotfiles")
    source.add_argument("--bundle",
                        help="directory of exported figure data, holding "
                             "amplitude.dat, profiles.dat and field_*.dat")
    source.add_argument("--json", help="amplitude record from a previous run")

    parser.add_argument("--wavelength", type=float, default=None)
    parser.add_argument("--kr0", type=float, default=0.7)
    parser.add_argument("--r0", type=float, default=1.0)
    parser.add_argument("--sigma", type=float, default=1.0)
    parser.add_argument("--rho", type=float, default=1.0)
    parser.add_argument("--mu", type=float, default=0.02)
    parser.add_argument("--r-min", type=float, default=0.825, dest="r_min",
                        help="lower limit of the interface figure")
    parser.add_argument("--r-max", type=float, default=1.175, dest="r_max",
                        help="upper limit of the interface figure; fixed "
                             "rather than scaled so that cases compare")
    parser.add_argument("--linear-bound", type=float, default=0.1,
                        dest="linear_bound",
                        help="amplitude, as a fraction of r0, beyond which "
                             "the record is treated as nonlinear")
    parser.add_argument("--rho-outer", type=float, default=1.225e-3,
                        dest="rho_outer",
                        help="ambient density, for the two-fluid reference")
    parser.add_argument("--t-fit", type=float, dest="t_fit", default=None,
                        help="start of the fit window; default 0.65 t_max")
    parser.add_argument("--width", type=float, default=0.8,
                        help="sliding window width for the local rate")
    parser.add_argument("--prefix", default="nddataPLT")
    parser.add_argument("--ls-field", dest="ls_field", default="L0101")
    parser.add_argument("--field", default=None,
                        help="also contour this plotfile variable; with "
                             "--bundle, defaults to every field exported")
    parser.add_argument("--field-times", type=float, nargs="+",
                        dest="field_times", default=None)
    parser.add_argument("--n-field", type=int, default=3, dest="n_field",
                        help="number of panels in the field figure; the "
                             "instants are spaced by equal amplitude")
    parser.add_argument("--interface-times", type=float, nargs="+",
                        dest="interface_times", default=None,
                        help="instants for the interface figure; "
                             "default six spread across the record")
    parser.add_argument("--n-interface", type=int, default=6,
                        dest="n_interface",
                        help="number of interface curves when the instants "
                             "are not named")
    parser.add_argument("--outdir", default="figures")
    args = parser.parse_args(argv)

    apply_style()
    os.makedirs(args.outdir, exist_ok=True)

    profiles = None
    fields = {}
    wavelength = args.wavelength

    if args.json:
        times, amps, means = _load_json(args.json)
    elif args.bundle:
        record = os.path.join(args.bundle, "amplitude.dat")
        if not os.path.isfile(record):
            parser.error(f"no amplitude.dat in {args.bundle}; write one with "
                         "python/dump_amplitude.py")
        times, amps, means = bundle_io.read_amplitude(record)
        profile_path = os.path.join(args.bundle, "profiles.dat")
        if os.path.isfile(profile_path):
            profiles = bundle_io.read_profiles(profile_path)
            if wavelength is None:
                wavelength = profiles["wavelength"]
        wanted = [args.field] if args.field else None
        for name, path in bundle_io.find_field_bundles(args.bundle).items():
            if wanted is None or name in wanted:
                fields[name] = bundle_io.read_field(path)
    else:
        if args.wavelength is None:
            parser.error("--wavelength is required together with --rundir")
        times, amps, means = _load_rundir(args.rundir, args.wavelength,
                                          args.ls_field, args.prefix)

    # Bounded on both sides. The upper limit is where the amplitude leaves the linear range; a window set as a fixed fraction of the record, as this once was, lands in the nonlinear wreckage whenever a run is carried past linearity, and returns a rate belonging to neither regime.
    t_lo, t_hi, nonlinear, wide = transient.select_window(
        times, amps, r0=args.r0, bound=args.linear_bound)
    t_fit = args.t_fit if args.t_fit is not None else t_lo
    if args.t_fit is not None:
        t_hi = times.max()
    wavenumber = args.kr0 / args.r0
    sigma_inv = float(dispersion.growth_rate(wavenumber, r0=args.r0,
                                             sigma=args.sigma, rho=args.rho))
    sigma_vis = float(dispersion.growth_rate_viscous_exact(
        wavenumber, r0=args.r0, sigma=args.sigma, rho=args.rho, mu=args.mu))
    tag = f"kr0_{args.kr0:g}".replace(".", "p")

    # The inviscid rate bounds any admissible estimate, since viscosity can only retard a capillary instability. Passing it lets the diagnostic reject an extrapolation that violates it instead of reporting an impossible number.
    diagnosis = transient.assess(times, amps, width=args.width,
                                 sigma_reference=sigma_vis, r0=args.r0,
                                 bound=args.linear_bound,
                                 sigma_max=sigma_inv)

    # Every figure and every fitted number uses the linear part only.
    keep = times <= t_hi
    times, amps = times[keep], amps[keep]
    if means.size:
        means = means[keep]

    slope, r2 = figure_growth(times, amps, t_fit, sigma_inv, sigma_vis,
                              args.kr0, tag, args.outdir)
    figure_local_sigma(times, amps, args.width, sigma_inv, sigma_vis,
                       args.kr0, tag, args.outdir,
                       settled=diagnosis["settled"])

    # The interface figure carries more curves than the field figure. Six instants show the mode emerging from the transient, growing through the linear regime and steepening as it approaches the nonlinear bound, which is the sequence the figure exists to convey; three would show the endpoints and lose the shape of the approach. The field figure stays at three because each instant is a full panel and a row of six is unreadable at column width. Both are spaced from the first sample rather than from zero, so that no curve is drawn from a snapshot that does not exist.
    interface_times = args.interface_times or list(
        np.linspace(times.min(), times.max(), args.n_interface))
    # Spaced by equal steps in amplitude rather than in time. The mode grows exponentially, so equally spaced instants give a geometric progression of amplitudes: the early panels differ from the initial condition by a fraction of a line width and look identical, while the whole change is crowded into the last one. Stepping through equal amplitudes instead makes every panel visibly different from its neighbour. The first panel is kept at the initial instant as a reference against which the rest are read.
    instants = args.field_times
    if instants is None:
        n_panel = max(int(args.n_field), 2)
        finite = np.isfinite(amps) & (amps > 0)
        a_end = float(amps[finite].max()) if finite.any() else 1.0
        instants = [float(times.min())]
        for step in range(1, n_panel):
            target = a_end * step / (n_panel - 1)
            index = int(np.argmin(np.abs(amps - target)))
            instants.append(float(times[index]))

    if args.rundir:
        z, p_times, radii = _profiles_from_rundir(args.rundir, args.prefix,
                                                  args.ls_field)
        keep_p = p_times <= t_hi
        figure_interface(z, p_times[keep_p], radii[keep_p], interface_times,
                         wavelength, tag, args.outdir,
                         r_window=(args.r_min, args.r_max))
        if args.field:
            r, fz, f_times, planes = _field_from_rundir(
                args.rundir, args.field, instants, args.prefix)
            picked = [int(np.argmin(np.abs(p_times - t))) for t in f_times]
            panel_amps = np.interp(f_times, times, amps)
            figure_field(args.field, r, fz, f_times, planes, z,
                         radii[picked], tag, args.outdir,
                         amplitudes=panel_amps)
    elif profiles is not None:
        # Truncated with everything else: the profile bundle spans the whole calculation, and a figure of the linear mode must not include the break-up that follows it.
        keep_p = profiles["times"] <= t_hi
        figure_interface(profiles["z"], profiles["times"][keep_p],
                         profiles["radius"][keep_p],
                         interface_times, wavelength, tag, args.outdir,
                         r_window=(args.r_min, args.r_max))
        for name, data in fields.items():
            # The bundle holds whatever slices were exported, and the export spaces them evenly in time. Drawing all of them reproduces the crowding this figure was changed to avoid, so the panels are chosen here from what is available, by the same equal-amplitude rule used for a run directory. Exporting more slices than are drawn therefore costs a little disk and buys a better choice.
            slice_amps = np.interp(data["times"], times, amps)
            chosen = [0]
            n_panel = max(int(args.n_field), 2)
            if args.field_times is not None:
                chosen = [int(np.argmin(np.abs(data["times"] - t)))
                          for t in args.field_times]
            elif data["times"].size > n_panel:
                a_end = float(slice_amps.max())
                for step in range(1, n_panel):
                    target = a_end * step / (n_panel - 1)
                    index = int(np.argmin(np.abs(slice_amps - target)))
                    if index not in chosen:
                        chosen.append(index)
            else:
                chosen = list(range(data["times"].size))
            chosen = sorted(set(chosen))
            if len(chosen) < n_panel and data["times"].size >= n_panel:
                print(f"note: {name} bundle holds {data['times'].size} "
                      f"slices; {len(chosen)} distinct panels available. "
                      f"Re-export with more slices for {n_panel}.")

            panel_times = data["times"][chosen]
            picked = [int(np.argmin(np.abs(profiles["times"] - t)))
                      for t in panel_times]
            figure_field(name, data["r"], data["z"], panel_times,
                         data["data"][chosen], profiles["z"],
                         profiles["radius"][picked], tag, args.outdir,
                         amplitudes=slice_amps[chosen])

    if nonlinear:
        print(f"linear range ends at t = {t_hi:.3f}; later snapshots are "
              f"excluded from every fit and figure")
    print(f"snapshots        {len(times)} (t up to {times.max():.3f})")
    print(f"amplitude        {amps[0]:.5f} -> {amps[-1]:.5f}"
          f"  (a/r0 = {amps[-1] / args.r0:.3f})")
    if means.size:
        drift = 100.0 * (means[-1] - means[0]) / means[0]
        print(f"mean radius      {drift:+.3f} per cent")
    print(f"sigma (t > {t_fit:.2f})  {slope:.4f}   r^2 = {r2:.6f}")
    print(f"viscous estimate {sigma_vis:.4f}"
          f"   ({100 * (slope - sigma_vis) / sigma_vis:+.2f} per cent)")
    print(f"inviscid         {sigma_inv:.4f}"
          f"   ({100 * (slope - sigma_inv) / sigma_inv:+.2f} per cent)")
    corrected = diagnosis.get("sigma_inf", float("nan"))
    if np.isfinite(corrected):
        print(f"transient-corrected {corrected:.4f}"
              f"   ({100 * (corrected - sigma_vis) / sigma_vis:+.2f} per cent "
              f"vs exact viscous)")
    print()
    print(transient.format_report(diagnosis, sigma_reference=sigma_vis))
    print(f"figures written to {args.outdir}/")


if __name__ == "__main__":
    main()
