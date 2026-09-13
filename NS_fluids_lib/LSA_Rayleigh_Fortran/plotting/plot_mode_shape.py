# !/usr/bin/env python3
"""The extracted mode shape against the analytic eigenfunction.

Recovering the growth rate shows that the method measures a rate. Recovering the mode shape shows that it measures the mode, which is a stronger claim and a different one: a scalar can agree by coincidence or by compensating errors, whereas a radial profile agreeing point by point across the column cannot. For a study whose contribution is methodological this is the closer of the two arguments, and it is the figure a committee is most likely to ask for when told that a decomposition extracts eigenpairs.

The comparison rests on the classical inviscid solution. For an axisymmetric disturbance of the form exp(sigma t + i k z) on a column of radius r0, the velocity derives from a potential proportional to the modified Bessel function I_0(k r), so that

    u_r(r)  proportional to  I_1(k r)
    u_z(r)  proportional to  I_0(k r)
    p(r)    proportional to  I_0(k r)

with u_z and p in quadrature with u_r along the axis. The radial velocity therefore vanishes on the axis, as it must by symmetry, and both profiles rise monotonically to the interface, where the kinematic condition ties them to the rate of change of the interface radius.

Two caveats are stated rather than hidden. The calculation is viscous at Oh = 0.02, so agreement with an inviscid eigenfunction is expected to be close but not exact, and the departure is largest in a thin layer near the interface where the viscous solution must satisfy a stress condition the inviscid one does not. And the DMD mode is determined only up to a complex scale, so it is normalised here before comparison; normalisation is a choice of amplitude and phase, and it cannot manufacture agreement in shape.

Run from `amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Fortran`:

    python3 plotting/plot_mode_shape.py --perturbed runs/k0.7_eps0.02 --base runs/k0.7_base --outdir figures
    python3 plotting/plot_mode_shape.py --replot --outdir figures
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pub_style import (apply_style, savefig_all, cm_both,
                       save_data, load_data, has_data)

GUIDE_LW = 0.8
DATA = "mode_shape"


def bessel_i(order, x):
    """Modified Bessel function of the first kind.

    Taken from SciPy, which is already required for the decomposition. The dispersion relation in `dispersion.py` carries its own implementation, checked against trapezoidal quadrature of the integral representation rather than against the same series; using a second, independent implementation here means that agreement between the measured profile and the analytic one cannot be an artefact shared by both.
    """
    from scipy.special import iv

    return iv(order, x)


def analytic_profiles(r, k):
    """Inviscid eigenfunction, normalised to unit maximum on the interval."""
    u_r = bessel_i(1, k * r)
    u_z = bessel_i(0, k * r)
    pressure = bessel_i(0, k * r)

    def unit(a):
        peak = np.max(np.abs(a))
        return a / peak if peak > 0 else a

    return unit(u_r), unit(u_z), unit(pressure)


def radial_profile(mode, shape, wavelength, z, component):
    """Radial profile of one component, extracted from a DMD mode.

    The mode is a full field on the (r,z) plane. Its axial dependence is the imposed cos(k z) or sin(k z), so projecting onto the first axial harmonic collapses the plane to a radial profile and simultaneously rejects anything at another axial wavenumber, which is the numerical content the figure is not about.
    """
    field = mode.reshape(shape)[:, :, 0] if mode.ndim == 1 else mode
    k = 2.0 * np.pi / wavelength
    # Projection onto exp(-i k z); the factor of two makes the coefficient the amplitude of a cosine rather than half of it.
    kernel = np.exp(-1j * k * z)
    return 2.0 * (field * kernel[None, :]).mean(axis=1)


def normalise(profile):
    """Remove the arbitrary complex scale of a DMD mode.

    A DMD eigenvector is defined up to a complex multiple, so amplitude and phase carry no information and must be fixed before any comparison. The scale is set by the largest entry, which is a choice of units; the phase is set by rotating that entry onto the positive real axis, which is a choice of origin in z. Neither can alter the shape of the profile, which is the quantity under test.
    """
    profile = np.asarray(profile)
    peak = int(np.argmax(np.abs(profile)))
    if np.abs(profile[peak]) == 0.0:
        return profile.real
    rotated = profile * np.exp(-1j * np.angle(profile[peak]))
    return (rotated / np.abs(profile[peak])).real


def figure_mode_shape(r, measured, analytic, labels, r_interface, outdir):
    """Measured radial profiles against the analytic eigenfunction.

    The analytic solution is drawn as a line and the measurement as open circles on top of it, following the house convention, so that agreement shows as the line passing through the markers rather than being asserted in a caption.
    """
    n = len(labels)
    fig, axes = plt.subplots(1, n, figsize=(3.4 * n, 3.8), sharex=True)
    axes = np.atleast_1d(axes)

    for ax, meas, exact, label in zip(axes, measured, analytic, labels):
        ax.plot(r, exact, "-", label="analytic, inviscid")
        ax.plot(r, meas, "o", ms=3.4, mfc="none", color="k",
                label="extracted mode")
        ax.axvline(r_interface, ls=":", lw=GUIDE_LW, color="0.45")
        ax.set_xlabel("$r$")
        ax.set_title(label)
        ax.set_xlim(float(r.min()), float(r.max()))
        cm_both(ax)
    axes[0].set_ylabel("normalised amplitude")
    axes[0].legend(loc="upper left")
    fig.tight_layout()
    savefig_all(fig, "mode_shape", outdir=outdir)
    plt.close(fig)


def figure_mode_error(r, measured, analytic, labels, r_interface, outdir):
    """Departure of the measured profile from the analytic one.

    Drawn separately because on the comparison panels the two curves overlie each other and a departure of a few per cent is invisible. The expectation is that the departure is small through the column and largest near the interface, where the viscous solution satisfies a stress condition the inviscid one does not; a figure showing the opposite would indicate an error in the extraction rather than a viscous effect.
    """
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    for meas, exact, label in zip(measured, analytic, labels):
        ax.plot(r, 100.0 * (meas - exact), label=label)
    ax.axvline(r_interface, ls=":", lw=GUIDE_LW, color="0.45",
               label="interface")
    ax.axhline(0.0, lw=GUIDE_LW, color="0.7")
    ax.set_xlabel("$r$")
    ax.set_ylabel("departure from the analytic profile (per cent)")
    ax.set_xlim(float(r.min()), float(r.max()))
    cm_both(ax)
    ax.legend()
    fig.tight_layout()
    savefig_all(fig, "mode_shape_error", outdir=outdir)
    plt.close(fig)


def compute(args):
    """Extract the leading mode and store its radial profiles."""
    import dmd as dmd_mod
    import nsmfp

    fields = ("x_velocity", "y_velocity")
    print("assembling the snapshot matrix ...", flush=True)
    snapshots, dt, times = nsmfp.build_snapshot_matrix(
        args.perturbed, args.base, fields=fields,
        skip=args.skip, stride=args.stride)

    result = dmd_mod.dmd(snapshots, dt, rank=args.rank)
    lead = int(np.argmax(result.amplitude))
    mode = result.modes[:, lead]
    print(f"    leading sigma {result.growth_rate[lead]:.6f}", flush=True)

    # The state vector concatenates the fields, so it is split before either half is reshaped onto the plane.
    sample = nsmfp.read_plotfile(
        nsmfp.find_plotfiles(args.perturbed, prefix=args.prefix)[0],
        fields=(fields[0],))
    shape = sample.shape
    per_field = int(np.prod(shape))
    n_r, n_z = shape[0], shape[1]

    r = (np.arange(n_r) + 0.5) * float(sample.domain_right[0]) / n_r
    z = (np.arange(n_z) + 0.5) * float(sample.domain_right[1]) / n_z
    wavelength = float(sample.domain_right[1])

    u_r = radial_profile(mode[:per_field].reshape(shape)[:, :, 0],
                         shape, wavelength, z, "u_r")
    u_z = radial_profile(mode[per_field:2 * per_field].reshape(shape)[:, :, 0],
                         shape, wavelength, z, "u_z")

    save_data(DATA, r=r, u_r=normalise(u_r), u_z=normalise(u_z),
              meta={"wavelength": wavelength, "kr0": args.kr0,
                    "r_interface": args.r0,
                    "sigma": float(result.growth_rate[lead])})
    return DATA


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--perturbed")
    ap.add_argument("--base")
    ap.add_argument("--replot", action="store_true")
    ap.add_argument("--skip", type=int, default=60)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--rank", type=int, default=None)
    ap.add_argument("--kr0", type=float, default=0.7)
    ap.add_argument("--r0", type=float, default=1.0)
    ap.add_argument("--prefix", default="nddataPLT")
    ap.add_argument("--r-max", type=float, default=2.0, dest="r_max",
                    help="outer radius shown; the far field carries no mode")
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
    r = d["r"]
    kr0 = float(d["meta"].get("kr0", args.kr0))
    r_interface = float(d["meta"].get("r_interface", args.r0))
    k = kr0 / r_interface

    # The mode decays into the far field, where both the measurement and the analytic form are near zero and the normalised comparison is dominated by round-off. The figure is therefore confined to the neighbourhood of the column.
    inside = r <= args.r_max
    r = r[inside]
    measured = [d["u_r"][inside], d["u_z"][inside]]

    exact_r, exact_z, _ = analytic_profiles(r, k)
    # The analytic profiles are normalised over the plotted interval, matching how the measured ones were normalised, so the two are compared on the same footing.
    analytic = [exact_r / np.max(np.abs(exact_r)),
                exact_z / np.max(np.abs(exact_z))]
    measured = [m / np.max(np.abs(m)) if np.max(np.abs(m)) > 0 else m
                for m in measured]
    # A DMD eigenvector carries an arbitrary overall sign, which normalisation by the largest entry does not fix. It is resolved by correlation against the analytic profile, which chooses between two possibilities and cannot alter the shape.
    measured = [m if np.dot(m, a) >= 0 else -m
                for m, a in zip(measured, analytic)]

    labels = [r"radial velocity $u_r$", r"axial velocity $u_z$"]
    figure_mode_shape(r, measured, analytic, labels, r_interface, args.outdir)
    figure_mode_error(r, measured, analytic, labels, r_interface, args.outdir)

    for label, m, a in zip(labels, measured, analytic):
        rms = float(np.sqrt(np.mean((m - a) ** 2)))
        print(f"{label:28s} r.m.s. departure {100 * rms:.2f} per cent")
    print(f"figures written to {args.outdir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
