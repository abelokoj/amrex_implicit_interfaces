"""
Interface-mode amplitude: the correct observable for an interfacial instability.

WHY THIS MODULE EXISTS (learned from running the solver, not from theory)
------------------------------------------------------------------------
The obvious NS-MFP observable is the norm of the perturbation field,
|Q'(t)| = |Q_pert(t) - Q_base(t)|.  For the Rayleigh-Plateau problem that
choice fails in two different ways depending on which field you pick, and both failures were only visible once real solver output was available.

1. VELOCITY FIELDS give Q'(0) = 0 exactly.
The perturbation lives on the interface, so both twins start from rest and the first snapshot of the perturbation field is identically zero.  That is fatal for the DMD amplitude ranking, which follows Ranjan et al. in projecting each mode onto the FIRST snapshot: every amplitude comes out zero and the mode ordering becomes meaningless.  Velocities are still a fine state vector, but only after the startup transient is skipped.

2. THE RAW LEVEL SET dilutes the signal.
L0101 is a signed distance function, so
       Q' = LS_pert - LS_base ~ eps*sin(k z)
across the WHOLE domain, not just near the interface.  Only the near-interface part participates in the instability; the rest is a static offset.  In a real run the raw level-set norm stayed flat to 0.6% while the interface amplitude genuinely grew by 37%.  A flat norm therefore does NOT mean a stable interface - it can mean the observable is wrong.

The fix is to measure the thing the dispersion relation actually predicts: the amplitude of the k-Fourier component of the interface radius r(z, t). For r(z,t) = r0 + a(t) cos(k z) linear theory gives a(t) = a(0) e^{sigma t}, so a simple log-linear fit of a(t) returns sigma directly, immune to both the zero-initial-condition problem and the static-offset dilution.

This is used alongside the DMD, not instead of it: DMD still supplies the spectrum and the mode shapes, while this gives a robust scalar growth rate to check it against.
"""

# ORIENTATION
#
# This module measures the shape of the interface, which is the quantity the growth rate is fitted to.
#
# The solver represents the interface implicitly, as the zero contour of a level-set field: the field is negative inside the liquid, positive outside, and the interface is wherever it crosses zero. To find the interface radius at one axial station, walk outwards along that row of cells until the sign changes, then interpolate linearly between the two straddling cell centres to locate the crossing. Doing that at every station gives the profile r(z).
#
# The measured quantity is then the amplitude of the first Fourier harmonic of r(z), obtained by projecting the profile onto cos(kz). This is the right observable for a specific reason. The velocity difference between the twins is exactly zero at the first instant, because both runs start from rest, so a velocity-based measure has no signal to begin from. The raw level-set difference is worse still, since it spreads the discrepancy over the whole domain and dilutes it. The interface amplitude is what the theory actually predicts, and it is what is measured here.


from __future__ import annotations

import numpy as np

import nsmfp


def interface_radius(plotfile_dir, level_set_field="L0101", level=0,
                     r_max=None):
    """Interface radius r(z) from the zero contour of the level set.

    Linear interpolation of the first sign change in each axial column. Returns (time, z, radius).  Columns with no sign change (interface out of the domain) get NaN rather than a silently wrong value.
    """
    pf = nsmfp.read_plotfile(plotfile_dir, fields=(level_set_field,),
                             level=level)
    ls = pf.data.reshape(pf.shape)[:, :, 0]           # (nr, nz)
    nr, nz = ls.shape

    if r_max is None:
        r_max = float(pf.domain_right[0])
    z_max = float(pf.domain_right[1])

    r = (np.arange(nr) + 0.5) * r_max / nr
    z = (np.arange(nz) + 0.5) * z_max / nz

    radius = np.full(nz, np.nan)
    for j in range(nz):
        col = ls[:, j]
        idx = np.where(np.diff(np.sign(col)))[0]
        if idx.size:
            i = idx[0]
            denom = col[i] - col[i + 1]
            f = col[i] / denom if denom != 0 else 0.0
            radius[j] = r[i] + f * (r[i + 1] - r[i])
    return float(pf.time), z, radius


def mode_amplitude(z, radius, wavelength, harmonic=1):
    """Amplitude of the `harmonic`-th Fourier component of r(z).

    a = 2 |<(r - <r>) exp(-i 2 pi n z / lambda)>|

    The factor 2 makes `a` the amplitude of a cos(k z) perturbation, matching the convention of the dispersion relation.  The mean is removed first so that a drifting mean radius (mass conservation error, for instance) does not leak into the mode amplitude.
    """
    m = np.isfinite(radius)
    if m.sum() < 4:
        return float("nan")
    zz, rr = z[m], radius[m]
    rr = rr - rr.mean()
    phase = np.exp(-1j * 2.0 * np.pi * harmonic * zz / wavelength)
    return float(2.0 * np.abs(np.sum(rr * phase) / zz.size))


def amplitude_series(run_dir, wavelength, level_set_field="L0101",
                     harmonic=1, level=0, prefix="nddataPLT"):
    """Interface mode amplitude a(t) over every plotfile in a run.

    Returns (times, amplitudes, mean_radius), each sorted by time. `mean_radius` is returned because a drifting mean is the usual signal of a mass-conservation problem, which would invalidate the growth rate.
    """
    paths = nsmfp.find_plotfiles(run_dir, prefix=prefix)
    if not paths:
        raise FileNotFoundError(
            f"no {prefix}* directories under {run_dir!r}. The solver writes "
            "nddataPLT* when ns.visual_nddata_format=1; pass prefix= if yours "
            "differ.")
    times, amps, means = [], [], []
    for p in paths:
        t, z, rad = interface_radius(p, level_set_field=level_set_field,
                                     level=level)
        times.append(t)
        amps.append(mode_amplitude(z, rad, wavelength, harmonic=harmonic))
        means.append(float(np.nanmean(rad)))
    order = np.argsort(times)
    return (np.asarray(times)[order], np.asarray(amps)[order],
            np.asarray(means)[order])


def fit_growth_rate(times, amps, t_min=None, t_max=None):
    """Log-linear fit of a(t); returns (sigma, intercept, r_squared)."""
    times = np.asarray(times, float)
    amps = np.asarray(amps, float)
    m = np.isfinite(amps) & (amps > 0)
    if t_min is not None:
        m &= times >= t_min
    if t_max is not None:
        m &= times <= t_max
    if m.sum() < 3:
        raise ValueError("need at least 3 positive samples in the fit window")
    t, y = times[m], np.log(amps[m])
    slope, intercept = np.polyfit(t, y, 1)
    pred = slope * t + intercept
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(slope), float(intercept), r2


def local_growth_rate(times, amps, width=0.6, centres=None):
    """Instantaneous growth rate d(ln a)/dt over sliding windows.

    This is the diagnostic that tells you whether the run has reached asymptotic modal growth.  If sigma_local is still CLIMBING at the end of the record, the flow is still in the startup transient and any single fitted growth rate underestimates the true one - which is exactly what the first real runs of this problem showed.
    """
    times = np.asarray(times, float)
    amps = np.asarray(amps, float)
    if centres is None:
        lo, hi = times.min() + width / 2, times.max() - width / 2
        if hi <= lo:
            return np.array([]), np.array([])
        centres = np.linspace(lo, hi, 8)

    out_c, out_s = [], []
    for c in centres:
        m = (times >= c - width / 2) & (times <= c + width / 2)
        m &= np.isfinite(amps) & (amps > 0)
        if m.sum() >= 4:
            s = np.polyfit(times[m], np.log(amps[m]), 1)[0]
            out_c.append(c)
            out_s.append(s)
    return np.asarray(out_c), np.asarray(out_s)


def has_converged(times, amps, width=0.6, rel_tol=0.05):
    """Heuristic: has the local growth rate stopped climbing?

    Compares the last two sliding-window estimates.  Returns (converged, sigma_last, relative_change).  A True here is necessary but not sufficient - also check grid convergence.
    """
    c, s = local_growth_rate(times, amps, width=width)
    if s.size < 2:
        return False, float("nan"), float("nan")
    rel = abs(s[-1] - s[-2]) / max(abs(s[-1]), 1e-30)
    return bool(rel <= rel_tol), float(s[-1]), float(rel)
