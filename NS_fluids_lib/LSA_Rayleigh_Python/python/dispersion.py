"""
Exact linear dispersion relation for the Rayleigh-Plateau capillary instability.

This is a Python port / generalisation of

    NS_fluids_lib/run2d/rayleigh_capillary_growth.F90

which implements the inviscid Rayleigh (1878) result quoted in Popinet (2009), section 6.5.  Keeping the two implementations side by side lets us verify that the Python reference used for validation reproduces the Fortran utility that ships with the solver (see tests/test_dispersion.py).

Inviscid dispersion relation for an axisymmetric liquid column of radius r0, density rho, surface tension sigma, perturbed as

    r(z,t) = r0 + eps * exp(sigma_g t) * cos(k z)

is

    sigma_g^2 = (sigma / (rho r0^3)) * x * I1(x) / I0(x) * (1 - x^2),
    x = k r0.

The column is unstable (sigma_g^2 > 0) for x < 1, i.e. for disturbance wavelengths longer than the circumference 2*pi*r0.  The most-amplified wavenumber is x ~= 0.697 ("k r0 = .7 is the critical value" in the Fortran header comment), giving lambda_max ~= 2*pi*r0/0.697 ~= 9.02 r0.

Viscosity is *not* included here.  The AMReX runs are viscous, so the measured growth rate is expected to sit below this inviscid envelope; the gap grows with the Ohnesorge number Oh = mu / sqrt(rho sigma r0).  `growth_rate_viscous` provides the standard long-wave viscous correction for a quantitative comparison at finite Oh.
"""

# ORIENTATION
#
# This module contains theory, not measurement. Nothing here reads a calculation; every function evaluates a closed-form or root-found prediction that a measured growth rate is compared against.
#
# The physical situation is a cylinder of liquid of radius r0. Surface tension makes it unstable to long-wavelength varicose disturbances, because a bulge has lower surface energy per unit volume than a neck: this is the Rayleigh-Plateau instability, and it is why a tap stream breaks into drops. A disturbance of axial wavenumber k grows as exp(sigma t), and sigma(k) is what these functions return.
#
# Three levels of description are provided, in increasing fidelity and cost.
#
#   growth_rate                 inviscid, closed form, one Bessel ratio
#   growth_rate_viscous         viscous, an explicit quadratic approximation
#   growth_rate_viscous_exact   viscous, exact, found by root-finding
#
# Read `growth_rate` first: everything else is a correction to it. The wavenumber enters only as the dimensionless product x = k*r0, and the instability exists only for x < 1, that is for wavelengths longer than the circumference. The fastest-growing disturbance is at x = 0.697, which is why a breaking jet produces drops of a fairly definite size.
#
# A note on units. Every function accepts dimensional arguments and returns a dimensional rate, but works internally in units where lengths are scaled by r0 and times by the capillary time sqrt(rho*r0^3/gamma). In those units the kinematic viscosity is numerically equal to the Ohnesorge number Oh, which is the single parameter measuring viscous against capillary effects. If you write a function like these yourself, non-dimensionalise first; it removes three parameters and makes the limits checkable.


from __future__ import annotations

import numpy as np
from scipy.special import iv  # modified Bessel function of the first kind


def growth_rate(k, r0=1.0, sigma=1.0, rho=1.0):
    """Inviscid Rayleigh-Plateau growth rate sigma_g(k).

    Parameters
    ----------
    k : float or array_like
        Axial wavenumber (rad / length).
    r0, sigma, rho : float
        Unperturbed column radius, surface tension coefficient, liquid density.

    Returns
    -------
    ndarray
        Growth rate.  Zero is returned in the stable band (x >= 1) rather than
        the imaginary oscillation frequency, so the result is always real; use
        `growth_rate_squared` if the stable branch is needed.
    """
    s2 = growth_rate_squared(k, r0=r0, sigma=sigma, rho=rho)
    return np.sqrt(np.maximum(s2, 0.0))


def growth_rate_squared(k, r0=1.0, sigma=1.0, rho=1.0):
    """sigma_g^2(k); negative in the stable band x > 1 (oscillatory)."""
    k = np.asarray(k, dtype=float)
    x = k * r0
    # x -> 0 limit: I1(x)/I0(x) ~ x/2, so x*I1/I0 ~ x^2/2 -> 0.  np.where alone would still evaluate the ratio at x = 0, which is 0/1 = 0 and therefore safe, but we guard anyway to keep the expression clean for x = 0 exactly.
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(x > 0.0, iv(1, x) / iv(0, x), 0.0)
    return (sigma / (rho * r0**3)) * x * ratio * (1.0 - x**2)


def most_unstable_wavenumber(r0=1.0, sigma=1.0, rho=1.0, n=200001):
    """Return (k_max, sigma_g_max) for the fastest-growing mode.

    Found by dense sampling plus a parabolic refinement, which is ample here because the curve is smooth and single-peaked on 0 < x < 1.
    """
    x = np.linspace(1e-9, 1.0, n)
    k = x / r0
    s = growth_rate(k, r0=r0, sigma=sigma, rho=rho)
    i = int(np.argmax(s))
    if 0 < i < n - 1:
        # parabolic vertex through the three points bracketing the maximum
        y0, y1, y2 = s[i - 1], s[i], s[i + 1]
        denom = y0 - 2.0 * y1 + y2
        delta = 0.5 * (y0 - y2) / denom if denom != 0.0 else 0.0
        dk = k[1] - k[0]
        k_max = k[i] + delta * dk
        return k_max, float(growth_rate(k_max, r0=r0, sigma=sigma, rho=rho))
    return float(k[i]), float(s[i])


def most_unstable_wavelength(r0=1.0, sigma=1.0, rho=1.0):
    """Wavelength of the fastest-growing mode, lambda = 2 pi / k_max."""
    k_max, _ = most_unstable_wavenumber(r0=r0, sigma=sigma, rho=rho)
    return 2.0 * np.pi / k_max


def ohnesorge(mu, rho=1.0, sigma=1.0, r0=1.0):
    """Oh = mu / sqrt(rho sigma r0).  Viscous effects matter once Oh ~ 1."""
    return mu / np.sqrt(rho * sigma * r0)


def growth_rate_viscous(k, r0=1.0, sigma=1.0, rho=1.0, mu=0.0):
    """Growth rate including the leading viscous correction.

    Uses the standard quadratic approximation to the viscous Rayleigh problem (e.g. Chandrasekhar, *Hydrodynamic and Hydromagnetic Stability*, ch. XII),

        s^2 + 3 nu k^2 s - sigma_inviscid^2 = 0,

    solved for the positive root

        s = -1.5 nu k^2 + sqrt( (1.5 nu k^2)^2 + sigma_inviscid^2 ).

    This reduces to the inviscid result as mu -> 0 and captures the correct Oh >> 1 asymptote s -> sigma_inviscid^2 / (3 nu k^2).  It is an approximation, not the exact viscous eigenvalue, and is provided so the measured DMD growth rates can be compared against a viscous reference at the modest Oh of the default input decks rather than only against the inviscid envelope.
    """
    nu = mu / rho
    s2_inv = growth_rate_squared(k, r0=r0, sigma=sigma, rho=rho)
    a = 1.5 * nu * np.asarray(k, dtype=float) ** 2
    # In the stable band s2_inv < 0 and the discriminant can go negative, meaning the root is complex: a decaying *oscillation* rather than a decaying non-oscillatory mode.  We report only the real growth part, so clamp at zero rather than returning NaN.
    disc = np.maximum(a**2 + s2_inv, 0.0)
    return -a + np.sqrt(disc)


def admissible_wavenumbers(L, r0=1.0, n_max=None):
    """Wavenumbers that fit an axially periodic box of length L.

    Only k_n = 2 pi n / L are representable.  Returns the unstable ones (k_n r0 < 1), which are the only modes the box can grow.
    """
    if n_max is None:
        n_max = int(np.floor(L / (2.0 * np.pi * r0))) + 2
    n = np.arange(1, max(n_max, 1) + 1)
    k = 2.0 * np.pi * n / L
    keep = k * r0 < 1.0
    return n[keep], k[keep]


def box_length_for_mode(k, r0=1.0, n=1):
    """Axial box length that holds exactly n wavelengths of wavenumber k."""
    return 2.0 * np.pi * n / k


if __name__ == "__main__":
    # Reproduces the table printed by rayleigh_capillary_growth.F90 (rho = r0 = sigma = 1), for eyeball comparison against the Fortran.
    k_max, s_max = most_unstable_wavenumber()
    print(f"# most unstable: k r0 = {k_max:.6f}  sigma_g = {s_max:.6f}  "
          f"lambda/r0 = {2*np.pi/k_max:.4f}")
    print(f"{'k':>10} {'sigma_g':>14}")
    for i in range(0, 101):
        k = i / 100.0
        print(f"{k:10.4f} {float(growth_rate(k)):14.8f}")


# ----------------------------------------------------------------------
# Exact references
#
# The approximation above is the one this study originally compared against, and it remains useful because it is explicit. It is not, however, the exact viscous eigenvalue, and a validation claim quoted to a fraction of a per cent cannot rest on a formula whose own error at the relevant Ohnesorge number is unquantified. The two functions below supply exact references: the Chandrasekhar viscous relation for a column in vacuum, and the inviscid two-fluid relation, which quantifies the ambient fluid the calculation actually contains and the single-fluid theory omits.
# ----------------------------------------------------------------------


def _bessel_ratio_i1p_i1(z):
    """I_1'(z) / I_1(z), evaluated stably at large argument.

    The root-find below drives the argument to large values as the viscosity falls, where I_0 and I_1 both overflow while their ratio stays near unity. The exponentially scaled forms carry the same factor in numerator and denominator, so it cancels and the ratio is obtained without ever forming the overflowing quantity.
    """
    from scipy.special import ive

    return ive(0, z) / ive(1, z) - 1.0 / z


def _viscous_residual(s, x, oh):
    """Residual of the exact viscous dispersion relation.

    Chandrasekhar, *Hydrodynamic and Hydromagnetic Stability*, chapter XII, in the non-dimensionalisation where lengths are scaled by the column radius and times by the capillary time sqrt(rho r0^3 / gamma), in which the kinematic viscosity equals the Ohnesorge number:

        (l^2 - x^2)/(l^2 + x^2) . x(1 - x^2) I_1(x)/I_0(x)
            = s^2 + 2 Oh x^2 s [ I_1'(x)/I_0(x)
                                 - 2 x l/(x^2 + l^2) . I_1(x) I_1'(l)
                                   / (I_0(x) I_1(l)) ]

    with l^2 = x^2 + s/Oh. As Oh tends to zero, l diverges, the prefactor on the left tends to unity and the bracket on the right stays bounded while its coefficient vanishes, so the relation reduces to the Rayleigh result. That limit is checked in the test suite rather than asserted here.
    """
    from scipy.special import iv

    l = np.sqrt(x * x + s / oh)
    i0x, i1x = iv(0, x), iv(1, x)
    i1p_over_i0 = (_bessel_ratio_i1p_i1(x) * i1x) / i0x
    bracket = (i1p_over_i0
               - (2.0 * x * l / (x * x + l * l)) * (i1x / i0x)
               * _bessel_ratio_i1p_i1(l))
    lhs = s * s + 2.0 * oh * x * x * s * bracket
    rhs = ((l * l - x * x) / (l * l + x * x)) * x * (1.0 - x * x) * i1x / i0x
    return lhs - rhs


def growth_rate_viscous_exact(k, r0=1.0, sigma=1.0, rho=1.0, mu=0.0):
    """Exact viscous growth rate, by root-finding on the dispersion relation.

    Returned in the same units as `growth_rate`, so it is directly comparable with a measured rate and with the approximation. The inviscid result bounds the root from above, since viscosity can only retard a capillary instability, and zero bounds it from below in the unstable band; the bracket is therefore guaranteed and no initial guess is required.

    In the stable band, x > 1, there is no growing root and zero is returned, matching the convention of the approximation.
    """
    from scipy.optimize import brentq

    k = np.asarray(k, dtype=float)
    scalar = k.ndim == 0
    k = np.atleast_1d(k)

    # The capillary time and the Ohnesorge number set the dimensionless problem, and the answer is scaled back on the way out.
    t_c = np.sqrt(rho * r0**3 / sigma)
    oh = ohnesorge(mu, rho=rho, sigma=sigma, r0=r0)

    out = np.zeros_like(k)
    for index, k_value in enumerate(k.ravel()):
        x = float(k_value * r0)
        if x <= 0.0 or x >= 1.0:
            continue
        upper = float(np.sqrt(x * (1.0 - x * x) * _i1_over_i0(x)))
        if oh <= 0.0:
            out.ravel()[index] = upper / t_c
            continue
        lo, hi = 1e-14, upper
        if _viscous_residual(lo, x, oh) * _viscous_residual(hi, x, oh) > 0.0:
            # No sign change means the inviscid bound is not above the root, which happens only if the relation has been evaluated outside its range of validity. Reporting NaN is preferable to returning an endpoint that would silently enter a figure.
            out.ravel()[index] = np.nan
            continue
        root = brentq(_viscous_residual, lo, hi, args=(x, oh),
                      xtol=1e-15, rtol=8.9e-16, maxiter=200)
        out.ravel()[index] = root / t_c

    return float(out[0]) if scalar else out


def _i1_over_i0(x):
    """I_1(x)/I_0(x), the combination that recurs in the relations above."""
    from scipy.special import iv

    return iv(1, x) / iv(0, x)


def growth_rate_two_fluid_inviscid(k, r0=1.0, sigma=1.0, rho=1.0,
                                   rho_outer=0.0):
    """Inviscid growth rate with an ambient fluid of finite density.

    The calculations of this study are two-fluid: the deck carries an ambient phase of density 1.225e-3 and viscosity 2.6e-4 against a column of density 1 and viscosity 2e-2. The single-fluid relation omits it entirely, and the size of that omission should be measured rather than assumed negligible.

    Matching potential flow inside and outside the column, with the disturbance decaying as K_0(kr) in the ambient phase, gives

        s^2 = x(1 - x^2) / [ I_0(x)/I_1(x) + (rho_2/rho_1) K_0(x)/K_1(x) ]

    which reduces to the Rayleigh result when the ambient density vanishes. The ambient fluid adds inertia and therefore always reduces the growth rate.
    """
    from scipy.special import iv, kv

    k = np.asarray(k, dtype=float)
    x = k * r0
    t_c = np.sqrt(rho * r0**3 / sigma)
    inside = iv(0, x) / iv(1, x)
    outside = (rho_outer / rho) * kv(0, x) / kv(1, x)
    s2 = np.where(x < 1.0, x * (1.0 - x * x) / (inside + outside), 0.0)
    return np.sqrt(np.maximum(s2, 0.0)) / t_c
