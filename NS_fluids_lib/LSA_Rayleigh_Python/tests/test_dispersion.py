"""Validate dispersion.py against the algorithm in rayleigh_capillary_growth.F90.

The repo utility computes I_n(x) with its own truncated power series rather than a library call, so we re-implement that series literally here and check (a) it agrees with scipy.special.iv, and (b) the assembled growth rate agrees with dispersion.growth_rate.  That pins the Python port to the Fortran the solver ships with, without needing a Fortran compiler.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

import dispersion  # noqa: E402


def mod_bessel_first_kind(x, n):
    """Literal transcription of `mod_bessel_first_kind` from the .F90 file.

        m_fact=1; mpn_fact=n!
        xpower=(x/2)^n
        I_n = sum_{i=0}^{100} xpower/(m_fact*mpn_fact)
        with xpower *= x^2/4, m_fact *= i, mpn_fact *= (i+n) each iteration.
    """
    m_fact = 1.0
    mpn_fact = 1.0
    for i in range(1, n + 1):
        mpn_fact *= i
    xpower = (0.5 * x) ** n
    I_n = 0.0
    for i in range(0, 101):
        if i > 0:
            m_fact *= i
            mpn_fact *= i + n
            xpower *= x * x / 4.0
        I_n += xpower / (m_fact * mpn_fact)
    return I_n


def fortran_growth_rate(k, density=1.0, r0=1.0, sigma=1.0):
    """Literal transcription of the growth-rate line in the .F90 program:

        c2 = sigma*I_1*x*(1-x*x)/(density*(r0**3)*I_0);  c = sqrt(c2)
    """
    x = k * r0
    I0 = mod_bessel_first_kind(x, 0)
    I1 = mod_bessel_first_kind(x, 1)
    c2 = sigma * I1 * x * (1.0 - x * x) / (density * (r0**3) * I0)
    return np.sqrt(c2) if c2 >= 0.0 else float("nan")


@pytest.mark.parametrize("x", [0.05, 0.1, 0.3, 0.5, 0.697, 0.7, 0.9, 0.99])
def test_bessel_series_matches_scipy(x):
    from scipy.special import iv

    assert mod_bessel_first_kind(x, 0) == pytest.approx(iv(0, x), rel=1e-13)
    assert mod_bessel_first_kind(x, 1) == pytest.approx(iv(1, x), rel=1e-13)


@pytest.mark.parametrize("k", [0.01, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 0.95])
def test_growth_rate_matches_fortran(k):
    """Python port must reproduce the Fortran formula to machine precision."""
    assert float(dispersion.growth_rate(k)) == pytest.approx(
        fortran_growth_rate(k), rel=1e-12, abs=1e-14
    )


def test_fortran_header_claim_k_r0_0p7():
    """The .F90 header asserts 'anecdotedly, k r0 = .7 is the critical value'
    for the most dangerous mode.  Confirm the true maximiser is ~0.697."""
    k_max, _ = dispersion.most_unstable_wavenumber(r0=1.0)
    assert k_max == pytest.approx(0.697, abs=0.005)


def test_wavelength_of_most_unstable_mode():
    """inputs.growthrate.LSA uses yblob = 9.0 for 'wave length=2 pi r0/.7'."""
    lam = dispersion.most_unstable_wavelength(r0=1.0)
    assert lam == pytest.approx(2.0 * np.pi / 0.7, rel=0.01)
    assert lam == pytest.approx(9.0, rel=0.01)


def test_stability_boundary():
    """Unstable strictly for k r0 < 1; marginal at k r0 = 1."""
    assert dispersion.growth_rate_squared(0.999) > 0.0
    assert dispersion.growth_rate_squared(1.0) == pytest.approx(0.0, abs=1e-12)
    assert dispersion.growth_rate_squared(1.001) < 0.0


def test_scaling_with_physical_parameters():
    """sigma_g scales as sqrt(sigma/(rho r0^3)) at fixed k r0."""
    k, r0 = 0.7, 1.0
    base = dispersion.growth_rate(k, r0=r0, sigma=1.0, rho=1.0)
    quad = dispersion.growth_rate(k, r0=r0, sigma=4.0, rho=1.0)
    assert quad == pytest.approx(2.0 * base, rel=1e-12)
    heavy = dispersion.growth_rate(k, r0=r0, sigma=1.0, rho=4.0)
    assert heavy == pytest.approx(0.5 * base, rel=1e-12)
    # fixed k*r0 = 0.7 while doubling r0 -> sigma_g scales as r0^{-3/2}
    big = dispersion.growth_rate(0.35, r0=2.0, sigma=1.0, rho=1.0)
    assert big == pytest.approx(base * 2.0 ** (-1.5), rel=1e-12)


def test_viscous_reduces_growth_and_limits_to_inviscid():
    k = 0.7
    inviscid = dispersion.growth_rate(k)
    assert dispersion.growth_rate_viscous(k, mu=0.0) == pytest.approx(inviscid)
    visc = dispersion.growth_rate_viscous(k, mu=0.02)
    assert 0.0 < visc < inviscid


def test_admissible_wavenumbers_in_periodic_box():
    """L = 72, r0 = 1 (the default deck) should admit many unstable modes,
    all with k r0 < 1, and n = 8 should be closest to the peak (72/9 = 8)."""
    n, k = dispersion.admissible_wavenumbers(L=72.0, r0=1.0)
    assert np.all(k < 1.0)
    k_max, _ = dispersion.most_unstable_wavenumber()
    n_best = n[int(np.argmin(np.abs(k - k_max)))]
    assert n_best == 8


def test_box_length_for_mode_round_trip():
    k_max, _ = dispersion.most_unstable_wavenumber()
    L = dispersion.box_length_for_mode(k_max, n=1)
    assert L == pytest.approx(9.0, rel=0.01)
