"""Tests for the exact viscous and two-fluid dispersion relations.

A reference curve that a validation claim is quoted against must itself be checked, and checked against something other than the formula it is meant to replace. These tests pin the exact viscous relation by its limits, where the answer is known independently: it must reduce to the Rayleigh result as viscosity vanishes, it must lie below the inviscid rate at every finite viscosity since viscosity can only retard a capillary instability, it must decrease monotonically with Ohnesorge number, and its most unstable wavenumber must drift towards longer waves. None of those follows from the implementation, so each is evidence that the implementation is right.

The relation is also required to agree with the quadratic approximation in the limit where the approximation is itself asymptotically correct, and to depart from it as viscosity rises. Agreement everywhere would mean the exact relation had not been implemented at all.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "python"))

import dispersion  # noqa: E402


def test_reduces_to_rayleigh_as_viscosity_vanishes():
    for kr0 in (0.2, 0.45, 0.697, 0.9):
        exact = dispersion.growth_rate_viscous_exact(kr0, mu=0.0)
        assert exact == pytest.approx(dispersion.growth_rate(kr0), rel=1e-12)


def test_approaches_rayleigh_continuously():
    """The limit must be approached, not merely attained at exactly zero."""
    kr0 = 0.7
    inviscid = dispersion.growth_rate(kr0)
    previous = 0.0
    for oh in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
        exact = dispersion.growth_rate_viscous_exact(kr0, mu=oh)
        gap = abs(exact - inviscid)
        if previous:
            assert gap < previous
        previous = gap
    assert previous < 1e-5


def test_viscosity_only_retards():
    for kr0 in (0.3, 0.5, 0.7, 0.9):
        inviscid = dispersion.growth_rate(kr0)
        for oh in (0.005, 0.02, 0.1, 0.5):
            assert dispersion.growth_rate_viscous_exact(kr0, mu=oh) < inviscid


def test_monotone_in_ohnesorge():
    rates = [dispersion.growth_rate_viscous_exact(0.7, mu=oh)
             for oh in (0.001, 0.01, 0.02, 0.05, 0.1, 0.3)]
    assert all(b < a for a, b in zip(rates, rates[1:]))


def test_peak_drifts_to_longer_waves():
    """Viscosity shifts the most unstable wavenumber below the inviscid 0.697."""
    grid = np.linspace(0.05, 0.99, 400)
    peaks = []
    for oh in (1e-6, 0.05, 0.2, 0.5):
        rates = dispersion.growth_rate_viscous_exact(grid, mu=oh)
        peaks.append(grid[int(np.nanargmax(rates))])
    assert peaks[0] == pytest.approx(0.697, abs=0.005)
    assert all(b <= a for a, b in zip(peaks, peaks[1:]))
    assert peaks[-1] < peaks[0]


def test_agrees_with_the_approximation_where_it_is_valid():
    """The approximation is asymptotically correct as Oh vanishes.

    The test is that the discrepancy shrinks with viscosity, not that it falls below a fixed tolerance. A fixed tolerance would either pass trivially or encode the discrepancy at one particular Ohnesorge number; the vanishing of the difference in the limit is the property that identifies the two expressions as descriptions of the same problem.
    """
    gaps = []
    for oh in (1e-3, 1e-4, 1e-5, 1e-6):
        exact = dispersion.growth_rate_viscous_exact(0.7, mu=oh)
        approx = dispersion.growth_rate_viscous(0.7, mu=oh)
        gaps.append(abs(exact - approx) / exact)
    assert all(b < a for a, b in zip(gaps, gaps[1:]))
    assert gaps[0] < 1e-4


def test_departs_from_the_approximation_as_viscosity_rises():
    """And is not merely a re-expression of it.

    At the Ohnesorge number of this study the two differ by about 0.13 per cent, which is small but is the same order as the precision the measurement is quoted to, and is the reason the exact relation is needed.
    """
    exact = dispersion.growth_rate_viscous_exact(0.7, mu=0.02)
    approx = dispersion.growth_rate_viscous(0.7, mu=0.02)
    assert exact != pytest.approx(approx, rel=1e-4)
    assert exact > approx
    assert abs(exact - approx) / exact < 0.01

    far = dispersion.growth_rate_viscous_exact(0.7, mu=0.5)
    far_approx = dispersion.growth_rate_viscous(0.7, mu=0.5)
    assert abs(far - far_approx) / far > abs(exact - approx) / exact


def test_stable_band_returns_zero():
    for kr0 in (1.0, 1.1, 2.0):
        assert dispersion.growth_rate_viscous_exact(kr0, mu=0.02) == 0.0


def test_accepts_an_array():
    grid = np.array([0.2, 0.5, 0.7])
    out = dispersion.growth_rate_viscous_exact(grid, mu=0.02)
    assert out.shape == grid.shape
    for kr0, value in zip(grid, out):
        assert value == pytest.approx(
            dispersion.growth_rate_viscous_exact(float(kr0), mu=0.02))


def test_two_fluid_reduces_to_single_fluid():
    for kr0 in (0.2, 0.5, 0.7):
        assert dispersion.growth_rate_two_fluid_inviscid(
            kr0, rho_outer=0.0) == pytest.approx(dispersion.growth_rate(kr0))


def test_ambient_fluid_retards_and_is_negligible_at_the_deck_ratio():
    """The ambient phase cannot account for a discrepancy of per cent order.

    The deck carries an ambient density of 1.225e-3 against a column of 1. The added inertia reduces the growth rate, but by a few hundredths of a per cent, so the single-fluid reference is not the source of the departure between theory and measurement.
    """
    kr0 = 0.7
    single = dispersion.growth_rate(kr0)
    both = dispersion.growth_rate_two_fluid_inviscid(kr0, rho_outer=1.225e-3)
    assert both < single
    assert abs(both - single) / single < 1e-3

    # A heavy ambient phase must produce a large effect, so that the small result above reflects the density ratio and not an inert parameter.
    heavy = dispersion.growth_rate_two_fluid_inviscid(kr0, rho_outer=1.0)
    assert abs(heavy - single) / single > 0.05
