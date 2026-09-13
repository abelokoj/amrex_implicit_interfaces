"""Verify the DMD implementation recovers known eigenvalues.

These are synthetic tests with analytically known answers, so a failure here means the DMD is wrong rather than the CFD being hard.  The final test is the one that matters most for this project: a synthetic Rayleigh-Plateau-like signal (a stationary, exponentially growing interface mode buried in decaying noise modes) must yield the planted growth rate.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

import dmd as dmd_mod  # noqa: E402
import dispersion  # noqa: E402


def _snapshots_from_modes(shapes, rates, freqs, amps, dt, n_snap):
    """Build X[:, j] = sum_m amps[m] shapes[:, m] exp((rate + i*2pi*f) j dt)."""
    t = np.arange(n_snap) * dt
    omega = np.asarray(rates) + 1j * 2.0 * np.pi * np.asarray(freqs)
    time = np.exp(np.outer(omega, t))                    # (n_modes, n_snap)
    return shapes @ (np.diag(amps) @ time)


def test_recovers_single_stationary_growing_mode():
    rng = np.random.default_rng(0)
    n_space, n_snap, dt = 200, 60, 0.05
    rate = 0.37
    shape = rng.standard_normal((n_space, 1))
    X = _snapshots_from_modes(shape, [rate], [0.0], [1.0], dt, n_snap)

    res = dmd_mod.dmd(X, dt)
    assert res.leading_growth_rate() == pytest.approx(rate, rel=1e-8)
    assert abs(res.frequency[0]) < 1e-9


def test_recovers_oscillatory_pair():
    """A real oscillating field is a conjugate pair of complex modes.

    The pair must be built from one complex shape and its conjugate: taking the real part then gives a genuinely 2-D spatial subspace (a cos - b sin), which is what a real oscillation looks like.  Two *independent* real shapes would instead collapse to rank 1 and carry no recoverable frequency.
    """
    rng = np.random.default_rng(1)
    n_space, n_snap, dt = 150, 200, 0.02
    rate, freq = -0.15, 1.3

    shape = rng.standard_normal((n_space, 1)) + 1j * rng.standard_normal((n_space, 1))
    shapes = np.hstack([shape, shape.conj()])
    X = _snapshots_from_modes(shapes, [rate, rate], [freq, -freq],
                              [1.0, 1.0], dt, n_snap).real.astype(complex)

    res = dmd_mod.dmd(X, dt)
    got = np.sort(np.abs(res.frequency[:2]))
    assert np.allclose(got, freq, rtol=1e-6)
    assert np.allclose(res.growth_rate[:2], rate, rtol=1e-6)


def test_recovers_multiple_mixed_modes():
    """Growing stationary + decaying oscillatory, as in a real spectrum."""
    rng = np.random.default_rng(2)
    n_space, n_snap, dt = 300, 400, 0.01
    rates = [0.42, -0.30, -0.30, -1.1]
    freqs = [0.0, 2.0, -2.0, 0.0]
    shapes = rng.standard_normal((n_space, 4))
    X = _snapshots_from_modes(shapes, rates, freqs, [1.0, 0.5, 0.5, 0.2],
                              dt, n_snap)

    res = dmd_mod.dmd(X, dt)
    # every planted eigenvalue must appear somewhere in the spectrum
    for r, f in zip(rates, freqs):
        d = np.abs(res.growth_rate - r) + np.abs(res.frequency - f)
        assert d.min() < 1e-6, f"missing mode rate={r} freq={f}"
    # and the leading stationary growth rate is the unstable one
    assert res.leading_growth_rate() == pytest.approx(0.42, rel=1e-6)


def test_growth_rate_sign_convention():
    """Positive sigma must mean growth: |x(T)| > |x(0)|."""
    rng = np.random.default_rng(3)
    dt, n_snap, rate = 0.1, 40, 0.25
    shape = rng.standard_normal((50, 1))
    X = _snapshots_from_modes(shape, [rate], [0.0], [1.0], dt, n_snap)
    assert np.linalg.norm(X[:, -1]) > np.linalg.norm(X[:, 0])
    assert dmd_mod.dmd(X, dt).leading_growth_rate() > 0


def test_robust_to_small_noise():
    rng = np.random.default_rng(4)
    n_space, n_snap, dt, rate = 400, 300, 0.02, 0.31
    shape = rng.standard_normal((n_space, 1))
    X = _snapshots_from_modes(shape, [rate], [0.0], [1.0], dt, n_snap).real
    X = X + 1e-6 * np.linalg.norm(X) / np.sqrt(X.size) * rng.standard_normal(X.shape)

    res = dmd_mod.dmd(X, dt, rank=1)
    assert res.leading_growth_rate() == pytest.approx(rate, rel=1e-3)


def test_rank_truncation_and_tolerance():
    rng = np.random.default_rng(5)
    shapes = rng.standard_normal((120, 3))
    X = _snapshots_from_modes(shapes, [0.2, -0.4, -0.9], [0.0, 1.0, -1.0],
                              [1.0, 1.0, 1.0], 0.05, 100)
    assert dmd_mod.dmd(X, 0.05, rank=2).rank == 2
    # data is exactly rank 3, so tolerance-based selection should find 3
    assert dmd_mod.dmd(X, 0.05, tol=1e-8).rank == 3


def test_invalid_inputs_raise():
    X = np.zeros((10, 5))
    with pytest.raises(ValueError):
        dmd_mod.dmd(X[:, :1], 0.1)          # too few snapshots
    with pytest.raises(ValueError):
        dmd_mod.dmd(X, -1.0)                 # bad dt
    with pytest.raises(ValueError):
        dmd_mod.dmd(np.full((10, 5), np.nan), 0.1)
    with pytest.raises(ValueError):
        dmd_mod.dmd(np.zeros((2, 2, 2)), 0.1)


def test_convergence_study_is_stable_for_clean_data():
    rng = np.random.default_rng(6)
    dt, rate = 0.01, 0.28
    shape = rng.standard_normal((100, 1))
    X = _snapshots_from_modes(shape, [rate], [0.0], [1.0], dt, 500).real
    runs = dmd_mod.convergence_study(X, dt, sizes=[50, 100, 250, 500],
                                     starts=(0,), rank=1)
    assert len(runs) == 4
    for r in runs:
        assert r["growth_rate"] == pytest.approx(rate, rel=1e-6)


def test_synthetic_rayleigh_plateau_growth_rate():
    """End-to-end check with a physically shaped Rayleigh-Plateau signal.

    Plant the *exact* inviscid growth rate for k r0 = 0.7 into a stationary sinusoidal interface perturbation, add faster-decaying contamination, and confirm DMD returns the dispersion-relation value.
    """
    r0, k = 1.0, 0.7
    sigma_exact = float(dispersion.growth_rate(k, r0=r0))

    nz, n_snap, dt = 256, 300, 0.05
    z = np.linspace(0.0, 2.0 * np.pi / k, nz, endpoint=False)

    interface = np.cos(k * z)[:, None]                 # the capillary mode
    decaying = np.stack([np.cos(3 * k * z), np.sin(5 * k * z)], axis=1)
    shapes = np.hstack([interface, decaying])

    X = _snapshots_from_modes(
        shapes,
        [sigma_exact, -2.0, -3.5],
        [0.0, 0.0, 0.7],
        [1e-5, 5e-6, 5e-6],
        dt, n_snap,
    ).real

    res = dmd_mod.dmd(X, dt)
    measured = res.leading_growth_rate()
    assert measured == pytest.approx(sigma_exact, rel=1e-6)
    # and it should be the physically expected O(0.3) value
    assert 0.2 < measured < 0.4
