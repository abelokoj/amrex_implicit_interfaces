"""Tests for the proportionality (linearity) assessment.

Synthetic amplitude sweeps are constructed with a *known* regime - clean linear, nonlinearly contaminated, or round-off dominated - and the assessment must both flag the failure and attribute it to the right cause, since the corrective action (reduce vs increase eps) is opposite in the two cases.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

import proportionality as prop  # noqa: E402
import nsmfp  # noqa: E402


def make_run(eps, times, rate=0.3, quad=0.0, noise=0.0, seed=0):
    """|Q'|(t) = eps e^{rate t} + quad (eps e^{rate t})^2 + round-off noise."""
    rng = np.random.default_rng(seed)
    linear = eps * np.exp(rate * times)
    norms = linear + quad * linear**2
    if noise:
        norms = norms + noise * rng.standard_normal(times.size)
    norms = np.abs(norms)
    r, _, r2 = nsmfp.fit_exponential_growth(times, norms)
    return prop.AmplitudeRun(eps=eps, times=times, norms=norms,
                             growth_rate=r, r_squared=r2)


def test_clean_linear_sweep_is_linear():
    t = np.linspace(0, 5, 120)
    runs = [make_run(e, t) for e in (1e-6, 1e-5, 1e-4)]
    rep = prop.assess_linearity(runs)
    assert rep.is_linear
    assert rep.verdict == "linear"
    assert rep.scaled_spread < 1e-9
    assert rep.growth_rate_mean == pytest.approx(0.3, rel=1e-6)


def test_growth_rate_recovered_independent_of_amplitude():
    t = np.linspace(0, 4, 100)
    runs = [make_run(e, t, rate=0.42) for e in (1e-7, 1e-6, 1e-5)]
    for r in runs:
        assert r.growth_rate == pytest.approx(0.42, rel=1e-6)
        assert r.r_squared == pytest.approx(1.0, abs=1e-9)


def test_nonlinear_contamination_is_detected_and_attributed():
    """Large amplitudes with a quadratic term: must say 'reduce eps'."""
    t = np.linspace(0, 5, 120)
    # quadratic term matters only for the largest amplitude, and only late
    runs = [make_run(1e-6, t, quad=50.0),
            make_run(1e-3, t, quad=50.0),
            make_run(1e-2, t, quad=50.0)]
    rep = prop.assess_linearity(runs)
    assert not rep.is_linear
    assert rep.verdict == "nonlinear"
    assert "REDUCE" in rep.detail
    assert prop.recommend_amplitude(rep) < min(r.eps for r in runs)


def test_round_off_is_detected_and_attributed():
    """Tiny amplitudes with additive noise: must say 'increase eps'."""
    t = np.linspace(0, 5, 120)
    runs = [make_run(1e-14, t, noise=1e-14, seed=1),
            make_run(1e-13, t, noise=1e-14, seed=2),
            make_run(1e-12, t, noise=1e-14, seed=3)]
    rep = prop.assess_linearity(runs)
    assert not rep.is_linear
    assert rep.verdict == "round-off"
    assert "INCREASE" in rep.detail
    assert prop.recommend_amplitude(rep) > max(r.eps for r in runs)


def test_recommend_amplitude_picks_largest_linear():
    t = np.linspace(0, 5, 60)
    runs = [make_run(e, t) for e in (1e-7, 1e-6, 1e-5)]
    rep = prop.assess_linearity(runs)
    assert prop.recommend_amplitude(rep) == 1e-5


def test_requires_at_least_two_amplitudes():
    t = np.linspace(0, 5, 60)
    with pytest.raises(ValueError):
        prop.assess_linearity([make_run(1e-5, t)])


def test_rejects_nonpositive_amplitude():
    with pytest.raises(ValueError):
        prop.collect_amplitude_runs({0.0: "somewhere"}, "base", verbose=False)


def test_report_renders():
    t = np.linspace(0, 5, 60)
    rep = prop.assess_linearity([make_run(e, t) for e in (1e-6, 1e-5)])
    s = str(rep)
    assert "Proportionality" in s and "LINEAR" in s


def test_fit_exponential_growth_windowing():
    """Restricting the window should exclude a contaminated tail."""
    t = np.linspace(0, 10, 400)
    clean = 1e-6 * np.exp(0.3 * t)
    contaminated = clean + np.where(t > 6, 1e-3 * (t - 6) ** 3, 0.0)
    r_all, _, _ = nsmfp.fit_exponential_growth(t, contaminated)
    r_win, _, r2 = nsmfp.fit_exponential_growth(t, contaminated, t_max=5.0)
    assert r_win == pytest.approx(0.3, rel=1e-6)
    assert abs(r_all - 0.3) > abs(r_win - 0.3)
    assert r2 > 0.999


def test_uniform_dt_check_rejects_nonuniform_sampling():
    with pytest.raises(ValueError, match="not uniform"):
        nsmfp._uniform_dt(np.array([0.0, 0.1, 0.2, 0.35]))
    assert nsmfp._uniform_dt(np.array([0.0, 0.1, 0.2, 0.3])) == pytest.approx(0.1)


def test_uniform_dt_rejects_non_increasing():
    with pytest.raises(ValueError, match="increasing"):
        nsmfp._uniform_dt(np.array([0.0, 0.1, 0.05]))
