"""Tests for the transient-convergence diagnostic.

The diagnostic answers one question: may a growth rate from this record be quoted? A wrong answer in either direction is costly. Passing a short record lets a biased rate into a thesis, which is the failure this module was written after. Failing a converged one wastes cluster time on a run that was already long enough.

Both directions are therefore tested against records built from a known transient, where the correct verdict is not a matter of judgement. A record that carries a large decaying contamination must be rejected; one carried far enough for that contamination to vanish must be accepted; and the run length the diagnostic prescribes must be long enough that a record of that length passes.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "python"))

import transient  # noqa: E402

SIGMA = 0.3294          # the asymptotic rate planted in every synthetic record
BETA = 0.45             # decay rate of the contaminating transient


def synthetic(t_end, sigma=SIGMA, beta=BETA, amp=0.35, n=160, eps=2e-2):
    """A record whose local growth rate approaches `sigma` from below.

    Integrating d(log a)/dt = sigma - amp exp(-beta t) gives the amplitude in closed form, so the record has exactly the local growth rate the diagnostic is meant to recover, and the test is not circular.
    """
    t = np.linspace(0.0, t_end, n)
    log_a = np.log(eps) + sigma * t + (amp / beta) * (np.exp(-beta * t) - 1.0)
    return t, np.exp(log_a)


def test_rejects_a_record_stopped_in_the_transient():
    t, a = synthetic(6.7)
    result = transient.assess(t, a)
    assert not result["settled"]
    assert "TOO SHORT" in result["message"]
    assert result["relative_bias"] > 0.01


def test_accepts_a_record_carried_to_the_plateau():
    # The amplitude is reduced along with the lengthening, since at the reference amplitude a record this long would leave the linear range and be truncated. That coupling between run length and amplitude is the point of the two-sided window.
    t, a = synthetic(30.0, eps=2e-6)
    result = transient.assess(t, a)
    assert result["settled"]
    assert result["relative_bias"] < transient.DEFAULT_TOLERANCE


def test_recovers_the_planted_asymptote_and_decay_rate():
    t, a = synthetic(20.0)
    result = transient.assess(t, a)
    assert result["sigma_inf"] == pytest.approx(SIGMA, rel=0.02)
    assert result["beta"] == pytest.approx(BETA, rel=0.15)


def test_the_prescribed_run_length_is_long_enough():
    """The remedy must actually work, which is the point of prescribing it."""
    t, a = synthetic(6.7, eps=2e-6)
    short = transient.assess(t, a)
    assert not short["settled"]

    t_long, a_long = synthetic(short["t_required"] * 1.05, eps=2e-6)
    assert transient.assess(t_long, a_long)["settled"]


def test_a_short_record_biases_the_fitted_rate_low():
    """The direction of the bias matters, not only its size.

    A decaying transient subtracts from the local rate, so an exponential fitted inside it must underestimate. Were the bias in the other direction, agreement with theory would have been evidence of nothing.
    """
    t, a = synthetic(6.7)
    late = t >= 0.65 * t.max()
    fitted = np.polyfit(t[late], np.log(a[late]), 1)[0]
    assert fitted < SIGMA
    assert (SIGMA - fitted) / SIGMA > 0.01


def test_n_periods_conversion():
    """The prescription must be expressed in the deck generator's own unit."""
    assert transient.n_periods_required(15.28, 0.3294) == pytest.approx(5.03,
                                                                        abs=0.02)
    assert np.isnan(transient.n_periods_required(float("nan"), 0.3294))


def test_short_record_is_reported_rather_than_guessed():
    t = np.linspace(0.0, 0.2, 6)
    # Scaled below the linear bound, so the record is refused for being short rather than for being nonlinear.
    result = transient.assess(t, 1e-3 * np.exp(0.3 * t), width=0.8)
    assert not result["settled"]
    assert "too short to assess" in result["message"]


def test_local_growth_rate_recovers_a_pure_exponential():
    t = np.linspace(0.0, 10.0, 200)
    centres, rates = transient.local_growth_rate(t, np.exp(0.3 * t), width=1.0)
    assert centres.size > 0
    assert np.allclose(rates, 0.3, atol=1e-10)


# ----------------------------------------------------------------------
# The other failure: a run carried past the linear range.
#
# The first version of this module tested only whether a run was too short. Applied to a record carried three times further than linearity permits, it fitted its model to the nonlinear blow-up and reported an asymptote above the inviscid limit, which is physically impossible. These tests pin the second bound.
# ----------------------------------------------------------------------
def with_breakup(t_end=15.13, sigma=SIGMA, beta=BETA, amp=0.33, eps=2e-2,
                 n=125):
    """A record that grows, plateaus, then saturates and oscillates.

    Past an amplitude of order the radius the column breaks and the fragments relax, so the amplitude stops growing and swings instead. That is what the measured records do, and it is what defeats a one-sided window.
    """
    t = np.linspace(0.0, t_end, n)
    log_a = np.log(eps) + sigma * t + (amp / beta) * (np.exp(-beta * t) - 1.0)
    a = np.exp(log_a)
    late = a > 0.8
    if late.any():
        a[late] = 0.8 + 0.7 * np.abs(np.sin(1.9 * (t[late] - t[late][0])))
    return t, a


def test_linear_end_finds_the_bound():
    t, a = with_breakup()
    t_hi, reached = transient.linear_end(t, a, r0=1.0, bound=0.1)
    assert reached
    index = int(np.argmin(np.abs(t - t_hi)))
    assert a[index] <= 0.1
    assert a[index + 1] > 0.1


def test_linear_end_reports_a_wholly_linear_record():
    t, a = synthetic(6.0, eps=1e-4)
    t_hi, reached = transient.linear_end(t, a, r0=1.0, bound=0.1)
    assert not reached
    assert t_hi == pytest.approx(t[-1])


def test_window_excludes_the_nonlinear_part():
    t, a = with_breakup()
    t_lo, t_hi, reached, wide = transient.select_window(t, a)
    assert reached
    assert t_hi < t[-1]
    inside = (t >= t_lo) & (t <= t_hi)
    assert (a[inside] <= 0.1).all()


def test_window_beats_the_one_sided_rule_on_a_broken_record():
    """The point of the change, stated as a test.

    A window set at a fixed fraction of the record straddles break-up and fits a straight line through data that is not exponential at all, which shows up as a collapsed coefficient of determination. The two-sided window stays inside the linear part and fits well.
    """
    t, a = with_breakup()

    t_lo, t_hi, _, _ = transient.select_window(t, a)
    inside = (t >= t_lo) & (t <= t_hi)
    y = np.log(a[inside])
    fitted = np.polyval(np.polyfit(t[inside], y, 1), t[inside])
    r2_two_sided = 1.0 - np.sum((y - fitted) ** 2) / np.sum((y - y.mean()) ** 2)

    old = t >= 0.65 * t.max()
    y_old = np.log(a[old])
    fitted_old = np.polyval(np.polyfit(t[old], y_old, 1), t[old])
    r2_one_sided = 1.0 - (np.sum((y_old - fitted_old) ** 2)
                          / np.sum((y_old - y_old.mean()) ** 2))

    assert r2_two_sided > 0.99
    assert r2_one_sided < r2_two_sided


def test_asymptote_stays_physical_on_a_broken_record():
    """The specific bug: an extrapolated rate above the inviscid limit.

    Fitted to the whole record, including break-up, the model returned 0.437 against an inviscid bound of 0.343. Truncating to the linear part first must keep the answer below that bound.
    """
    t, a = with_breakup()
    result = transient.assess(t, a)
    assert np.isfinite(result["sigma_inf"])
    assert result["sigma_inf"] < 0.3434


def test_asymptote_recovers_the_planted_rate_from_a_broken_record():
    """And is not merely bounded, but close.

    This is what makes the correction usable: where no clean window exists, the extrapolated asymptote is still within a fraction of a per cent of the rate actually planted.
    """
    t, a = with_breakup()
    result = transient.assess(t, a)
    assert result["sigma_inf"] == pytest.approx(SIGMA, rel=0.01)


def test_no_usable_window_is_reported_as_such():
    """Both bounds biting at once needs a different remedy from either alone."""
    t, a = with_breakup()
    result = transient.assess(t, a)
    assert not result["settled"]
    assert "NO USABLE WINDOW" in result["message"]
    # The remedy must be to reduce the amplitude, not to lengthen the run.
    assert "Reduce eps" in result["message"]
    assert "Advance the calculation" not in result["message"]


def test_an_amplitude_too_large_from_the_start_is_refused():
    t = np.linspace(0.0, 5.0, 60)
    a = 0.5 * np.exp(0.3 * t)
    result = transient.assess(t, a)
    assert not result["settled"]
    assert "TOO LARGE" in result["message"]
    assert not np.isfinite(result["sigma_inf"])


def test_a_small_amplitude_run_still_gets_a_clean_window():
    """Reducing eps is the prescribed remedy, so it must work.

    At a hundredth of the reference amplitude the transient decays long before nonlinearity, so the record settles and the direct fit recovers the planted rate without any correction.
    """
    t, a = with_breakup(t_end=30.0, eps=2e-4, n=300)
    result = transient.assess(t, a)
    assert result["settled"]
    t_lo, t_hi, _, wide = transient.select_window(t, a)
    assert wide
    inside = (t >= t_lo) & (t <= t_hi)
    slope = np.polyfit(t[inside], np.log(a[inside]), 1)[0]
    assert slope == pytest.approx(SIGMA, rel=0.02)
