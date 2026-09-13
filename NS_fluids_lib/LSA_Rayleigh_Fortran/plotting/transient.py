"""Has the growth rate settled, or is the run still in its transient?

This module exists because a measurement in this study was reported without it and was wrong by about four per cent.

A perturbed calculation does not grow at the modal rate from the first instant. The imposed shape is a single Fourier mode of the interface, but the velocity field is not the corresponding eigenvector, so the initial state is a combination of the growing mode and a set of decaying ones. Only after the decaying part has died away does the amplitude grow at the eigenvalue. Fitting an exponential before that point measures a mixture, and because the decaying part reduces the apparent rate, the fit is biased low.

The bias is not obvious in the amplitude curve, which looks like a clean exponential on a logarithmic axis long before it is one. It is obvious in the local growth rate, which rises towards a plateau; if the curve is still climbing at the end of the record, the run was stopped too early and no rate from it should be quoted. That is the check performed here, and it is the check that was missing.

The transient is modelled as a single decaying exponential,

    sigma(t) = sigma_inf - A exp(-beta t),

which is the leading behaviour once the slowest decaying mode dominates the remainder. Three quantities follow. The asymptote sigma_inf estimates what the rate would settle to; the residual A exp(-beta t_end) states how much transient is still present at the end of the record, which is the bias in the reported rate; and the time at which that residual falls below a chosen tolerance gives the run length actually required.

The model is deliberately simple and should not be trusted far beyond the data. Extrapolating an asymptote from a curve that is still steeply rising is an extrapolation, and its purpose here is diagnostic rather than quantitative: it answers whether the run is long enough, and if not roughly how much longer it must be. The rate that goes into the thesis should come from a run whose local growth rate has visibly flattened, not from this extrapolation.
"""

# ORIENTATION
#
# If you read one module in this project to understand why measurements go wrong, read this one. It exists because a rate was reported from a run that had not finished settling, and was wrong by four per cent.


from __future__ import annotations

import numpy as np

#: Residual transient, as a fraction of the asymptote, that the fit window is required to open below. One per cent is a deliberate compromise rather than an ideal. Demanding less pushes the opening of the window towards the end of the record and leaves too few samples to fit, and the resulting narrow window trades a transient bias for a sampling one; demanding more admits contamination that the fit averages straight into the answer. At one per cent the bias on the fitted rate is a few tenths of a per cent, which is below the discretisation error at every resolution used here, so the measurement is limited by the mesh rather than by the window.
DEFAULT_TOLERANCE = 1.0e-2


#: Amplitude, as a fraction of the column radius, beyond which the linear description is not claimed to hold. The proportionality test in this study holds to about a tenth of the radius; past that the second-order term is no longer negligible, harmonics appear, and the mean radius drifts.
LINEAR_BOUND = 0.1


def linear_end(times, amplitudes, r0=1.0, bound=LINEAR_BOUND):
    """Last instant at which the amplitude is still small enough to be linear.

    This is the upper limit of any window a growth rate may be fitted over, and omitting it is how a correct diagnosis of one failure produces another. A run continued past this point does not merely stop being informative: the amplitude saturates and then oscillates as the column breaks and the fragments relax, so the local growth rate first exceeds the eigenvalue and then swings wildly. A window chosen as a fixed fraction of the record then straddles the break-up, and the fit returns a number that belongs to neither regime.

    Returns the time itself and whether the bound was reached at all. A record that never reaches it is entirely linear, which is the ordinary case for a short run and for every calculation in the stable band.
    """
    times = np.asarray(times, dtype=float)
    amplitudes = np.asarray(amplitudes, dtype=float)
    limit = bound * r0

    beyond = np.where(np.isfinite(amplitudes) & (amplitudes > limit))[0]
    if beyond.size == 0:
        return float(times[-1]), False
    first = int(beyond[0])
    if first == 0:
        # The very first sample already exceeds the bound, so no part of the record is linear. The caller must be told rather than handed the initial time as though it were a window.
        return float(times[0]), True
    return float(times[first - 1]), True


def select_window(times, amplitudes, r0=1.0, bound=LINEAR_BOUND,
                  keep=0.45, width=0.8, tolerance=DEFAULT_TOLERANCE):
    """Fit window bounded on both sides: after the transient, before nonlinearity.

    The upper bound comes from `linear_end`. The lower bound is set from the measured decay of the transient wherever that can be fitted, and only otherwise from a fixed fraction of the linear record.

    That distinction is the difference between a rate biased by a per cent and one biased by three. A least-squares fit over a window returns something close to the average growth rate across it, so a window that opens while the transient still contributes several per cent carries that contamination into the answer no matter how well the run has settled by the time it ends. Choosing the opening time from the fitted decay, as the instant at which the residual transient falls below the tolerance, makes the bias on the fitted rate comparable to the tolerance rather than to the transient at the start of the window.

    Where the two bounds leave too few samples, the window is widened back towards the fraction rule rather than returned empty, and the caller is told by `wide_enough` that the result carries more transient than requested.
    """
    times = np.asarray(times, dtype=float)
    amplitudes = np.asarray(amplitudes, dtype=float)
    t_hi, reached = linear_end(times, amplitudes, r0=r0, bound=bound)
    t_start = float(times[0])

    # The fraction rule is the fallback, and also the floor: the transient rule must not open the window earlier than this, only later.
    t_lo = t_start + keep * (t_hi - t_start)

    linear = times <= t_hi
    fit = fit_amplitude_model(times[linear], amplitudes[linear])
    if fit is not None and fit["beta"] > 0.0 and fit["sigma_inf"] > 0.0:
        # Residual transient relative to the asymptote is (A/sigma_inf) exp(-beta t); invert for the instant at which it falls below the tolerance.
        ratio = fit["amp"] / fit["sigma_inf"]
        if ratio > tolerance:
            t_settled = float(np.log(ratio / tolerance) / fit["beta"])
            t_lo = max(t_lo, min(t_settled, t_hi))

    wide_enough = (t_hi - t_lo) >= 2.0 * width
    if not wide_enough:
        # Too narrow to fit over, so the window is widened back; the flag records that the answer carries more transient than the tolerance allows.
        t_lo = min(t_lo, max(t_start, t_hi - 2.0 * width))
    return t_lo, t_hi, reached, wide_enough


def local_growth_rate(times, amplitudes, width, n_windows=8):
    """Growth rate fitted over sliding windows of the record.

    Each window gives one estimate, obtained by fitting a straight line to the logarithm of the amplitude. Plotted against the window centre these estimates form the curve whose flattening signals that the transient has gone.
    """
    times = np.asarray(times, dtype=float)
    amplitudes = np.asarray(amplitudes, dtype=float)

    lo = times.min() + 0.5 * width
    hi = times.max() - 0.5 * width
    if hi <= lo:
        return np.array([]), np.array([])

    centres, rates = [], []
    for centre in np.linspace(lo, hi, n_windows):
        inside = (times >= centre - 0.5 * width) & (times <= centre + 0.5 * width)
        inside &= np.isfinite(amplitudes) & (amplitudes > 0.0)
        if inside.sum() >= 4:
            # A straight line through log a has slope sigma, since a = a0 exp(sigma t) implies log a = log a0 + sigma t.
            slope = np.polyfit(times[inside], np.log(amplitudes[inside]), 1)[0]
            centres.append(centre)
            rates.append(slope)
    return np.asarray(centres), np.asarray(rates)


def fit_approach(centres, rates):
    """Fit sigma(t) = sigma_inf - A exp(-beta t) to the local rate.

    Returns the three parameters, or None if the fit does not converge, which happens when the record is too short or too noisy to constrain three unknowns. Returning None rather than an unconverged answer matters, because the caller uses this to decide whether a measurement may be quoted.
    """
    from scipy.optimize import curve_fit

    centres = np.asarray(centres, dtype=float)
    rates = np.asarray(rates, dtype=float)
    if centres.size < 4:
        return None

    def model(t, sigma_inf, amp, beta):
        return sigma_inf - amp * np.exp(-beta * t)

    # The starting guess matters: the asymptote is at least the largest observed rate, the transient amplitude is roughly the range covered, and the decay rate is of order the reciprocal of the record length.
    guess = [float(rates.max()) * 1.05,
             float(rates.max() - rates.min()),
             # ndarray.ptp was removed in NumPy 2; the free function is the portable spelling.
             1.0 / max(float(np.ptp(centres)), 1.0)]
    try:
        params, _ = curve_fit(model, centres, rates, p0=guess, maxfev=20000)
    except (RuntimeError, ValueError):
        return None

    sigma_inf, amp, beta = (float(v) for v in params)
    if beta <= 0.0 or not np.isfinite(sigma_inf):
        # A non-positive decay rate means the model has fitted a growing transient, which is not the situation this diagnostic describes.
        return None
    return sigma_inf, amp, beta



def fit_amplitude_model(times, amplitudes, sigma_max=None):
    """Fit the transient model to the amplitude record itself.

    Integrating d(log a)/dt = sigma_inf - A exp(-beta t) gives

        log a(t) = log a_0 + sigma_inf t - (A/beta)(1 - exp(-beta t)),

    which is fitted here to every sample in the linear window. That is a better estimator than fitting the same model to a handful of sliding-window growth rates, and the difference is not cosmetic. Each sliding-window rate is itself a fit, so using eight of them discards all but eight numbers from a record of several dozen, and the differencing that produces them amplifies noise. Fitted to the amplitudes directly, the same three physical parameters are constrained by every sample, and the extrapolation is correspondingly better conditioned.

    This matters most in exactly the case that motivated it: a window short enough that the local growth rate has not flattened. Fitted to the rates, the model then extrapolated to values above the inviscid limit, which is impossible. Fitted to the amplitudes it does not.

    `sigma_max` is the largest physically admissible rate, normally the inviscid one, since viscosity can only retard a capillary instability. It is not imposed on the fit but is used to judge the result: an estimate above it is returned with `physical` false, because a bound violated is evidence that the extrapolation has failed rather than something to be silently clipped.
    """
    from scipy.optimize import curve_fit

    times = np.asarray(times, dtype=float)
    amplitudes = np.asarray(amplitudes, dtype=float)
    good = np.isfinite(amplitudes) & (amplitudes > 0.0)
    if good.sum() < 8:
        return None
    t, y = times[good], np.log(amplitudes[good])

    def model(tt, log_a0, sigma_inf, amp, beta):
        return log_a0 + sigma_inf * tt - (amp / beta) * (1.0 - np.exp(-beta * tt))

    # A straight line through the later half gives the scale of the rate, and the shortfall at early times gives the scale of the transient.
    half = t > 0.5 * (t[0] + t[-1])
    slope_guess = np.polyfit(t[half], y[half], 1)[0] if half.sum() > 2 else 0.3
    guess = [float(y[0]), float(slope_guess), float(max(slope_guess, 0.1)),
             1.0 / max(float(np.ptp(t)), 1.0)]
    try:
        params, _ = curve_fit(model, t, y, p0=guess, maxfev=40000,
                              bounds=([-np.inf, 0.0, 0.0, 1e-3],
                                      [np.inf, np.inf, np.inf, 50.0]))
    except (RuntimeError, ValueError):
        return None

    log_a0, sigma_inf, amp, beta = (float(v) for v in params)
    residual = y - model(t, *params)
    ss_res = float(np.sum(residual ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    physical = True if sigma_max is None else sigma_inf <= sigma_max * 1.001
    return {"sigma_inf": sigma_inf, "amp": amp, "beta": beta,
            "log_a0": log_a0, "r2": r2, "physical": physical}


def residual_transient(amp, beta, t_end):
    """Transient still present at the end of the record, in absolute units."""
    return float(amp * np.exp(-beta * t_end))


def required_time(amp, beta, sigma_inf, tolerance=DEFAULT_TOLERANCE):
    """Time at which the remaining transient falls below a tolerance.

    Obtained by inverting A exp(-beta t) = tolerance . sigma_inf.
    """
    target = tolerance * abs(sigma_inf)
    if target <= 0.0 or amp <= 0.0:
        return float("nan")
    return float(np.log(amp / target) / beta)


def assess(times, amplitudes, width=0.8, tolerance=DEFAULT_TOLERANCE,
           sigma_reference=None, r0=1.0, bound=LINEAR_BOUND,
           sigma_max=None):
    """Full diagnosis of whether a record is long enough to quote a rate.

    The returned dictionary carries `settled`, which is the answer to the question the caller is really asking. Everything else supports it: the local rate curve for plotting, the fitted asymptote, the bias remaining in the record, and the run length that would remove it.
    """
    times = np.asarray(times, dtype=float)
    amplitudes = np.asarray(amplitudes, dtype=float)

    # The record is cut back to its linear part before anything is fitted. Assessing the settling of a transient on a record that includes break-up fits the model to the wrong thing entirely, and the earlier version of this function did exactly that: it returned an asymptote above the inviscid limit, which is impossible, because the nonlinear rise dominated the fit.
    t_hi, reached = linear_end(times, amplitudes, r0=r0, bound=bound)
    linear = times <= t_hi
    times_lin = times[linear]
    amps_lin = amplitudes[linear]

    centres, rates = local_growth_rate(times_lin, amps_lin, width)

    result = {
        "centres": centres,
        "rates": rates,
        "settled": False,
        "sigma_inf": float("nan"),
        "beta": float("nan"),
        "residual": float("nan"),
        "relative_bias": float("nan"),
        "t_required": float("nan"),
        "nonlinear_reached": reached,
        "t_linear_end": t_hi,
        "message": "",
    }
    if reached and t_hi <= float(times[0]):
        result["message"] = (
            "THE PERTURBATION IS TOO LARGE. The amplitude exceeds "
            f"{bound:g} r0 from the first sample, so no part of this record "
            "is linear and no growth rate can be fitted. Reduce eps.")
        return result
    if centres.size < 4:
        if reached:
            result["message"] = (
                "THE PERTURBATION IS TOO LARGE. The amplitude leaves the "
                f"linear range at t = {t_hi:.2f}, which is too early for the "
                "transient to have decayed, so there is no window in which a "
                "rate may be fitted. Reduce eps rather than lengthening the "
                "run.")
        else:
            result["message"] = ("record too short to assess; fewer than four "
                                 "sliding windows fit inside it")
        return result

    # The final rise, expressed against the spread of the whole curve, is a model-free indicator. A curve that has flattened shows a last increment far smaller than its total rise.
    total_rise = float(rates.max() - rates.min())
    last_step = float(rates[-1] - rates[-2])
    result["last_increment"] = last_step
    result["total_rise"] = total_rise

    # Fitted to the amplitudes, which uses every sample; the rate-curve fit is kept only as a fallback for a record too short for the four-parameter form.
    direct = fit_amplitude_model(times_lin, amps_lin, sigma_max=sigma_max)
    if direct is not None:
        fit = (direct["sigma_inf"], direct["amp"], direct["beta"])
        result["model_r2"] = direct["r2"]
        result["physical"] = direct["physical"]
    else:
        fit = fit_approach(centres, rates)
        result["physical"] = True
    if fit is None:
        result["message"] = ("the approach to a plateau could not be fitted; "
                             "inspect the local growth rate by eye")
        return result

    sigma_inf, amp, beta = fit
    # The end of the linear part, not the end of the record, since the rate is only ever quoted from inside the linear part.
    t_end = float(times_lin.max())
    residual = residual_transient(amp, beta, t_end)

    result.update({
        "sigma_inf": sigma_inf,
        "beta": beta,
        "residual": residual,
        "relative_bias": residual / abs(sigma_inf) if sigma_inf else np.nan,
        "t_required": required_time(amp, beta, sigma_inf, tolerance),
        "settled": residual / abs(sigma_inf) < tolerance if sigma_inf else False,
    })

    if not result.get("physical", True):
        # An estimate above the inviscid bound is not a result to be clipped; it means the window is too short to support the extrapolation at all.
        result["settled"] = False
        result["message"] = (
            f"THE CORRECTION IS NOT TRUSTWORTHY. The extrapolated rate, "
            f"{result['sigma_inf']:.4f}, exceeds the inviscid bound of "
            f"{sigma_max:.4f}, which no viscous flow can. The linear window "
            f"is too short to constrain the extrapolation. Reduce eps so that "
            f"the amplitude stays below {bound:g} r0 until about t = "
            f"{result['t_required']:.1f}, and quote the direct fit rather "
            f"than this correction in the meantime.")
        return result

    if result["settled"]:
        note = ""
        if reached:
            note = (f" The record continues past t = {t_hi:.2f}, where the "
                    "amplitude leaves the linear range; that part is excluded.")
        result["message"] = (
            f"the local growth rate has settled: the remaining transient is "
            f"{100 * result['relative_bias']:.2f} per cent of the asymptote, "
            f"below the {100 * tolerance:.2f} per cent tolerance." + note)
    elif reached:
        # Both bounds bite at once, which is the situation a longer run cannot fix.
        result["message"] = (
            f"NO USABLE WINDOW. The transient still contributes about "
            f"{100 * result['relative_bias']:.2f} per cent at t = {t_hi:.2f}, "
            f"where the amplitude already leaves the linear range. Lengthening "
            f"the run will not help, because the amplitude grows past "
            f"{bound:g} r0 before the transient decays. Reduce eps: the "
            f"transient needs until t = {result['t_required']:.2f}, so eps "
            f"must be small enough that the amplitude is still below "
            f"{bound:g} r0 then.")
    else:
        result["message"] = (
            f"THE RUN IS TOO SHORT. The local growth rate is still rising, and "
            f"about {100 * result['relative_bias']:.2f} per cent of transient "
            f"remains at t = {t_end:.2f}, which biases the fitted rate low by "
            f"roughly that amount. Advance the calculation to t = "
            f"{result['t_required']:.2f} before quoting a rate.")

    if sigma_reference is not None and np.isfinite(sigma_inf):
        result["reference_gap"] = (sigma_inf - sigma_reference) / sigma_reference

    return result


def n_periods_required(t_required, sigma_reference):
    """Run length in e-folding times, the unit the deck generator expects.

    `make_inputs` measures the stop time in e-folding times of the reference rate, so a required time is converted here rather than in the caller, where the convention would have to be restated.
    """
    if not np.isfinite(t_required) or sigma_reference <= 0.0:
        return float("nan")
    return float(t_required * sigma_reference)


def format_report(result, sigma_reference=None):
    """The diagnosis as plain text, for printing beneath a measurement."""
    lines = ["-" * 66, "transient convergence", "-" * 66]
    if np.isfinite(result.get("sigma_inf", np.nan)):
        lines.append(f"{'extrapolated asymptote':<28} "
                     f"{result['sigma_inf']:.6f}")
        lines.append(f"{'transient decay rate beta':<28} "
                     f"{result['beta']:.4f}")
        lines.append(f"{'transient left in the record':<28} "
                     f"{100 * result['relative_bias']:.2f} per cent")
        lines.append(f"{'time needed for convergence':<28} "
                     f"{result['t_required']:.2f}")
        if sigma_reference:
            lines.append(f"{'as e-folding times':<28} "
                         f"{n_periods_required(result['t_required'], sigma_reference):.2f}")
    lines.append("")
    lines.append(result["message"])
    lines.append("-" * 66)
    return "\n".join(lines)
