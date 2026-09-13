"""
Proportionality (linearity) test for the NS-MFP subspace.

Why this is load-bearing rather than optional
---------------------------------------------
NS-MFP marches the *nonlinear* solver.  Its output is a Jacobian-vector product only insofar as the O(eps^2) term in

    dQ'/dt = Q' dF/dQ + Q'^2 d2F/dQ2 + H.O.T.            (Ranjan Eq. 7)

is negligible.  Ranjan et al. state this plainly: "linearity is ensured through simple tests based on proportionality between input and output ... Although the value of eps_0 should be evaluated for each problem".  Nothing in the algorithm detects a violation; if eps_0 is too large the DMD still returns a clean-looking spectrum, but it is the spectrum of a contaminated operator.

This test therefore is not a verification afterthought but the evidence that the operator being decomposed is the linearised one.  It also distinguishes this approach from methods that infer snapshots from nonlinear simulation data, where no such proportionality is available.

The test
--------
Run the same perturbation shape at several amplitudes eps_i.  Under linearity,

    Q'(t; eps_i) = (eps_i / eps_j) Q'(t; eps_j)     for all i, j, t,

so (a) the norm |Q'(t; eps)| must be proportional to eps at every time, and (b) the normalised fields Q'/eps must collapse onto a single curve, and (c) the extracted growth rate must be independent of eps.

Failures have distinct signatures:
  * eps too LARGE  -> normalised curves fan out at late time as the nonlinear
                      term bites; growth rate drifts with eps.
  * eps too SMALL  -> curves are noisy from the start and the collapse
                      degrades at *early* time; round-off, not physics.
Both are reported separately by `assess_linearity`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

import report

import nsmfp


@dataclass
class AmplitudeRun:
    """One member of the amplitude sweep."""

    eps: float
    times: np.ndarray
    norms: np.ndarray
    growth_rate: float = float("nan")
    r_squared: float = float("nan")
    snapshots: np.ndarray = field(default=None, repr=False)


@dataclass
class LinearityReport:
    runs: list
    scaled_spread: float          # max relative spread of |Q'|/eps
    growth_rate_spread: float     # max-min of fitted growth rates
    growth_rate_mean: float
    verdict: str
    detail: str

    @property
    def is_linear(self):
        return self.verdict == "linear"

    def __str__(self):
        rows = []
        for r in self.runs:
            final = r.norms[-1] / r.eps if r.eps else float("nan")
            rows.append([f"{r.eps:.3e}", f"{r.growth_rate:.6f}",
                         f"{r.r_squared:.5f}", f"{final:.3e}"])
        body = report.table(
            rows, ["eps", "growth rate", "r^2", "|Q'|/eps final"],
            colalign=["right"] * 4, title="Proportionality (linearity) test")
        summary = report.kv([
            ("scaled-norm spread", report.pct(self.scaled_spread, 4)),
            ("growth-rate spread", f"{self.growth_rate_spread:.3e}"),
            ("mean growth rate", f"{self.growth_rate_mean:.3e}"),
            ("verdict", self.verdict.upper()),
        ])
        return "\n".join([body, "", summary, "", self.detail])


def collect_amplitude_runs(cases, unperturbed_dir, fields=nsmfp.DEFAULT_FIELDS,
                           level=0, t_min=None, t_max=None, verbose=True):
    """Read an amplitude sweep.

    Parameters
    ----------
    cases : mapping {eps: perturbed_run_directory}
        Each directory holds the plt* output of a run that differs from the
        others *only* in perturbation amplitude.
    unperturbed_dir : str
        The shared twin (base-state) run.
    """
    runs = []
    for eps in sorted(cases):
        if eps <= 0:
            raise ValueError(f"perturbation amplitude must be positive, got {eps}")
        d = cases[eps]
        if verbose:
            print(f"[eps={eps:g}] reading {d}")
        X, dt, times = nsmfp.build_snapshot_matrix(
            d, unperturbed_dir, fields=fields, level=level, verbose=False)
        norms = np.linalg.norm(X, axis=0)
        try:
            rate, _, r2 = nsmfp.fit_exponential_growth(
                times, norms, t_min=t_min, t_max=t_max)
        except ValueError:
            rate, r2 = float("nan"), float("nan")
        runs.append(AmplitudeRun(eps=eps, times=times, norms=norms,
                                 growth_rate=rate, r_squared=r2, snapshots=X))
    return runs


def assess_linearity(runs, scaled_tol=0.02, rate_tol=0.02,
                     early_frac=0.25):
    """Judge whether an amplitude sweep is in the linear regime.

    Parameters
    ----------
    runs : list of AmplitudeRun scaled_tol : float
        Allowed relative spread of |Q'|/eps across amplitudes (default 2%).
    rate_tol : float
        Allowed spread of fitted growth rates, relative to their mean.
    early_frac : float
        Fraction of the record treated as "early" when attributing a failure
        to round-off (early degradation) rather than nonlinearity (late
        degradation).
    """
    if len(runs) < 2:
        raise ValueError("need at least two amplitudes to test proportionality")

    n = min(len(r.times) for r in runs)
    if n < 2:
        raise ValueError("runs have too few snapshots in common")

    # |Q'|/eps for each run on the common time window
    scaled = np.vstack([r.norms[:n] / r.eps for r in runs])
    ref = scaled.mean(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.abs(scaled - ref) / np.where(ref > 0, ref, np.nan)
    rel = np.nan_to_num(rel, nan=0.0)

    scaled_spread = float(np.nanmax(rel))
    i_early = max(1, int(early_frac * n))
    # Attribution uses the window MEAN, not the max.  The early and late windows hold different numbers of samples, and the maximum of a noisy sequence grows with sample count, so a max-vs-max comparison would systematically favour whichever window is longer.
    early_spread = float(np.nanmean(rel[:, :i_early]))
    late_spread = float(np.nanmean(rel[:, i_early:])) if i_early < n else 0.0

    rates = np.array([r.growth_rate for r in runs], dtype=float)
    finite = rates[np.isfinite(rates)]
    if finite.size:
        rate_mean = float(finite.mean())
        rate_spread = float(finite.max() - finite.min())
    else:
        rate_mean, rate_spread = float("nan"), float("nan")

    rate_ok = (np.isfinite(rate_spread) and np.isfinite(rate_mean)
               and abs(rate_mean) > 0
               and rate_spread <= rate_tol * abs(rate_mean))
    scaled_ok = scaled_spread <= scaled_tol

    if scaled_ok and rate_ok:
        verdict = "linear"
        detail = ("Q' scales proportionally with eps and the growth rate is "
                  "amplitude-independent: the snapshots represent the "
                  "linearised operator.")
    elif late_spread > early_spread:
        verdict = "nonlinear"
        detail = ("Collapse degrades at LATE time (late spread "
                  f"{late_spread:.2%} > early {early_spread:.2%}): the "
                  "O(eps^2) term is active.  REDUCE the perturbation "
                  "amplitude and re-run.")
    else:
        verdict = "round-off"
        detail = ("Collapse degrades at EARLY time (early spread "
                  f"{early_spread:.2%} >= late {late_spread:.2%}): the "
                  "perturbation is near machine precision.  INCREASE the "
                  "amplitude, or tighten solver tolerances "
                  "(mac.mac_abs_tol, mg.bot_atol).")

    return LinearityReport(
        runs=runs,
        scaled_spread=scaled_spread,
        growth_rate_spread=rate_spread,
        growth_rate_mean=rate_mean,
        verdict=verdict,
        detail=detail,
    )


def recommend_amplitude(report):
    """Suggest the amplitude to use for production runs.

    Picks the largest amplitude that is still linear, because that maximises the signal-to-round-off margin without entering the nonlinear regime.
    """
    if report.verdict == "linear":
        return max(r.eps for r in report.runs)
    if report.verdict == "nonlinear":
        return min(r.eps for r in report.runs) / 10.0
    return max(r.eps for r in report.runs) * 10.0
