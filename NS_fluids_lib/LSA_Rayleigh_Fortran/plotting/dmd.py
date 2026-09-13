"""
Dynamic Mode Decomposition for stability-mode extraction.

Implements the SVD-based ("exact"/projected) DMD of

    Schmid, P.J. (2010) "Dynamic mode decomposition of numerical and
    experimental data", J. Fluid Mech. 656, 5-28,

which is the variant Ranjan, Unnikrishnan & Gaitonde (JCP 403, 2020) use in section 2.2 to extract stability modes from an NS-MFP snapshot subspace, with the amplitude ranking of Jovanovic, Schmid & Nichols (2014).

Conventions (following Ranjan et al. section 2.2)
------------------------------------------------
Given N snapshots separated by a constant sampling interval dt, DMD returns Ritz values `lam`.  The continuous-time eigenvalue of the Jacobian A is

    omega_c = log(lam) / dt,          (complex)

whose real part is the **growth rate** and whose imaginary part is the **circular frequency**:

    growth_rate    sigma   = Re(omega_c)
    circular freq  omega_c = Im(omega_c)
    linear freq    omega   = Im(omega_c) / (2 pi)

Note Ranjan et al. warn that this is "the reverse of the convention usually followed in the classical stability analysis"; here `growth_rate` is unambiguously the exponential rate e^{sigma t}, which is what the Rayleigh-Plateau dispersion relation predicts.

For the Rayleigh-Plateau base state the physically relevant leading mode is *stationary* (purely real eigenvalue, omega = 0): the capillary instability grows monotonically rather than oscillating.  Helpers `stationary_modes` and `leading_growth_rate` exist for exactly this case.
"""

# ORIENTATION
#
# Dynamic Mode Decomposition, which is the step that turns a sequence of flow fields into growth rates.
#
# The idea in one sentence: assume each snapshot is obtained from the previous one by multiplying by a fixed matrix A, then find the eigenvalues of A without ever forming it. If x_{n+1} = A x_n, then an eigenvalue lambda of A with |lambda| > 1 is a mode that grows from one snapshot to the next, and sigma = log(lambda)/dt converts that per-snapshot factor into a growth rate per unit time.
#
# The obstacle is size. Each snapshot here has thousands of entries, so A would have millions, and forming it is neither possible nor necessary. The standard trick, used below, is to project onto the handful of directions the snapshots actually span, found by a singular value decomposition, and to take the eigenvalues of the small projected matrix instead. Those eigenvalues are the ones we want; the large matrix never appears. This is what "Jacobian-free" means in this project.


from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

import report


@dataclass
class DMDResult:
    """Outcome of a DMD.  Arrays are ordered by descending `amplitude`."""

    ritz: np.ndarray            # discrete-time eigenvalues lambda
    omega: np.ndarray           # continuous-time eigenvalues log(lam)/dt
    modes: np.ndarray           # (n_space, n_modes) complex mode shapes
    amplitude: np.ndarray       # |b|, contribution to the first snapshot
    dt: float
    singular_values: np.ndarray = field(default=None, repr=False)
    rank: int = 0

    @property
    def growth_rate(self) -> np.ndarray:
        """sigma = Re(log(lambda)/dt); e^{sigma t} amplification."""
        return self.omega.real

    @property
    def circular_frequency(self) -> np.ndarray:
        """omega_c = Im(log(lambda)/dt)."""
        return self.omega.imag

    @property
    def frequency(self) -> np.ndarray:
        """Linear frequency omega = omega_c / (2 pi) (Strouhal number if the
        snapshots are non-dimensionalised by characteristic length/velocity)."""
        return self.omega.imag / (2.0 * np.pi)

    def sort_by(self, key="amplitude"):
        """Return index order: 'amplitude', 'growth' (descending sigma), or
        'magnitude' (descending |lambda|)."""
        if key == "amplitude":
            return np.argsort(-self.amplitude)
        if key == "growth":
            return np.argsort(-self.growth_rate)
        if key == "magnitude":
            return np.argsort(-np.abs(self.ritz))
        raise ValueError(f"unknown sort key {key!r}")

    def stationary_modes(self, freq_tol=1e-6, amp_frac=0.0):
        """Indices of non-oscillatory modes (|omega_c| <= freq_tol).

        The Rayleigh-Plateau instability is stationary, so the mode of interest lives here.  `amp_frac` optionally drops modes whose amplitude is below that fraction of the largest amplitude.
        """
        keep = np.abs(self.circular_frequency) <= freq_tol
        if amp_frac > 0.0 and self.amplitude.size:
            keep &= self.amplitude >= amp_frac * self.amplitude.max()
        return np.where(keep)[0]

    def leading_growth_rate(self, stationary_only=True, freq_tol=1e-6,
                            amp_frac=1e-3):
        """Largest growth rate among physically retained modes.

        For Rayleigh-Plateau use the default (`stationary_only=True`): the capillary mode is non-oscillatory, and restricting to it avoids picking up weak spurious oscillatory pairs.  Returns NaN if nothing survives the filter.
        """
        if stationary_only:
            idx = self.stationary_modes(freq_tol=freq_tol, amp_frac=amp_frac)
        else:
            idx = np.arange(self.omega.size)
            if amp_frac > 0.0 and self.amplitude.size:
                idx = idx[self.amplitude[idx] >= amp_frac * self.amplitude.max()]
        if idx.size == 0:
            return float("nan")
        return float(np.max(self.growth_rate[idx]))

    def summary(self, n=10, sort="amplitude"):
        """Human-readable table of the leading modes."""
        order = self.sort_by(sort)[:n]
        rows = [[str(j), f"{self.growth_rate[i]:.6f}",
                 f"{self.frequency[i]:.6f}", f"{abs(self.ritz[i]):.6f}",
                 f"{self.amplitude[i]:.3e}"]
                for j, i in enumerate(order)]
        title = (f"Decomposition: rank {self.rank}, dt {self.dt:.3e}, "
                 f"{self.omega.size} modes")
        return report.table(
            rows,
            ["mode", "growth sigma", "freq omega", "|lambda|", "amplitude"],
            colalign=["right"] * 5, title=title)


def _truncation_rank(s, rank=None, tol=1e-10):
    """Resolve the SVD truncation rank from an explicit rank or a tolerance."""
    if rank is not None and rank > 0:
        return min(int(rank), s.size)
    if s.size == 0:
        return 0
    return max(int(np.sum(s > tol * s[0])), 1)


def dmd(snapshots, dt, rank=None, tol=1e-10, remove_mean=False,
        exact_modes=True):
    """SVD-based DMD of a snapshot sequence.

    Parameters
    ----------
    snapshots : (n_space, n_snap) array
        Column j is the state at time j*dt.  For NS-MFP these are the
        *perturbation* fields Q', i.e. the Jacobian-vector products.
    dt : float
        Sampling interval between consecutive columns.
    rank : int, optional
        SVD truncation rank.  If None, chosen from `tol`.
    tol : float
        Relative singular-value cutoff used when `rank` is None.
    remove_mean : bool
        Subtract the temporal mean before decomposing.  Leave False for
        stability analysis: the mean carries the stationary mode we want.
    exact_modes : bool
        True  -> exact DMD modes (Tu et al. 2014), Phi = Y V S^-1 w.
        False -> projected modes, Phi = U w.

    Returns
    -------
    DMDResult
    """
    X = np.asarray(snapshots)
    if X.ndim != 2:
        raise ValueError(f"snapshots must be 2-D (n_space, n_snap), got {X.shape}")
    if X.shape[1] < 2:
        raise ValueError("need at least 2 snapshots")
    if not np.all(np.isfinite(X)):
        raise ValueError("snapshots contain non-finite values")
    if dt <= 0:
        raise ValueError("dt must be positive")

    X = X.astype(np.complex128, copy=True)
    if remove_mean:
        X -= X.mean(axis=1, keepdims=True)

    A, B = X[:, :-1], X[:, 1:]

    U, s, Vh = np.linalg.svd(A, full_matrices=False)
    r = _truncation_rank(s, rank=rank, tol=tol)
    U_r, s_r, V_r = U[:, :r], s[:r], Vh[:r].conj().T

    # Companion matrix in the POD basis:  S = U* B V S^-1
    S = U_r.conj().T @ B @ V_r @ np.diag(1.0 / s_r)
    mu, w = np.linalg.eig(S)

    if exact_modes:
        modes = B @ V_r @ np.diag(1.0 / s_r) @ w
    else:
        modes = U_r @ w

    # Continuous-time eigenvalues.  A zero Ritz value would mean a mode that annihilates in one step; map it to -inf growth rather than raising.
    with np.errstate(divide="ignore", invalid="ignore"):
        omega = np.log(mu.astype(np.complex128)) / dt
    omega = np.where(mu == 0, -np.inf + 0j, omega)

    # Amplitudes: contribution of each mode to the first snapshot (Jovanovic et al. 2014 ranking, in its simplest projection form).
    b, *_ = np.linalg.lstsq(modes, X[:, 0], rcond=None)

    order = np.argsort(-np.abs(b))
    return DMDResult(
        ritz=mu[order],
        omega=omega[order],
        modes=modes[:, order],
        amplitude=np.abs(b)[order],
        dt=dt,
        singular_values=s,
        rank=r,
    )


def convergence_study(snapshots, dt, sizes=None, starts=(0,), **kw):
    """Repeat the DMD over varying subspace sizes and starting snapshots.

    Ranjan et al. establish convergence by "varying the size and identity of the subspace" (section 5); this automates that sweep.  Returns a list of dicts with the leading stationary growth rate for each configuration.
    """
    X = np.asarray(snapshots)
    n_snap = X.shape[1]
    if sizes is None:
        sizes = [n for n in (50, 100, 200, 400, 800, 1600, n_snap)
                 if 2 < n <= n_snap]
        sizes = sorted(set(sizes))

    out = []
    for start in starts:
        for n in sizes:
            if start + n > n_snap:
                continue
            res = dmd(X[:, start:start + n], dt, **kw)
            out.append({
                "start": start,
                "n_snapshots": n,
                "rank": res.rank,
                "growth_rate": res.leading_growth_rate(),
                "result": res,
            })
    return out
