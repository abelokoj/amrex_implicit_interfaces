"""
NS-MFP snapshot assembly for the Rayleigh-Plateau capillary instability.

Implements the subspace-generation half of Ranjan, Unnikrishnan & Gaitonde (JCP 403, 2020) on top of the `amrex_implicit_interfaces` two-phase solver, without modifying the solver.

The B_f body force, and why the solver requires no modification
--------------------------------------------------------------
Ranjan et al. add a restraining body force

    B_f = dQbar/dt - F(Qbar)                                    (their Eq. 3)

so that base-state drift is cancelled and the marched perturbed solution linearises about a *fixed* Qbar.  Their Eq. (6) shows the role of B_f is exactly to annihilate the bracketed term

    [ dQbar/dt - F(Qbar) - B_f ] = 0.

We obtain the same cancellation without a body force by running an unperturbed *twin* alongside the perturbed case from the identical restart, with identical discretisation, time step and grid, and subtracting them snapshot by snapshot:

    Q'(t_n) = Q_perturbed(t_n) - Q_unperturbed(t_n).

Both runs experience the same base-state drift, so the difference removes it to machine precision, and Q' obeys the linearised dynamics to O(eps^2) exactly as in their Eq. (7)-(8).  This is the same "subtract the base flow from the total flow" step drawn in their Fig. 1; the drift that B_f exists to cancel is here cancelled by a computed twin rather than by a stored forcing term.  The cost is one extra base-state march, which is cheap relative to the perturbed sweep, and the benefit is that no C++/Fortran modification of the solver is required.

Two consequences worth stating explicitly, because they matter for the thesis:

1. This is valid only if the twin and perturbed runs are bitwise-reproducibly
discretised - same grid, same fixed_dt, same AMR behaviour, same MPI decomposition.  Use `ns.fixed_dt`, `amr.max_level=0` (or a frozen regrid schedule) and an identical processor count for both runs. `check_run_compatibility` enforces the metadata part of this.

2. Unlike the single-phase compressible case of Ranjan et al., the dominant
Rayleigh-Plateau perturbation lives on the *interface*, so the perturbation must be applied to the level set, not only to a state variable.  The linearisation of the surface-tension/curvature term is what makes the multiphase extension non-trivial; it is exercised implicitly here because the nonlinear solver applies the full curvature operator to the perturbed interface.

Reading plotfiles
-----------------
The runs must be configured with `ns.visual_nddata_format=1`, which makes the solver call AMReX `WriteMultiLevelPlotfile` and emit standard `plt*` files. Those are read with yt's AMReX frontend rather than a hand-rolled parser.
"""

# ORIENTATION
#
# This module builds the snapshot matrix, which is the input the decomposition consumes. Two ideas are worth understanding before reading the code.
#
# First, why there are two calculations. Linear stability analysis asks how a small disturbance evolves on top of a steady base flow. Here the base flow is not perfectly steady: the numerical scheme makes the unperturbed column drift slowly. The published method removes that drift by adding a body force to the equations, which means modifying the solver. This project instead runs the solver twice, once perturbed and once not, and subtracts one from the other snapshot by snapshot. Whatever drift both share cancels exactly, and the solver is never touched. That subtraction is the core of what this module does.
#
# Second, why a covering grid. The subtraction requires that the two runs have the same array shape at the same instant, and that every snapshot has the same shape as every other, because the snapshots become columns of one matrix. Reading each plotfile onto a uniform grid guarantees that. With adaptive refinement switched off, which is what the input decks do, the uniform grid is the data rather than an interpolation of it.


from __future__ import annotations

import glob
import os
import re
from dataclasses import dataclass

import numpy as np

# Rayleigh-Plateau runs are axisymmetric (geometry.coord_sys=1), so the two coordinate directions in a 2-D plotfile are (r, z).
DEFAULT_FIELDS = ("x_velocity", "y_velocity")


def _require_yt():
    try:
        import yt
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ImportError(
            "Reading AMReX plotfiles requires yt.  Install with:\n"
            "    pip install yt\n"
            "and ensure the runs used ns.visual_nddata_format=1."
        ) from exc
    return yt


def find_plotfiles(directory, prefix="nddataPLT"):
    """Return plotfile paths under `directory`, sorted by step number.

    Sorting is numeric on the trailing digits, so plt00009 precedes plt00010 (a plain lexicographic sort would also work for zero-padded names, but not for unpadded ones).

    NOTE ON THE PREFIX (confirmed by running the solver): the solver does NOT honour `amr.plot_file`.  With ns.visual_nddata_format=1 it writes directories named `nddataPLT00000023` (plus `MOF_PLT*` and `.tec` files), not `plt00023`.  The default here is therefore `nddataPLT`; pass `prefix="plt"` only if your build behaves differently.
    """
    pattern = os.path.join(directory, prefix + "[0-9]*")
    paths = [p for p in glob.glob(pattern) if os.path.isdir(p)]

    def step_of(p):
        m = re.search(r"(\d+)$", os.path.basename(p))
        return int(m.group(1)) if m else -1

    return sorted(paths, key=step_of)


def plotfile_step(path):
    """Integer step number encoded in a plotfile name."""
    m = re.search(r"(\d+)$", os.path.basename(path))
    if not m:
        raise ValueError(f"cannot parse a step number from {path!r}")
    return int(m.group(1))


@dataclass
class PlotfileData:
    """A single plotfile flattened onto a uniform level-0 grid."""

    path: str
    time: float
    step: int
    fields: tuple
    data: np.ndarray          # 1-D, fields concatenated
    shape: tuple              # per-field grid shape
    domain_left: np.ndarray
    domain_right: np.ndarray


def read_plotfile(path, fields=DEFAULT_FIELDS, level=0):
    """Read `fields` from an AMReX plotfile onto a uniform covering grid.

    A covering grid at `level` is used so every snapshot has an identical, fixed-size state vector - a hard requirement for assembling a snapshot matrix.  With `amr.max_level=0` (recommended for LSA, and what the reference deck uses) this is exact rather than an interpolation.
    """
    yt = _require_yt()
    ds = yt.load(path)

    dims = ds.domain_dimensions.copy()
    # yt pads 2-D datasets to 3 dimensions with a single cell in z.
    if ds.dimensionality == 2:
        dims[2] = 1

    ref = int(ds.refine_by) ** int(level)
    grid_dims = dims * ref
    cg = ds.covering_grid(level=level, left_edge=ds.domain_left_edge,
                          dims=grid_dims)

    arrays = []
    for f in fields:
        arr = np.asarray(cg[f]).astype(np.float64)
        arrays.append(arr.ravel(order="C"))

    shape = tuple(int(d) for d in grid_dims)
    return PlotfileData(
        path=path,
        time=float(ds.current_time),
        step=plotfile_step(path),
        fields=tuple(fields),
        data=np.concatenate(arrays),
        shape=shape,
        domain_left=np.asarray(ds.domain_left_edge, dtype=float),
        domain_right=np.asarray(ds.domain_right_edge, dtype=float),
    )


def check_run_compatibility(a: PlotfileData, b: PlotfileData, time_tol=1e-9):
    """Raise unless two plotfiles may legitimately be subtracted.

    Guards the silent-garbage failure mode where the perturbed and twin runs drifted onto different grids, times or field sets and the difference is meaningless rather than a linearised perturbation.
    """
    if a.fields != b.fields:
        raise ValueError(f"field mismatch: {a.fields} vs {b.fields}")
    if a.shape != b.shape:
        raise ValueError(
            f"grid mismatch: {a.shape} ({a.path}) vs {b.shape} ({b.path}). "
            "Perturbed and unperturbed runs must use the same grid; set "
            "amr.max_level=0 or freeze the regrid schedule."
        )
    if a.data.size != b.data.size:
        raise ValueError(f"state size mismatch: {a.data.size} vs {b.data.size}")
    if not np.allclose(a.domain_left, b.domain_left) or not np.allclose(
            a.domain_right, b.domain_right):
        raise ValueError("domain extent mismatch between the two runs")
    scale = max(1.0, abs(a.time), abs(b.time))
    if abs(a.time - b.time) > time_tol * scale:
        raise ValueError(
            f"time mismatch: {a.time} ({a.path}) vs {b.time} ({b.path}). "
            "Both runs must use the same ns.fixed_dt and plot cadence."
        )


def build_snapshot_matrix(perturbed_dir, unperturbed_dir,
                          fields=DEFAULT_FIELDS, level=0,
                          skip=0, stride=1, max_snapshots=None,
                          verbose=True, prefix="nddataPLT"):
    """Assemble the NS-MFP subspace  Q'(t_n) = Q_pert(t_n) - Q_unpert(t_n).

    Parameters
    ----------
    perturbed_dir, unperturbed_dir : str
        Directories holding the `plt*` output of the two twin runs.
    fields : sequence of str
        Which plotfile variables enter the state vector.  Ranjan et al. note
        (section 2.2) that any subset may be used since DMD is decoupled from
        subspace generation, and that using one well-chosen variable can cut
        memory up to 80%.  For Rayleigh-Plateau the velocity components carry
        the capillary mode cleanly.
    skip : int
        Discard this many leading snapshots to let the initial impulse shed
        its non-modal transient before the subspace is formed.
    stride : int
        Use every `stride`-th snapshot; multiplies the effective sampling dt.
    max_snapshots : int, optional
        Cap on the number of snapshots retained.

    Returns
    -------
    (X, dt, times) : (n_space, n_snap) array, float, ndarray
    """
    p_paths = find_plotfiles(perturbed_dir, prefix=prefix)
    u_paths = find_plotfiles(unperturbed_dir, prefix=prefix)
    if not p_paths:
        raise FileNotFoundError(f"no nddataPLT* directories under {perturbed_dir!r}")
    if not u_paths:
        raise FileNotFoundError(f"no nddataPLT* directories under {unperturbed_dir!r}")

    # Pair strictly by step number: a missing plotfile on one side must not silently shift the pairing.
    u_by_step = {plotfile_step(p): p for p in u_paths}
    pairs = [(p, u_by_step[plotfile_step(p)])
             for p in p_paths if plotfile_step(p) in u_by_step]
    if len(pairs) < 2:
        raise ValueError(
            f"only {len(pairs)} matching step(s) between the two runs; "
            "need at least 2.  Check that both used the same amr.plot_int."
        )

    pairs = pairs[skip::stride]
    if max_snapshots is not None:
        pairs = pairs[:max_snapshots]
    if len(pairs) < 2:
        raise ValueError("fewer than 2 snapshots remain after skip/stride")

    cols, times = [], []
    for i, (pp, up) in enumerate(pairs):
        pd = read_plotfile(pp, fields=fields, level=level)
        ud = read_plotfile(up, fields=fields, level=level)
        check_run_compatibility(pd, ud)
        cols.append(pd.data - ud.data)
        times.append(pd.time)
        if verbose and (i % 25 == 0 or i == len(pairs) - 1):
            print(f"  [{i + 1}/{len(pairs)}] t={pd.time:.6g} "
                  f"|Q'|={np.linalg.norm(cols[-1]):.6e}")

    X = np.column_stack(cols)
    times = np.asarray(times)
    dt = _uniform_dt(times)

    if verbose:
        print(f"subspace: {X.shape[0]} dof x {X.shape[1]} snapshots, dt={dt:g}")
    return X, dt, times


def _uniform_dt(times, rtol=1e-6):
    """Sampling interval, insisting the snapshots are evenly spaced.

    DMD assumes a constant dt; non-uniform sampling would bias every eigenvalue, so this is checked rather than averaged over silently.
    """
    d = np.diff(times)
    if d.size == 0:
        raise ValueError("need at least 2 snapshot times")
    if np.any(d <= 0):
        raise ValueError("snapshot times are not strictly increasing")
    dt = float(d.mean())
    if np.max(np.abs(d - dt)) > rtol * max(dt, 1e-300):
        raise ValueError(
            f"snapshot spacing is not uniform (spread "
            f"{np.max(np.abs(d - dt)):.3e} vs dt {dt:.3e}); DMD requires a "
            "constant sampling interval."
        )
    return dt


def perturbation_norms(perturbed_dir, unperturbed_dir,
                       fields=DEFAULT_FIELDS, level=0, verbose=False,
                       prefix="nddataPLT"):
    """L2 norm of Q'(t) versus time.

    On a linear, exponentially growing mode this is a straight line on a log plot with slope equal to the growth rate, so it is both a sanity check on the DMD result and the raw material for the proportionality (linearity) test.
    """
    X, dt, times = build_snapshot_matrix(
        perturbed_dir, unperturbed_dir, fields=fields, level=level,
        verbose=verbose, prefix=prefix)
    return times, np.linalg.norm(X, axis=0)


def fit_exponential_growth(times, norms, t_min=None, t_max=None):
    """Least-squares slope of log|Q'| versus t over an optional window.

    Returns (growth_rate, intercept, r_squared).  A poor r_squared is the usual signal that the window still contains the initial transient, or that nonlinearity has set in at the late end.
    """
    times = np.asarray(times, dtype=float)
    norms = np.asarray(norms, dtype=float)
    m = np.isfinite(norms) & (norms > 0)
    if t_min is not None:
        m &= times >= t_min
    if t_max is not None:
        m &= times <= t_max
    if m.sum() < 2:
        raise ValueError("need at least 2 positive samples in the fit window")

    t, y = times[m], np.log(norms[m])
    slope, intercept = np.polyfit(t, y, 1)
    pred = slope * t + intercept
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(slope), float(intercept), r2
