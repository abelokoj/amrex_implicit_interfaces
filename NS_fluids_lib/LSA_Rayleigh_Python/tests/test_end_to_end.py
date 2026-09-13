"""End-to-end test: real AMReX plotfiles -> twin subtraction -> DMD -> sigma.

Synthetic plotfiles are written in the genuine on-disk AMReX format and read back through the production code path (nsmfp.read_plotfile -> yt -> covering_grid).  Nothing here is mocked, so a format or ordering mistake in the reader shows up as a wrong growth rate rather than passing silently.

The planted signal mimics the real experiment:

    Q_unperturbed(t) = drift(t)                      (parasitic currents)
    Q_perturbed(t)   = drift(t) + eps e^{sigma t} phi(z)

so the twin subtraction must remove `drift` exactly and leave a clean exponentially growing capillary mode whose rate DMD recovers.  `drift` is made deliberately large compared with the perturbation - as spurious currents often are in a surface-tension calculation - so the test fails if the subtraction is skipped or mispaired.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
sys.path.insert(0, os.path.dirname(__file__))

import dispersion  # noqa: E402
import dmd as dmd_mod  # noqa: E402
import nsmfp  # noqa: E402
import synthetic_plotfile as sp  # noqa: E402

yt = pytest.importorskip("yt")
yt.set_log_level(50)

NX, NY = 10, 32
KR0 = 0.7
R0 = 1.0
LAMBDA = 2.0 * np.pi / (KR0 / R0)
SIGMA_EXACT = float(dispersion.growth_rate(KR0 / R0, r0=R0))
EPS = 1e-3
DT = 0.05
NSNAP = 60


def _z(ny=NY):
    return (np.arange(ny) + 0.5) * LAMBDA / ny


def _mode_shape(nx=NX, ny=NY):
    """Interface-localised sinusoid: cos(k z), tapered in r."""
    z = _z(ny)
    r = (np.arange(nx) + 0.5) * (4.0 * R0) / nx
    taper = np.exp(-((r - R0) ** 2) / 0.25)
    return np.outer(taper, np.cos(KR0 * z / R0))


def _drift(t, nx=NX, ny=NY):
    """Base-state drift (parasitic currents): large, slow, non-modal."""
    z = _z(ny)
    r = (np.arange(nx) + 0.5) * (4.0 * R0) / nx
    return 0.5 * np.outer(np.sin(3.0 * r), np.cos(2.0 * z)) * (1.0 + 0.1 * t)


@pytest.fixture(scope="module")
def twin_runs(tmp_path_factory):
    """Write an unperturbed twin and a perturbed run as real plotfiles."""
    root = tmp_path_factory.mktemp("runs")
    base_dir = str(root / "base")
    pert_dir = str(root / "pert")

    steps = list(range(0, NSNAP * 10, 10))
    times = [i * DT for i in range(NSNAP)]
    shape_x = _mode_shape()
    shape_y = 0.6 * _mode_shape()

    def base_fn(t, name, nx, ny):
        return _drift(t, nx, ny)

    def pert_fn(t, name, nx, ny):
        growth = EPS * np.exp(SIGMA_EXACT * t)
        shape = shape_x if name == "x_velocity" else shape_y
        return _drift(t, nx, ny) + growth * shape

    sp.write_series(base_dir, base_fn, steps, times, shape=(NX, NY),
                    domain_right=(4.0 * R0, LAMBDA))
    sp.write_series(pert_dir, pert_fn, steps, times, shape=(NX, NY),
                    domain_right=(4.0 * R0, LAMBDA))
    return pert_dir, base_dir


def test_find_plotfiles_sorts_numerically(twin_runs):
    pert, _ = twin_runs
    paths = nsmfp.find_plotfiles(pert)
    assert len(paths) == NSNAP
    steps = [nsmfp.plotfile_step(p) for p in paths]
    assert steps == sorted(steps)


def test_read_plotfile_shape_and_time(twin_runs):
    pert, _ = twin_runs
    pd = nsmfp.read_plotfile(nsmfp.find_plotfiles(pert)[3])
    assert pd.data.size == NX * NY * 2          # two fields concatenated
    assert pd.time == pytest.approx(3 * DT)
    assert pd.step == 30


def test_twin_subtraction_removes_the_drift(twin_runs):
    """|Q'| must reflect only the planted perturbation, not the drift."""
    pert, base = twin_runs
    X, dt, times = nsmfp.build_snapshot_matrix(pert, base, verbose=False)
    assert dt == pytest.approx(DT)
    assert X.shape[1] == NSNAP

    expected0 = EPS * np.linalg.norm(
        np.concatenate([_mode_shape().ravel(), 0.6 * _mode_shape().ravel()]))
    assert np.linalg.norm(X[:, 0]) == pytest.approx(expected0, rel=1e-10)

    # the drift itself is far larger than the perturbation, so this is a meaningful check that it really cancelled
    drift_norm = np.linalg.norm(_drift(0.0).ravel())
    assert drift_norm > 100 * np.linalg.norm(X[:, 0])


def test_end_to_end_growth_rate_matches_dispersion_relation(twin_runs):
    """The headline check: recovered sigma equals the exact Rayleigh value."""
    pert, base = twin_runs
    X, dt, _ = nsmfp.build_snapshot_matrix(pert, base, verbose=False)
    res = dmd_mod.dmd(X, dt)
    measured = res.leading_growth_rate()
    assert measured == pytest.approx(SIGMA_EXACT, rel=1e-6)
    assert 0.3 < measured < 0.4          # the physically expected magnitude


def test_end_to_end_norm_fit_agrees_with_dmd(twin_runs):
    pert, base = twin_runs
    times, norms = nsmfp.perturbation_norms(pert, base)
    rate, _, r2 = nsmfp.fit_exponential_growth(times, norms)
    assert rate == pytest.approx(SIGMA_EXACT, rel=1e-8)
    assert r2 > 0.999999


def test_mispaired_steps_are_rejected(twin_runs, tmp_path):
    """Deleting a twin plotfile must not silently shift the pairing."""
    pert, base = twin_runs
    import shutil

    trimmed = str(tmp_path / "base_trimmed")
    shutil.copytree(base, trimmed)
    # remove one interior snapshot from the twin
    victim = nsmfp.find_plotfiles(trimmed)[5]
    shutil.rmtree(victim)

    # Pairing is by step number, so the survivors still align; the missing step is dropped, which makes the sampling non-uniform and must raise.
    with pytest.raises(ValueError, match="not uniform"):
        nsmfp.build_snapshot_matrix(pert, trimmed, verbose=False)


def test_grid_mismatch_is_rejected(tmp_path):
    """Twins on different grids must not be subtracted."""
    a = str(tmp_path / "a")
    b = str(tmp_path / "b")
    sp.write_series(a, lambda t, n, nx, ny: np.ones((nx, ny)),
                    [0, 10], [0.0, DT], shape=(8, 16))
    sp.write_series(b, lambda t, n, nx, ny: np.ones((nx, ny)),
                    [0, 10], [0.0, DT], shape=(8, 32))
    with pytest.raises(ValueError, match="grid mismatch"):
        nsmfp.build_snapshot_matrix(a, b, verbose=False)


def test_time_mismatch_is_rejected(tmp_path):
    """Same steps but different times means the decks disagreed on dt."""
    a = str(tmp_path / "a")
    b = str(tmp_path / "b")
    sp.write_series(a, lambda t, n, nx, ny: np.ones((nx, ny)),
                    [0, 10], [0.0, 0.05], shape=(8, 16))
    sp.write_series(b, lambda t, n, nx, ny: np.ones((nx, ny)),
                    [0, 10], [0.0, 0.07], shape=(8, 16))
    with pytest.raises(ValueError, match="time mismatch"):
        nsmfp.build_snapshot_matrix(a, b, verbose=False)


def test_missing_directory_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        nsmfp.build_snapshot_matrix(str(tmp_path / "nope"),
                                    str(tmp_path / "also_nope"), verbose=False)


def test_skip_and_stride(twin_runs):
    """Striding must scale dt, and still give the same growth rate."""
    pert, base = twin_runs
    X, dt, _ = nsmfp.build_snapshot_matrix(pert, base, skip=10, stride=3,
                                           verbose=False)
    assert dt == pytest.approx(3 * DT)
    assert X.shape[1] == len(range(10, NSNAP, 3))
    assert dmd_mod.dmd(X, dt).leading_growth_rate() == pytest.approx(
        SIGMA_EXACT, rel=1e-6)


def test_single_field_subspace_still_works(twin_runs):
    """Ranjan et al. note one well-chosen variable can suffice (and cuts
    memory up to 80%).  Check the pipeline honours a single-field subspace."""
    pert, base = twin_runs
    X, dt, _ = nsmfp.build_snapshot_matrix(pert, base, fields=("x_velocity",),
                                           verbose=False)
    assert X.shape[0] == NX * NY
    assert dmd_mod.dmd(X, dt).leading_growth_rate() == pytest.approx(
        SIGMA_EXACT, rel=1e-6)


def test_identical_runs_give_zero_perturbation(twin_runs):
    """Subtracting a run from itself must give exactly zero - the guard the
    driver uses to detect a perturbed deck that forgot to set radblob."""
    _, base = twin_runs
    X, _, _ = nsmfp.build_snapshot_matrix(base, base, verbose=False)
    assert np.allclose(X, 0.0)
