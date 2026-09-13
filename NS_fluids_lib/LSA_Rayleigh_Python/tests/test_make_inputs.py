"""Tests for input-deck generation.

The critical invariant is the twin invariant: a perturbed deck and its unperturbed twin must differ in exactly one functional setting (`ns.radblob`).  If anything else differs, the snapshot subtraction that stands in for the NS-MFP body force B_f is no longer a linearised perturbation, and the whole pipeline silently produces a wrong answer.
"""

import math
import os
import re
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

import make_inputs  # noqa: E402
import dispersion  # noqa: E402


def functional_lines(text):
    """Deck lines that actually configure the solver (drop comments/blanks)."""
    out = []
    for line in text.splitlines():
        s = line.split("#", 1)[0].strip()
        if s:
            out.append(s)
    return out


def get_param(text, key):
    for line in functional_lines(text):
        m = re.match(rf"^{re.escape(key)}\s*=\s*(.+)$", line)
        if m:
            return m.group(1).strip()
    raise KeyError(key)


@pytest.fixture
def twin_pair():
    base = make_inputs.plan_run(0.7, 0.0, "twin")
    pert = make_inputs.plan_run(0.7, 1e-3, "perturbed")
    return make_inputs.render(base), make_inputs.render(pert)


def test_twin_invariant_only_radblob_differs(twin_pair):
    base_txt, pert_txt = twin_pair
    b, p = functional_lines(base_txt), functional_lines(pert_txt)
    assert len(b) == len(p), "decks have different numbers of settings"
    diffs = [(x, y) for x, y in zip(b, p) if x != y]
    assert len(diffs) == 1, f"expected exactly one functional difference, got {diffs}"
    assert diffs[0][0].startswith("ns.radblob")
    assert diffs[0][1].startswith("ns.radblob")


def test_radblob_carries_the_amplitude(twin_pair):
    base_txt, pert_txt = twin_pair
    assert float(get_param(base_txt, "ns.radblob")) == 0.0
    assert float(get_param(pert_txt, "ns.radblob")) == pytest.approx(1e-3)


def test_no_unsubstituted_tokens(twin_pair):
    for txt in twin_pair:
        assert "@" not in txt, "template token left unsubstituted"


def test_wavelength_matches_domain_and_yblob():
    """yblob (perturbation wavelength) must equal the axial domain height,
    or the imposed sinusoid is not periodic in the box."""
    for kr0 in (0.2, 0.5, 0.7, 0.9):
        txt = make_inputs.render(make_inputs.plan_run(kr0, 1e-3, "p"))
        yblob = float(get_param(txt, "ns.yblob"))
        hi = get_param(txt, "geometry.prob_hi").split()
        assert float(hi[1]) == pytest.approx(yblob, rel=1e-9)
        assert yblob == pytest.approx(2 * math.pi / kr0, rel=1e-9)


def test_axial_direction_is_periodic():
    txt = make_inputs.render(make_inputs.plan_run(0.7, 1e-3, "p"))
    assert get_param(txt, "geometry.is_periodic").split() == ["0", "1"]
    # z BCs must be interior/periodic (0) to match
    assert get_param(txt, "ns.lo_bc").split()[1] == "0"
    assert get_param(txt, "ns.hi_bc").split()[1] == "0"


def test_probtype_and_axis_dir_are_the_capillary_case():
    txt = make_inputs.render(make_inputs.plan_run(0.7, 1e-3, "p"))
    assert get_param(txt, "ns.probtype") == "41"
    assert get_param(txt, "ns.axis_dir") == "4"
    assert get_param(txt, "geometry.coord_sys") == "1"   # RZ


def test_plotfile_format_is_amrex_plt():
    """visual_nddata_format must be 1, else no plt* files are written and
    the Python post-processing has nothing to read."""
    txt = make_inputs.render(make_inputs.plan_run(0.7, 1e-3, "p"))
    assert get_param(txt, "ns.visual_nddata_format") == "1"


def test_builtin_power_iteration_lsa_is_disabled():
    """NS-MFP does not use the solver's own Krylov/power-iteration driver."""
    txt = make_inputs.render(make_inputs.plan_run(0.7, 1e-3, "p"))
    assert get_param(txt, "amr.LSA_activate") == "0"
    assert get_param(txt, "amr.LSA_nsteps_krylov_subspace_method") == "0"


def test_fixed_dt_is_set_and_below_capillary_limit():
    for kr0 in (0.2, 0.7, 0.9):
        plan = make_inputs.plan_run(kr0, 1e-3, "p")
        assert plan["fixed_dt"] < plan["dt_cap"]
        txt = make_inputs.render(plan)
        assert float(get_param(txt, "ns.fixed_dt")) == pytest.approx(
            plan["fixed_dt"], rel=1e-9)


def test_single_amr_level_for_grid_identity():
    txt = make_inputs.render(make_inputs.plan_run(0.7, 1e-3, "p"))
    assert get_param(txt, "amr.max_level") == "0"


def test_max_step_is_multiple_of_plot_int():
    """Both twins must end on a written plotfile at the same step."""
    for kr0 in (0.2, 0.5, 0.7, 0.9):
        p = make_inputs.plan_run(kr0, 1e-3, "p")
        assert p["max_step"] % p["plot_int"] == 0
        assert p["n_snapshots"] >= 2


def test_growth_rate_annotation_matches_dispersion_module():
    for kr0 in (0.3, 0.7):
        p = make_inputs.plan_run(kr0, 1e-3, "p")
        assert p["sigma_exact"] == pytest.approx(
            float(dispersion.growth_rate(kr0)), rel=1e-12)


def test_stable_band_still_produces_a_finite_run():
    """k r0 > 1 is stable: no e-folding time exists, so the planner must
    fall back on the capillary time scale rather than dividing by ~0."""
    p = make_inputs.plan_run(1.3, 1e-3, "p")
    assert math.isfinite(p["stop_time"]) and p["stop_time"] > 0
    assert p["max_step"] > 0
    assert p["sigma_exact"] == 0.0


def test_cell_counts_are_block_friendly():
    for kr0 in (0.2, 0.45, 0.7, 0.95):
        p = make_inputs.plan_run(kr0, 1e-3, "p")
        assert p["nr"] % 8 == 0 and p["nz"] % 8 == 0


def test_cells_are_near_square():
    """Radial and axial spacing should be comparable, so the interface is
    resolved isotropically."""
    for kr0 in (0.2, 0.5, 0.7, 0.9):
        p = make_inputs.plan_run(kr0, 1e-3, "p")
        assert 0.5 < p["dr"] / p["dz"] < 2.0


def test_resolution_option_is_honoured():
    coarse = make_inputs.plan_run(0.7, 1e-3, "p", cells_per_r0=8)
    fine = make_inputs.plan_run(0.7, 1e-3, "p", cells_per_r0=32)
    assert fine["nr"] > coarse["nr"]
    assert fine["dr"] < coarse["dr"]


def test_write_deck_round_trip(tmp_path):
    plan = make_inputs.plan_run(0.7, 1e-3, "p")
    path = make_inputs.write_deck(plan, str(tmp_path), "inputs.test")
    assert os.path.exists(path)
    with open(path) as f:
        assert get_param(f.read(), "ns.probtype") == "41"


def test_cli_sweep_eps_writes_one_shared_twin(tmp_path, capsys):
    make_inputs.main(["sweep-eps", "--kr0", "0.7",
                      "--eps", "1e-4", "1e-3",
                      "--outdir", str(tmp_path)])
    names = sorted(os.listdir(tmp_path))
    assert "inputs.k0.7.base" in names
    assert "inputs.k0.7.eps0.0001" in names
    assert "inputs.k0.7.eps0.001" in names
    # exactly one base deck is generated regardless of amplitude count
    assert sum(n.endswith(".base") for n in names) == 1


def test_cli_sweep_k_writes_a_twin_per_wavenumber(tmp_path):
    make_inputs.main(["sweep-k", "--kr0", "0.4", "0.7",
                      "--eps", "1e-3", "--outdir", str(tmp_path)])
    names = sorted(os.listdir(tmp_path))
    assert sum(n.endswith(".base") for n in names) == 2
