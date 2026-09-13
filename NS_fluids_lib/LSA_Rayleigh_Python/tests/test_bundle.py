"""Tests for the portable figure-data bundle.

Three properties are checked. The writers and readers here round-trip, including the sentinel that marks an axial station where no interface was found. The tolerant real parser accepts the Fortran output form that omits the exponent letter, which is the one way a valid bundle can arrive unreadable. And a bundle written by the Fortran tree is read correctly here, which is the property that matters: the two implementations write the same format, and neither is trusted to validate its own output alone.

The Fortran interoperability test is skipped when the Fortran bundle is absent, so that the Python suite remains runnable on a machine with no Fortran compiler. Running `make test` in the Fortran tree writes the file this looks for.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "python"))

import bundle  # noqa: E402

#: Where `bin/test_bundle` of the Fortran tree leaves its output.
FORTRAN_BUNDLE = "/tmp/lsa_bundle_test"


def test_profile_round_trip(tmp_path):
    z = np.linspace(0.05, 8.9, 24)
    times = np.array([0.0, 1.5, 3.0])
    radius = np.vstack([1.0 + 0.02 * np.exp(0.3 * t) * np.cos(0.7 * z)
                        for t in times])

    path = str(tmp_path / "profiles.dat")
    bundle.write_profiles(path, 8.976, z, times, radius)
    got = bundle.read_profiles(path)

    assert got["wavelength"] == pytest.approx(8.976)
    assert got["z"] == pytest.approx(z)
    assert got["times"] == pytest.approx(times)
    assert got["radius"] == pytest.approx(radius)


def test_missing_interface_becomes_nan(tmp_path):
    """A station with no interface must leave a gap, not a value.

    The writer stores the sentinel and the reader returns NaN. Were the sentinel returned unchanged, a plot would draw a line to r = -1 and the figure would show an excursion that never happened.
    """
    z = np.linspace(0.1, 1.0, 5)
    radius = np.array([[1.0, 1.1, np.nan, 1.1, 1.0]])
    path = str(tmp_path / "profiles.dat")
    bundle.write_profiles(path, 1.0, z, [0.0], radius)

    with open(path) as fh:
        assert "-1.0000000000000000E+00" in fh.read()

    got = bundle.read_profiles(path)
    assert np.isnan(got["radius"][0, 2])
    assert got["radius"][0, 0] == pytest.approx(1.0)


def test_field_round_trip_preserves_index_order(tmp_path):
    """The radial index must vary fastest, and survive the round trip.

    The values are distinct in all three indices, so a transposed or mis-strided write cannot pass unnoticed; a shape-only check would.
    """
    r = np.arange(1.0, 6.0)
    z = np.arange(1.0, 9.0)
    times = np.array([0.0, 2.0])
    data = np.empty((2, r.size, z.size))
    for it in range(2):
        for i in range(r.size):
            for j in range(z.size):
                data[it, i, j] = 100 * (it + 1) + 10 * (j + 1) + (i + 1)

    path = str(tmp_path / "field_v.dat")
    bundle.write_field(path, "v", r, z, times, data)
    got = bundle.read_field(path)

    assert got["field"] == "v"
    assert got["data"].shape == (2, r.size, z.size)
    assert got["data"] == pytest.approx(data)


def test_tolerant_real_parses_the_fortran_form():
    """Fortran drops the exponent letter when the exponent needs three digits.

    The writers here use a three-digit descriptor and never emit that form, but a bundle written by an earlier build may carry it, and such a file is otherwise intact. Rejecting it would discard usable figure data.
    """
    assert bundle._real("4.6598606362008991-310") == pytest.approx(4.66e-310,
                                                                   rel=1e-6)
    assert bundle._real("-1.5000000000000000+002") == pytest.approx(-150.0)
    assert bundle._real("-1.25E+00") == pytest.approx(-1.25)
    with pytest.raises(ValueError):
        bundle._real("not-a-number-at-all")


def test_amplitude_record_is_sorted_by_time(tmp_path):
    """A restart can emit snapshots out of time order, and a fit needs them ordered."""
    path = str(tmp_path / "amplitude.dat")
    with open(path, "w") as fh:
        fh.write("# time  amplitude  mean_radius\n")
        for t, a in [(2.0, 0.3), (0.0, 0.1), (1.0, 0.2)]:
            fh.write(f"{t:24.16E}{a:24.16E}{1.0:24.16E}\n")

    times, amps, means = bundle.read_amplitude(path)
    assert times == pytest.approx([0.0, 1.0, 2.0])
    assert amps == pytest.approx([0.1, 0.2, 0.3])
    assert means.size == 3


def test_find_field_bundles(tmp_path):
    for name in ["field_y_velocity.dat", "field_pressure.dat",
                 "profiles.dat", "amplitude.dat"]:
        (tmp_path / name).write_text("")
    found = bundle.find_field_bundles(str(tmp_path))
    assert sorted(found) == ["pressure", "y_velocity"]


def test_rejects_a_bundle_of_the_wrong_kind(tmp_path):
    z = np.linspace(0.0, 1.0, 4)
    path = str(tmp_path / "profiles.dat")
    bundle.write_profiles(path, 1.0, z, [0.0], np.ones((1, 4)))
    with pytest.raises(ValueError):
        bundle.read_field(path)


@pytest.mark.skipif(
    not os.path.isfile(os.path.join(FORTRAN_BUNDLE, "exported.dat")),
    reason="run 'make test' in LSA_Rayleigh_Fortran to produce the bundle")
def test_reads_a_bundle_written_by_fortran():
    """The Fortran tree writes the same format this reader consumes.

    The Fortran test plants r = r0 + a cos(k z) in genuine plotfiles, recovers the profile through its own reader and exports it. Recomputing the closed form here checks the whole path across the language boundary, rather than checking that the file merely parses.
    """
    got = bundle.read_profiles(os.path.join(FORTRAN_BUNDLE, "exported.dat"))
    k = 2.0 * np.pi / got["wavelength"]
    expected = 1.0 + 0.05 * np.cos(k * got["z"])
    assert got["radius"][0] == pytest.approx(expected, abs=1e-12)


@pytest.mark.skipif(
    not os.path.isfile(os.path.join(FORTRAN_BUNDLE, "field_test_var.dat")),
    reason="run 'make test' in LSA_Rayleigh_Fortran to produce the bundle")
def test_reads_a_fortran_field_bundle_in_the_right_order():
    got = bundle.read_field(os.path.join(FORTRAN_BUNDLE,
                                         "field_test_var.dat"))
    assert got["field"] == "test_var"
    # The Fortran test stores 100*it + 10*j + i with all indices counted from one, so this entry pins the striding as well as the values.
    assert got["data"][1, 2, 4] == pytest.approx(253.0)
