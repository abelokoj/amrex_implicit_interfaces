"""Portable figure data: interface profiles and field slices as plain text.

An amplitude record is enough to draw the growth curve and the local growth rate, but two of the four figures reported for a case need more than that. The interface figure needs the radius profile r(z) at several instants, and the field figure needs a plotfile variable on the (r,z) plane together with the interface contour laid over it. Both of those were previously obtainable only by reading the plotfiles, which requires yt and requires the plotfiles themselves, and a sweep of sixteen cases leaves tens of gigabytes of plotfiles on the machine that ran it. Transferring them to a workstation to redraw a legend is not practical.

This module defines a small text format that carries exactly what the figures need and nothing else. A profile bundle for a typical case is a few hundred kilobytes and a field bundle with three panels a few hundred more, against tens of gigabytes of plotfiles, so the whole figure set travels in a directory small enough to copy over a login session. The format is written by both implementations, by `python/dump_profiles.py` here and by `bin/dump_profiles` in the Fortran tree, and the two are byte-comparable so that a bundle from either can be plotted by the same script.

The format is deliberately crude: a keyword line, then free-format numbers. Fortran writes it with ordinary formatted output and reads it with list-directed input, and Python parses it by splitting on whitespace. Nothing here depends on a library that might be absent on a cluster.

A profile bundle, `profiles.dat`:

    # LSA profile bundle v1
    WAVELENGTH   8.9760000000000000E+00
    NZ           72
    NBLOCKS      120
    Z
      <NZ numbers, free format>
    BLOCK  1   <time>
      <NZ radius values, free format>
    BLOCK  2   <time>
      ...
    END

A field bundle, `field_<name>.dat`:

    # LSA field bundle v1
    FIELD        y_velocity
    NR           32
    NZ           72
    NBLOCKS      3
    R
      <NR numbers>
    Z
      <NZ numbers>
    BLOCK  1   <time>
      <NR*NZ numbers, i varying fastest within each j>
    END

A radius of -1 in a profile block marks an axial station at which no interface was found, which happens when the interface leaves the radial extent of the domain. Those entries are returned as NaN so that a plot leaves a gap rather than drawing a spurious excursion to the axis.
"""

# ORIENTATION
#
# This module is a file format and nothing more: it writes numbers to text and reads them back. There is no physics in it.
#
# The reason it exists is practical. The solver writes plotfiles, which are directories of binary data holding the entire flow field at one instant, and a parameter sweep leaves tens of gigabytes of them on whichever machine ran it. Drawing a figure needs a tiny fraction of that: a radius profile here, one velocity field there. Rather than move the plotfiles, this module extracts the fraction the figures need into a few hundred kilobytes of plain text that can be copied anywhere.
#
# The format is a keyword line followed by free-format numbers, chosen so that both Fortran and Python can read and write it with nothing but their own standard facilities. If you are designing an interchange format yourself, this is usually the right instinct: the simplest thing that both ends can parse without a library will outlive anything cleverer.


from __future__ import annotations

import os

import numpy as np

PROFILE_MAGIC = "# LSA profile bundle v1"
FIELD_MAGIC = "# LSA field bundle v1"

#: Sentinel written for an axial station at which no interface was located. The value is negative because a radius cannot be, so it cannot collide with data.
NO_INTERFACE = -1.0


def _fmt_block(values, per_line=5):
    """Format a sequence of reals as free-format lines.

    Five to a line keeps the file readable in a pager while remaining trivial for list-directed input to consume, which ignores line breaks entirely.
    """
    out, values = [], np.asarray(values, dtype=float).ravel()
    for start in range(0, values.size, per_line):
        chunk = values[start:start + per_line]
        out.append("".join(f"{v:24.16E}" for v in chunk))
    return "\n".join(out) + "\n"


# ----------------------------------------------------------------------
# Writing
# ----------------------------------------------------------------------
def write_profiles(path, wavelength, z, times, radii):
    """Write a profile bundle.

    Parameters
    ----------
    path : str
        Destination file, conventionally `profiles.dat`.
    wavelength : float
        Axial wavelength of the deck, recorded so that the plotting script need not be told it a second time.
    z : array of shape (nz,)
        Axial cell centres, shared by every block.
    times : sequence of float
        Snapshot times, in increasing order.
    radii : array of shape (nt, nz)
        Interface radius at each station. Non-finite entries are written as the NO_INTERFACE sentinel.
    """
    z = np.asarray(z, dtype=float)
    radii = np.asarray(radii, dtype=float)
    if radii.shape != (len(times), z.size):
        raise ValueError(f"radii shape {radii.shape} does not match "
                         f"({len(times)}, {z.size})")

    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write(PROFILE_MAGIC + "\n")
        fh.write(f"WAVELENGTH {wavelength:24.16E}\n")
        fh.write(f"NZ {z.size:d}\n")
        fh.write(f"NBLOCKS {len(times):d}\n")
        fh.write("Z\n")
        fh.write(_fmt_block(z))
        for index, (t, row) in enumerate(zip(times, radii), start=1):
            clean = np.where(np.isfinite(row), row, NO_INTERFACE)
            fh.write(f"BLOCK {index:d} {float(t):24.16E}\n")
            fh.write(_fmt_block(clean))
        fh.write("END\n")
    return path


def write_field(path, field, r, z, times, data):
    """Write a field bundle.

    Parameters
    ----------
    data : array of shape (nt, nr, nz)
        The variable on the (r,z) plane at each retained instant, indexed as data[t, i, j] to match the plotfile convention of the radial index varying fastest.
    """
    r = np.asarray(r, dtype=float)
    z = np.asarray(z, dtype=float)
    data = np.asarray(data, dtype=float)
    if data.shape != (len(times), r.size, z.size):
        raise ValueError(f"data shape {data.shape} does not match "
                         f"({len(times)}, {r.size}, {z.size})")

    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write(FIELD_MAGIC + "\n")
        fh.write(f"FIELD {field}\n")
        fh.write(f"NR {r.size:d}\n")
        fh.write(f"NZ {z.size:d}\n")
        fh.write(f"NBLOCKS {len(times):d}\n")
        fh.write("R\n")
        fh.write(_fmt_block(r))
        fh.write("Z\n")
        fh.write(_fmt_block(z))
        for index, (t, plane) in enumerate(zip(times, data), start=1):
            fh.write(f"BLOCK {index:d} {float(t):24.16E}\n")
            # Fortran order, so that the radial index varies fastest and the file matches what the Fortran writer emits from the same array.
            fh.write(_fmt_block(plane.ravel(order="F")))
        fh.write("END\n")
    return path


# ----------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------
def _real(token):
    """Parse a real, tolerating the Fortran form that omits the exponent letter.

    Written with a two-digit exponent descriptor, Fortran drops the letter when the exponent needs three digits, so 4.66e-310 emerges as `4.6598606362008991-310`. The writers in this project use a three-digit descriptor and never produce that form, but a bundle written by an earlier build may carry it, and rejecting such a file outright would lose figure data that is otherwise intact. The sign is sought after the first character so that a leading minus is not mistaken for the exponent marker.
    """
    try:
        return float(token)
    except ValueError:
        pass
    for i in range(1, len(token)):
        if token[i] in "+-" and token[i - 1] not in "eEdD":
            try:
                return float(token[:i] + "E" + token[i:])
            except ValueError:
                break
    raise ValueError(f"cannot parse {token!r} as a real")


class _Tokens:
    """A whitespace token stream with the keyword lines left in place.

    The format mixes keywords and numbers, so the reader walks tokens rather than lines. Keeping the position explicit makes a malformed file fail with the offending token named instead of with an index error thirty numbers later.
    """

    def __init__(self, text):
        self.items = text.split()
        self.pos = 0

    def next(self):
        if self.pos >= len(self.items):
            raise ValueError("bundle ended before the END keyword")
        item = self.items[self.pos]
        self.pos += 1
        return item

    def expect(self, keyword):
        got = self.next()
        if got != keyword:
            raise ValueError(f"expected keyword {keyword!r}, found {got!r}")

    def reals(self, count):
        values = np.empty(count, dtype=float)
        for i in range(count):
            values[i] = _real(self.next())
        return values


def _strip_comments(text):
    """Drop comment lines so that the magic line does not enter the stream."""
    return "\n".join(line for line in text.splitlines()
                     if not line.lstrip().startswith("#"))


def read_profiles(path):
    """Read a profile bundle.

    Returns a dict with keys `wavelength`, `z`, `times` and `radius`, where `radius` has shape (nt, nz) and carries NaN wherever no interface was located.
    """
    with open(path) as fh:
        raw = fh.read()
    if PROFILE_MAGIC.split()[2] not in raw:
        raise ValueError(f"{path} is not a profile bundle")
    tok = _Tokens(_strip_comments(raw))

    tok.expect("WAVELENGTH")
    wavelength = _real(tok.next())
    tok.expect("NZ")
    nz = int(tok.next())
    tok.expect("NBLOCKS")
    nblocks = int(tok.next())
    tok.expect("Z")
    z = tok.reals(nz)

    times = np.empty(nblocks, dtype=float)
    radius = np.empty((nblocks, nz), dtype=float)
    for block in range(nblocks):
        tok.expect("BLOCK")
        tok.next()                       # the block index, for readability
        times[block] = _real(tok.next())
        radius[block] = tok.reals(nz)
    tok.expect("END")

    radius = np.where(radius <= 0.0, np.nan, radius)
    return {"wavelength": wavelength, "z": z, "times": times,
            "radius": radius}


def read_field(path):
    """Read a field bundle.

    Returns a dict with keys `field`, `r`, `z`, `times` and `data`, where `data` has shape (nt, nr, nz).
    """
    with open(path) as fh:
        raw = fh.read()
    if FIELD_MAGIC.split()[2] not in raw:
        raise ValueError(f"{path} is not a field bundle")
    tok = _Tokens(_strip_comments(raw))

    tok.expect("FIELD")
    field = tok.next()
    tok.expect("NR")
    nr = int(tok.next())
    tok.expect("NZ")
    nz = int(tok.next())
    tok.expect("NBLOCKS")
    nblocks = int(tok.next())
    tok.expect("R")
    r = tok.reals(nr)
    tok.expect("Z")
    z = tok.reals(nz)

    times = np.empty(nblocks, dtype=float)
    data = np.empty((nblocks, nr, nz), dtype=float)
    for block in range(nblocks):
        tok.expect("BLOCK")
        tok.next()
        times[block] = _real(tok.next())
        flat = tok.reals(nr * nz)
        data[block] = flat.reshape((nr, nz), order="F")
    tok.expect("END")

    return {"field": field, "r": r, "z": z, "times": times, "data": data}


def read_amplitude(path):
    """Read the three-column amplitude record written by `dump_amplitude`.

    Kept here so that a plotting script needs one import to read every part of a bundle, and so that the column order is stated in exactly one place.
    """
    table = np.loadtxt(path, comments="#", ndmin=2)
    if table.shape[1] < 2:
        raise ValueError(f"{path} has fewer than two columns")
    times = table[:, 0]
    amps = table[:, 1]
    means = table[:, 2] if table.shape[1] > 2 else np.array([])
    order = np.argsort(times)
    return times[order], amps[order], (means[order] if means.size else means)


def find_field_bundles(directory):
    """Return the field bundles in a directory, keyed by variable name."""
    found = {}
    if not os.path.isdir(directory):
        return found
    for name in sorted(os.listdir(directory)):
        if name.startswith("field_") and name.endswith(".dat"):
            found[name[len("field_"):-len(".dat")]] = os.path.join(directory,
                                                                   name)
    return found
