"""Write minimal but genuine AMReX plotfiles, for testing the pipeline.

The production pipeline reads plotfiles produced by the AMReX solver.  We cannot run the solver in a test, so instead we synthesise plotfiles in the real on-disk format and read them back through the *same* code path (nsmfp.read_plotfile -> yt).  That exercises the reader, the twin subtraction, the snapshot assembly and the DMD together, and it fails loudly if the format assumptions are wrong.

Only what the tests need is supported: a single level, a single grid, 2-D, double precision, no ghost cells.  That matches the recommended LSA configuration (amr.max_level = 0).

Format references are yt's own parser, `yt.frontends.amrex.data_structures.BoxlibHierarchy._parse_index`, and the FAB header example quoted in `_cache_endianness`.
"""

from __future__ import annotations

import os

import numpy as np

# FAB real descriptor: (bytes-per-real, (format)), (bytes, (byte order)). yt's `_cache_endianness` decides endianness from the byte-order block: ascending "1 2 ... 8" is BIG endian, descending "8 7 ... 1" is LITTLE. The example quoted in yt's source is the big-endian one, so it must NOT be copied verbatim for native x86 output.
_FAB_REAL_DESCRIPTOR_LE = "((8, (64 11 52 0 1 12 0 1023)),(8, (8 7 6 5 4 3 2 1)))"
_FAB_REAL_DESCRIPTOR_BE = "((8, (64 11 52 0 1 12 0 1023)),(8, (1 2 3 4 5 6 7 8)))"


def write_plotfile(path, fields, time, step, domain_left=(0.0, 0.0),
                   domain_right=(1.0, 1.0), coord_sys=0):
    """Write one plotfile directory.

    Parameters
    ----------
    path : str
        Directory to create (e.g. ``.../plt00100``).
    fields : mapping {name: 2-D array}
        All arrays must share a shape (nx, ny), indexed [i, j].
    time : float step : int coord_sys : int
        0 cartesian, 1 cylindrical (RZ).  The solver uses 1 for this problem;
        tests default to 0 because yt rewrites the domain edges for
        cylindrical geometry, which would complicate shape bookkeeping.
    """
    names = list(fields)
    if not names:
        raise ValueError("need at least one field")
    arrays = [np.ascontiguousarray(fields[n], dtype="<f8") for n in names]
    shape = arrays[0].shape
    if len(shape) != 2:
        raise ValueError("only 2-D fields are supported")
    for a in arrays:
        if a.shape != shape:
            raise ValueError("all fields must share a shape")

    nx, ny = shape
    ncomp = len(names)
    os.makedirs(os.path.join(path, "Level_0"), exist_ok=True)

    dx = [(domain_right[0] - domain_left[0]) / nx,
          (domain_right[1] - domain_left[1]) / ny]
    box = f"((0,0) ({nx - 1},{ny - 1}) (0,0))"

    # ---- global Header -------------------------------------------------
    h = ["HyperCLaw-V1.1", str(ncomp)]
    h += names
    h += [
        "2",                                    # dimensionality
        repr(float(time)),
        "0",                                    # finest level
        f"{domain_left[0]!r} {domain_left[1]!r}",
        f"{domain_right[0]!r} {domain_right[1]!r}",
        "",                                     # ref ratios (none at level 0)
        box,                                    # global index space
        str(int(step)),                         # steps per level
        f"{dx[0]!r} {dx[1]!r}",                 # cell size per level
        str(int(coord_sys)),
        "0",                                    # bwidth; yt asserts this is 0
        f"0 1 {float(time)!r}",                 # level, ngrids, time
        str(int(step)),                         # level steps
        f"{domain_left[0]!r} {domain_right[0]!r}",   # grid extent, x
        f"{domain_left[1]!r} {domain_right[1]!r}",   # grid extent, y
        "Level_0/Cell",
    ]
    with open(os.path.join(path, "Header"), "w") as f:
        f.write("\n".join(h) + "\n")

    # ---- FAB data file ------------------------------------------------- One FAB: an ASCII header line, then each component in Fortran order.
    fab_header = (f"FAB {_FAB_REAL_DESCRIPTOR_LE}{box} {ncomp}\n").encode("ascii")
    data_path = os.path.join(path, "Level_0", "Cell_D_00000")
    with open(data_path, "wb") as f:
        f.write(fab_header)
        for a in arrays:
            f.write(np.asfortranarray(a).astype("<f8").tobytes(order="F"))

    # ---- level header Cell_H ------------------------------------------- yt skips lines 1-2, reads ncomp, skips nghost, then expects "(N 0" followed by N box specs, a ")" line, the grid count again, and then "FabOnDisk: <file> <offset>" per grid.  The offset points at the START of the FAB, i.e. at its ASCII header line, which yt then skips.
    ch = [
        "1",                    # version
        "1",                    # how
        str(ncomp),
        "0",                    # nghost
        "(1 0",
        box,
        ")",
        "1",
        "FabOnDisk: Cell_D_00000 0",
        "",
        f"1,{ncomp}",
    ]
    ch.append(",".join(repr(float(a.min())) for a in arrays) + ",")
    ch.append("")
    ch.append(f"1,{ncomp}")
    ch.append(",".join(repr(float(a.max())) for a in arrays) + ",")
    with open(os.path.join(path, "Level_0", "Cell_H"), "w") as f:
        f.write("\n".join(ch) + "\n")

    return path


def write_series(directory, field_fn, steps, times, shape=(16, 32),
                 field_names=("x_velocity", "y_velocity"),
                 domain_right=(1.0, 1.0), prefix="nddataPLT"):
    """Write a sequence of plotfiles ``<prefix>%05d``.

    The default prefix matches what the solver actually writes: with ns.visual_nddata_format=1 it emits nddataPLT00000023, not plt00023, and it does not honour amr.plot_file.  Tests use the realistic name so they exercise the same default the production path uses.

    `field_fn(t, name, nx, ny) -> 2-D array` supplies the data.
    """
    os.makedirs(directory, exist_ok=True)
    out = []
    nx, ny = shape
    for step, t in zip(steps, times):
        p = os.path.join(directory, f"{prefix}{step:05d}")
        fields = {n: field_fn(t, n, nx, ny) for n in field_names}
        write_plotfile(p, fields, time=t, step=step,
                       domain_right=domain_right)
        out.append(p)
    return out
