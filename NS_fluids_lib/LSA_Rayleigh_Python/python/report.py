"""Formatted table output for the reports.

Every report this package produces is a small table, and the reader of a thesis appendix or a batch log should be able to take in a growth rate, its residual and its departure from theory without counting columns. This module renders aligned tables, using `tabulate` when it is installed and falling back to an equivalent renderer when it is not.

The fallback exists because the analysis must run on machines where installing a package is inconvenient, and because a missing formatting dependency is not a reason for a calculation that took hours to fail at the point of printing its answer. Output is identical in structure either way, so a log file is comparable across machines regardless of which path produced it; the `simple` table format of `tabulate` was chosen for that reason, since it is what the fallback reproduces.

Plain ASCII is used rather than box-drawing characters. Output is read through terminals, log files, `less`, and text editors whose locale settings are not known in advance, and a table that degrades into replacement characters is worse than one that was never drawn.

The Fortran implementation writes the same shapes through `src/table_mod.f90`, so that the two implementations may be compared line by line.
"""

from __future__ import annotations

try:
    from tabulate import tabulate as _tabulate
    HAVE_TABULATE = True
except ImportError:
    _tabulate = None
    HAVE_TABULATE = False


def _fallback(rows, headers, colalign):
    """Render a table without `tabulate`, matching its `simple` format."""
    cells = [[str(c) for c in row] for row in rows]
    head = [str(h) for h in headers]
    ncol = len(head)
    widths = [len(head[j]) for j in range(ncol)]
    for row in cells:
        for j in range(min(ncol, len(row))):
            widths[j] = max(widths[j], len(row[j]))

    def line(vals):
        out = []
        for j in range(ncol):
            v = vals[j] if j < len(vals) else ""
            how = colalign[j] if colalign and j < len(colalign) else "right"
            out.append(v.ljust(widths[j]) if how == "left"
                       else v.rjust(widths[j]))
        return "  ".join(out).rstrip()

    total = sum(widths) + 2 * (ncol - 1)
    return "\n".join([line(head), "-" * total]
                     + [line(r) for r in cells] + ["-" * total])


def table(rows, headers, colalign=None, title=None):
    """Render a table as a string.

    Cells are supplied already formatted, which keeps the numeric formatting decisions with the code that knows what each number means rather than in the renderer. `colalign` defaults to left for the first column and right for the rest, which is the layout every report here wants.
    """
    ncol = len(headers)
    if colalign is None:
        colalign = ["left"] + ["right"] * (ncol - 1)
    cells = [[str(c) for c in row] for row in rows]

    if HAVE_TABULATE:
        body = _tabulate(cells, headers=list(headers), tablefmt="simple",
                         colalign=tuple(colalign), disable_numparse=True)
    else:
        body = _fallback(cells, headers, colalign)

    return f"\n{title}\n{body}" if title else body


def kv(pairs, title=None):
    """Render an aligned key-value block, for summarising a single case.

    Labels are padded to a common width so that the values form a column. This is the same idea as `table` with the header suppressed, and is used where there is one value per line rather than several cases to compare.
    """
    items = [(str(k), str(v)) for k, v in pairs]
    if not items:
        return ""
    wl = max(len(k) for k, _ in items)
    wv = max(len(v) for _, v in items)
    lines = [f"{k.ljust(wl)}  {v.rjust(wv)}" for k, v in items]
    if title:
        return "\n" + title + "\n" + "-" * (wl + 2 + wv) + "\n" + "\n".join(lines)
    return "\n".join(lines)


def pct(frac, places=2):
    """A percentage, given the fraction, with a trailing unit so that a column of them reads as percentages without one in the header."""
    return f"{100.0 * frac:.{places}f} %"
