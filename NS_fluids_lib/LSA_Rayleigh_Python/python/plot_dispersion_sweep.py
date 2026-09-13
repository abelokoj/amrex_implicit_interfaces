"""Measured dispersion relation, drawn from the amplitude records of a sweep.

The figure distributed with the archive, `fig4_dispersion`, shows a single measured point at the reference wavenumber against the analytic curves, because a single point was all that existed when it was written. This draws the whole measured curve from whatever records a sweep has produced, so that the comparison against theory is made across the band rather than at one wavenumber.

Two panels are produced. The upper shows the growth rate against the inviscid and viscous predictions. The lower shows the departure from the viscous prediction as a percentage, which is where the character of the error becomes legible: a discretisation error appears as a roughly constant offset across the band, whereas a mistake in the deck or in the fit window appears as scatter or as a trend that changes sign.

Wavenumbers are recovered from the record labels, which the sweep writes as `disp<NNN>` for the dispersion cases, with the reference case `res08` supplying the point at 0.7. Cases outside the dispersion sweep, such as the amplitude and resolution studies, are ignored.

Usage, from the working directory:

    python3 plotting/plot_dispersion_sweep.py
    python3 plotting/plot_dispersion_sweep.py --frac 0.6 --outdir figures
"""

import argparse
import glob
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import dispersion
import transient  # noqa: E402
import pub_style  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

#: Width of a reference or guide line, matched to the axis frame so that a guide never competes with a data curve. Data curves take lines.linewidth from pub_style, which is 1.0, and therefore carry no explicit width anywhere in this file.
GUIDE_LW = 0.8



def kr0_of(label):
    """Recover the wavenumber a case was run at from its record label."""
    if label.startswith("disp"):
        digits = "".join(c for c in label[4:] if c.isdigit())
        if digits:
            return int(digits) / 100.0
    if label == "res08":
        return 0.7
    return None


def fit(t, a, t_min):
    """Log-linear fit of a(t) over t >= t_min, returning the rate and r squared."""
    m = (t >= t_min) & (a > 0.0) & np.isfinite(a)
    if m.sum() < 3:
        return float("nan"), float("nan")
    x, y = t[m], np.log(a[m])
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (slope * x + intercept)
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - float(np.sum(resid ** 2)) / ss_tot if ss_tot > 0 else float("nan")
    return float(slope), r2


def collect(records, frac, t_min_explicit):
    """Fit every dispersion record, returning wavenumbers, rates and residuals."""
    out = []
    for path in records:
        label = os.path.basename(path)[:-4]
        kr0 = kr0_of(label)
        if kr0 is None:
            continue
        try:
            d = np.loadtxt(path)
        except Exception:
            continue
        if d.ndim != 2 or d.shape[0] < 3:
            continue
        t, a = d[:, 0], d[:, 1]
        t_min = t_min_explicit if t_min_explicit is not None else \
            t[0] + frac * (t[-1] - t[0])
        # Above kr0 = 1 the mode oscillates rather than growing, so a fitted rate there is not a growth rate and must not be drawn on this curve. Such cases are reported separately by plotting/fit_records.py.
        if dispersion.growth_rate_squared(kr0, 1.0, 1.0, 1.0) < 0.0:
            continue
        sigma, r2 = fit(t, a, t_min)
        if math.isfinite(sigma):
            out.append((kr0, sigma, r2, label))
    out.sort()
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="measured dispersion relation")
    ap.add_argument("records", nargs="*",
                    help="record files; defaults to runs/records/*.dat")
    ap.add_argument("--frac", type=float, default=0.5,
                    help="fraction of each record discarded as transient")
    ap.add_argument("--t-min", type=float, default=None,
                    help="explicit start of the fit window, overriding --frac")
    ap.add_argument("--mu", type=float, default=2.0e-2,
                    help="dynamic viscosity of the liquid, for the viscous curve")
    ap.add_argument("--outdir", default="figures")
    args = ap.parse_args(argv)

    files = args.records or sorted(glob.glob("runs/records/*.dat"))
    pts = collect(files, args.frac, args.t_min)
    if not pts:
        print("no dispersion records found; expected runs/records/disp*.dat",
              file=sys.stderr)
        return 1

    k = np.array([p[0] for p in pts])
    s = np.array([p[1] for p in pts])

    x = np.linspace(1e-3, 1.2, 600)
    inviscid = np.array([float(dispersion.growth_rate(v, 1.0, 1.0, 1.0))
                         for v in x])
    viscous = np.array([float(dispersion.growth_rate_viscous(
        v, 1.0, 1.0, 1.0, args.mu)) for v in x])
    vis_at = np.array([float(dispersion.growth_rate_viscous(
        v, 1.0, 1.0, 1.0, args.mu)) for v in k])

    pub_style.apply_style()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.0, 6.4), sharex=True,
                                   gridspec_kw={"height_ratios": [2.1, 1.0]})

    ax1.plot(x, inviscid, "-", label="inviscid, Rayleigh (1878)")
    ax1.plot(x, viscous, "--",
             label=rf"viscous estimate, $\mathrm{{Oh}} = {args.mu:g}$")
    ax1.plot(k, s, "o", ms=6, mfc="none", mew=1.6, color="k",
             label="measured, NS-MFP")
    k_peak, _ = dispersion.most_unstable_wavenumber(1.0, 1.0, 1.0)
    ax1.axvline(k_peak, ls=":", lw=GUIDE_LW, color="0.5")
    ax1.axhline(0.0, lw=GUIDE_LW, color="0.7")
    # The band above kr0 = 1 is stable: the mode oscillates and does not grow. Marking it keeps the reader from expecting measured points there.
    x_hi = max(1.2, float(k.max()) * 1.1)
    ax1.axvspan(1.0, x_hi, color="0.92", zorder=0)
    ax1.annotate("stable, oscillatory", xy=(1.02, 0.30), fontsize=8,
                 rotation=90, va="center", color="0.35")
    ax1.set_ylabel(r"$\sigma$")
    ax1.legend(loc="upper left", fontsize=8)

    # Departure from the viscous prediction, which is the appropriate comparison at finite Ohnesorge number. A flat offset indicates discretisation error; scatter indicates something else.
    good = vis_at > 1e-12
    dev = np.full_like(s, np.nan)
    dev[good] = (s[good] - vis_at[good]) / vis_at[good] * 100.0
    ax2.plot(k[good], dev[good], "s", ms=5, mfc="none", mew=1.4, color="k")
    ax2.axhline(0.0, lw=GUIDE_LW, color="0.4")
    if np.isfinite(dev).sum() >= 2:
        mean_dev = float(np.nanmean(dev[k <= 0.7]))
        ax2.axhline(mean_dev, ls="--", lw=GUIDE_LW, color="0.6")
        ax2.annotate(rf"mean ${mean_dev:.2f}\,\%$ for $kr_0 \leq 0.7$",
                     xy=(0.04, 0.08), xycoords="axes fraction", fontsize=8)
    ax2.set_xlabel("$k r_0$")
    ax2.set_ylabel(r"departure from viscous, \%")
    ax2.set_xlim(0.0, max(1.2, float(k.max()) * 1.1))

    for ax in (ax1, ax2):
        pub_style.cm_both(ax)

    fig.tight_layout()
    os.makedirs(args.outdir, exist_ok=True)
    pub_style.savefig_all(fig, "fig5_dispersion_measured", outdir=args.outdir)
    plt.close(fig)

    print(f"{'k r0':>6}{'measured':>11}{'viscous':>11}{'departure':>12}{'r^2':>9}")
    for (kk, ss, rr, _), vv in zip(pts, vis_at):
        d = (ss - vv) / vv * 100.0 if vv > 1e-12 else float("nan")
        print(f"{kk:>6.2f}{ss:>11.6f}{vv:>11.6f}{d:>11.2f}%{rr:>9.5f}")
    print(f"\nwrote {args.outdir}/fig5_dispersion_measured.[png|pdf|svg]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
