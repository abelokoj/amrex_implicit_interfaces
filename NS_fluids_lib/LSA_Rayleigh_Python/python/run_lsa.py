# !/usr/bin/env python3
"""
NS-MFP + DMD stability analysis driver for the Rayleigh-Plateau problem.

Subcommands
-----------
analyse      one (perturbed, twin) pair -> growth rate, compared with theory linearity    amplitude sweep -> proportionality test verdict dispersion   wavenumber sweep -> measured vs exact dispersion curve converge     subspace-size / start-index convergence of the growth rate

Examples
--------
python run_lsa.py analyse --perturbed runs/k0.7_eps1e-3 \\
                            --base runs/k0.7_base --kr0 0.7

python run_lsa.py linearity --base runs/k0.7_base --kr0 0.7 \\
      --case 1e-4=runs/k0.7_eps1e-4 --case 1e-3=runs/k0.7_eps1e-3 \\
      --case 1e-2=runs/k0.7_eps1e-2

python run_lsa.py dispersion --root runs --eps 1e-3 \\
      --kr0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import dispersion  # noqa: E402
import dmd as dmd_mod  # noqa: E402
import nsmfp  # noqa: E402
import proportionality as prop  # noqa: E402
import report  # noqa: E402


def _theory(kr0, args):
    """Exact inviscid and estimated viscous growth rates for this run."""
    k = kr0 / args.r0
    inviscid = float(dispersion.growth_rate(
        k, r0=args.r0, sigma=args.sigma, rho=args.rho_l))
    viscous = float(dispersion.growth_rate_viscous(
        k, r0=args.r0, sigma=args.sigma, rho=args.rho_l, mu=args.mu_l))
    return inviscid, viscous


def _analyse_pair(perturbed, base, args, verbose=True):
    """Build the subspace for one twin pair and extract the growth rate.

    Two independent estimates are produced deliberately:
      * `norm_fit`  - slope of log|Q'(t)|, a direct and assumption-light
                      measurement that only sees the dominant mode;
      * `dmd`       - the leading stationary DMD eigenvalue, which also
                      yields the mode shape and the rest of the spectrum.
    They should agree once the transient has shed.  Disagreement is informative: it usually means the fit window still contains the non-modal transient, or that no single mode dominates yet.
    """
    X, dt, times = nsmfp.build_snapshot_matrix(
        perturbed, base,
        fields=tuple(args.fields), level=args.level,
        skip=args.skip, stride=args.stride,
        max_snapshots=args.max_snapshots, verbose=verbose)

    norms = np.linalg.norm(X, axis=0)
    if np.allclose(norms, 0.0):
        raise ValueError(
            "the perturbation field is identically zero: the two runs "
            "produced identical output.  Check that the perturbed deck "
            "really has ns.radblob != 0.")

    t_min = args.t_min if args.t_min is not None else times[len(times) // 3]
    rate_fit, _, r2 = nsmfp.fit_exponential_growth(
        times, norms, t_min=t_min, t_max=args.t_max)

    res = dmd_mod.dmd(X, dt, rank=args.rank, tol=args.tol)
    rate_dmd = res.leading_growth_rate()

    return dict(times=times, norms=norms, X=X, dt=dt,
                rate_fit=rate_fit, r2=r2, rate_dmd=rate_dmd, result=res,
                t_min=t_min)


def _print_comparison(kr0, rate_dmd, rate_fit, r2, inviscid, viscous):
    print(report.table(
        [["DMD: leading stationary mode", f"{rate_dmd:.6f}", ""],
         ["log-norm slope fit", f"{rate_fit:.6f}", f"r^2 = {r2:.5f}"],
         ["exact inviscid (Rayleigh)", f"{inviscid:.6f}", "theory"],
         ["viscous estimate", f"{viscous:.6f}", "theory"]],
        ["quantity", "growth rate", "note"],
        colalign=["left", "right", "left"],
        title=f"Rayleigh-Plateau growth rate at k r0 = {kr0:g}"))

    rows = []
    if inviscid > 0:
        rows.append(["measured vs inviscid",
                     report.pct((rate_dmd - inviscid) / inviscid), ""])
    if viscous > 0:
        rows.append(["measured vs viscous",
                     report.pct((rate_dmd - viscous) / viscous), ""])
    if np.isfinite(rate_dmd) and np.isfinite(rate_fit):
        spread = abs(rate_dmd - rate_fit) / max(abs(rate_dmd), 1e-30)
        note = "consistent" if spread < 0.05 else \
            "DISAGREE: check the fit window"
        rows.append(["DMD vs norm fit", report.pct(spread), note])
    if rows:
        print(report.table(rows, ["comparison", "departure", "verdict"],
                           colalign=["left", "right", "left"]))


def cmd_analyse(args):
    inviscid, viscous = _theory(args.kr0, args)
    out = _analyse_pair(args.perturbed, args.base, args)

    print()
    print(out["result"].summary(n=args.n_modes))
    _print_comparison(args.kr0, out["rate_dmd"], out["rate_fit"], out["r2"],
                      inviscid, viscous)

    if args.save:
        payload = dict(
            kr0=args.kr0, dt=out["dt"],
            growth_rate_dmd=out["rate_dmd"],
            growth_rate_fit=out["rate_fit"], r_squared=out["r2"],
            growth_rate_exact_inviscid=inviscid,
            growth_rate_viscous_estimate=viscous,
            n_snapshots=int(out["X"].shape[1]),
            n_dof=int(out["X"].shape[0]),
            rank=int(out["result"].rank),
            times=out["times"].tolist(),
            norms=out["norms"].tolist(),
            ritz_real=out["result"].ritz.real.tolist(),
            ritz_imag=out["result"].ritz.imag.tolist(),
            growth_rates=out["result"].growth_rate.tolist(),
            frequencies=out["result"].frequency.tolist(),
            amplitudes=out["result"].amplitude.tolist(),
        )
        with open(args.save, "w") as f:
            json.dump(payload, f, indent=2)
        print(f"\nwrote {args.save}")

    if args.plot:
        _plot_single(out, args.kr0, inviscid, args.plot)


def cmd_linearity(args):
    cases = {}
    for spec in args.case:
        if "=" not in spec:
            raise SystemExit(f"--case expects EPS=DIR, got {spec!r}")
        eps_s, d = spec.split("=", 1)
        cases[float(eps_s)] = d

    runs = prop.collect_amplitude_runs(
        cases, args.base, fields=tuple(args.fields), level=args.level,
        t_min=args.t_min, t_max=args.t_max, verbose=True)
    report = prop.assess_linearity(runs, scaled_tol=args.scaled_tol,
                                   rate_tol=args.rate_tol)
    print()
    print(report)

    inviscid, viscous = _theory(args.kr0, args)
    print(f"\nexact inviscid growth rate at k r0={args.kr0:g}: {inviscid:.6f}")
    print(f"viscous estimate                        : {viscous:.6f}")
    print(f"recommended amplitude for production    : "
          f"{prop.recommend_amplitude(report):.3e}")

    if args.plot:
        _plot_linearity(runs, args.plot)
    if not report.is_linear:
        # Non-zero exit so a driving script notices the sweep is unusable.
        raise SystemExit(2)


def cmd_dispersion(args):
    rows = []
    for kr0 in args.kr0:
        pert = args.perturbed_pattern.format(kr0=kr0, eps=args.eps)
        base = args.base_pattern.format(kr0=kr0)
        pert = os.path.join(args.root, pert)
        base = os.path.join(args.root, base)
        if not os.path.isdir(pert) or not os.path.isdir(base):
            print(f"[k r0={kr0:g}] SKIP (missing {pert} or {base})")
            continue
        print(f"[k r0={kr0:g}] reading ...")
        try:
            out = _analyse_pair(pert, base, args, verbose=False)
        except Exception as exc:                      # keep the sweep going
            print(f"[k r0={kr0:g}] FAILED: {exc}")
            continue
        inviscid, viscous = _theory(kr0, args)
        rows.append(dict(kr0=kr0, measured=out["rate_dmd"],
                         measured_fit=out["rate_fit"],
                         inviscid=inviscid, viscous=viscous))
        print(f"[k r0={kr0:g}] DMD {out['rate_dmd']:.6f}  "
              f"exact {inviscid:.6f}")

    if not rows:
        raise SystemExit("no runs could be analysed")

    print()
    print(f"{'k r0':>7} {'DMD':>12} {'norm fit':>12} "
          f"{'inviscid':>12} {'viscous':>12} {'err vs visc':>12}")
    for r in rows:
        err = (100 * (r["measured"] - r["viscous"]) / r["viscous"]
               if r["viscous"] > 0 else float("nan"))
        print(f"{r['kr0']:7.3f} {r['measured']:12.6f} {r['measured_fit']:12.6f} "
              f"{r['inviscid']:12.6f} {r['viscous']:12.6f} {err:11.2f}%")

    peak = max(rows, key=lambda r: r["measured"])
    k_exact, s_exact = dispersion.most_unstable_wavenumber(
        r0=args.r0, sigma=args.sigma, rho=args.rho_l)
    print(f"\nmost unstable sampled k r0 : {peak['kr0']:.3f} "
          f"(sigma = {peak['measured']:.6f})")
    print(f"theoretical peak           : k r0 = {k_exact * args.r0:.3f} "
          f"(sigma = {s_exact:.6f})")

    if args.save:
        with open(args.save, "w") as f:
            json.dump(rows, f, indent=2)
        print(f"\nwrote {args.save}")
    if args.plot:
        _plot_dispersion(rows, args, args.plot)


def cmd_converge(args):
    X, dt, times = nsmfp.build_snapshot_matrix(
        args.perturbed, args.base, fields=tuple(args.fields),
        level=args.level, skip=args.skip, stride=args.stride, verbose=True)
    runs = dmd_mod.convergence_study(X, dt, sizes=args.sizes,
                                     starts=tuple(args.starts),
                                     rank=args.rank, tol=args.tol)
    inviscid, viscous = _theory(args.kr0, args)

    print()
    print(f"{'start':>7} {'n_snap':>8} {'rank':>6} {'growth':>13} "
          f"{'err vs visc':>12}")
    for r in runs:
        err = (100 * (r["growth_rate"] - viscous) / viscous
               if viscous > 0 else float("nan"))
        print(f"{r['start']:7d} {r['n_snapshots']:8d} {r['rank']:6d} "
              f"{r['growth_rate']:13.6f} {err:11.2f}%")

    rates = np.array([r["growth_rate"] for r in runs], dtype=float)
    rates = rates[np.isfinite(rates)]
    if rates.size > 1:
        spread = rates.max() - rates.min()
        print(f"\nspread over configurations: {spread:.3e} "
              f"({spread / max(abs(rates.mean()), 1e-30):.2%} of the mean)")
        print("Converged if this spread is small compared with the "
              "difference from theory.")
    print(f"exact inviscid: {inviscid:.6f}   viscous estimate: {viscous:.6f}")


# ----------------------------------------------------------------- plots
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _plot_single(out, kr0, inviscid, path):
    plt = _mpl()
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))

    ax[0].semilogy(out["times"], out["norms"], lw=1.6, label=r"$\|Q'\|_2$")
    m = out["times"] >= out["t_min"]
    ax[0].semilogy(out["times"][m],
                   out["norms"][m][0] * np.exp(
                       out["rate_fit"] * (out["times"][m] - out["times"][m][0])),
                   "--", lw=1.4, label=f"fit $\\sigma$={out['rate_fit']:.4f}")
    ax[0].semilogy(out["times"][m],
                   out["norms"][m][0] * np.exp(
                       inviscid * (out["times"][m] - out["times"][m][0])),
                   ":", lw=1.4, label=f"exact $\\sigma$={inviscid:.4f}")
    ax[0].set_xlabel("t")
    ax[0].set_ylabel(r"$\|Q'\|_2$")
    ax[0].set_title(f"Perturbation growth, $kr_0$={kr0:g}")
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3)

    r = out["result"]
    sc = ax[1].scatter(r.growth_rate, r.frequency, s=22,
                       c=np.log10(np.maximum(r.amplitude, 1e-300)),
                       cmap="viridis")
    ax[1].axvline(0.0, color="k", lw=0.8, ls="--")
    ax[1].axvline(inviscid, color="crimson", lw=1.0, ls=":",
                  label="exact inviscid")
    ax[1].set_xlabel(r"growth rate $\sigma$")
    ax[1].set_ylabel(r"frequency $\omega$")
    ax[1].set_title("DMD spectrum")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)
    fig.colorbar(sc, ax=ax[1], label=r"$\log_{10}$ amplitude")

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")


def _plot_linearity(runs, path):
    plt = _mpl()
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for r in runs:
        ax[0].semilogy(r.times, r.norms, lw=1.4, label=f"$\\epsilon$={r.eps:g}")
        ax[1].semilogy(r.times, r.norms / r.eps, lw=1.4,
                       label=f"$\\epsilon$={r.eps:g}")
    ax[0].set_xlabel("t"); ax[0].set_ylabel(r"$\|Q'\|_2$")
    ax[0].set_title("Raw perturbation norm")
    ax[1].set_xlabel("t"); ax[1].set_ylabel(r"$\|Q'\|_2/\epsilon$")
    ax[1].set_title("Scaled: curves must collapse if linear")
    for a in ax:
        a.legend(fontsize=8); a.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(path, dpi=150)
    print(f"wrote {path}")


def _plot_dispersion(rows, args, path):
    plt = _mpl()
    kk = np.linspace(1e-3, 1.0, 400) / args.r0
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(kk * args.r0,
            dispersion.growth_rate(kk, r0=args.r0, sigma=args.sigma,
                                   rho=args.rho_l),
            "-", lw=1.6, label="exact (inviscid)")
    ax.plot(kk * args.r0,
            dispersion.growth_rate_viscous(kk, r0=args.r0, sigma=args.sigma,
                                           rho=args.rho_l, mu=args.mu_l),
            "--", lw=1.4, label="viscous estimate")
    ax.plot([r["kr0"] for r in rows], [r["measured"] for r in rows],
            "o", ms=7, label="NS-MFP + DMD")
    ax.set_xlabel(r"$k r_0$")
    ax.set_ylabel(r"growth rate $\sigma$")
    ax.set_title("Rayleigh-Plateau dispersion relation")
    ax.legend(fontsize=9); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(path, dpi=150)
    print(f"wrote {path}")


# ------------------------------------------------------------------ CLI
def _add_common(p):
    p.add_argument("--fields", nargs="+", default=list(nsmfp.DEFAULT_FIELDS),
                   help="plotfile variables forming the state vector")
    p.add_argument("--level", type=int, default=0)
    p.add_argument("--skip", type=int, default=0,
                   help="discard leading snapshots (shed the transient)")
    p.add_argument("--stride", type=int, default=1)
    p.add_argument("--max-snapshots", type=int, dest="max_snapshots")
    p.add_argument("--rank", type=int, default=None, help="DMD truncation rank")
    p.add_argument("--tol", type=float, default=1e-10,
                   help="singular-value cutoff when --rank is unset")
    p.add_argument("--t-min", type=float, dest="t_min",
                   help="start of the growth-rate fit window")
    p.add_argument("--t-max", type=float, dest="t_max")
    p.add_argument("--r0", type=float, default=1.0)
    p.add_argument("--sigma", type=float, default=1.0)
    p.add_argument("--rho-l", type=float, dest="rho_l", default=1.0)
    p.add_argument("--mu-l", type=float, dest="mu_l", default=0.02)
    p.add_argument("--plot", help="write a figure to this path")
    p.add_argument("--save", help="write results as JSON to this path")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("analyse", help="one twin pair -> growth rate")
    p.add_argument("--perturbed", required=True)
    p.add_argument("--base", required=True)
    p.add_argument("--kr0", type=float, required=True)
    p.add_argument("--n-modes", type=int, default=10, dest="n_modes")
    _add_common(p)
    p.set_defaults(func=cmd_analyse)

    p = sub.add_parser("linearity", help="amplitude sweep -> proportionality")
    p.add_argument("--base", required=True)
    p.add_argument("--kr0", type=float, required=True)
    p.add_argument("--case", action="append", required=True,
                   metavar="EPS=DIR", help="repeatable, e.g. --case 1e-3=runs/a")
    p.add_argument("--scaled-tol", type=float, default=0.02, dest="scaled_tol")
    p.add_argument("--rate-tol", type=float, default=0.02, dest="rate_tol")
    _add_common(p)
    p.set_defaults(func=cmd_linearity)

    p = sub.add_parser("dispersion", help="wavenumber sweep -> curve")
    p.add_argument("--root", default="runs")
    p.add_argument("--kr0", type=float, nargs="+", required=True)
    p.add_argument("--eps", type=float, default=1e-3)
    p.add_argument("--perturbed-pattern", default="k{kr0:g}_eps{eps:g}",
                   dest="perturbed_pattern")
    p.add_argument("--base-pattern", default="k{kr0:g}_base",
                   dest="base_pattern")
    _add_common(p)
    p.set_defaults(func=cmd_dispersion)

    p = sub.add_parser("converge", help="subspace convergence study")
    p.add_argument("--perturbed", required=True)
    p.add_argument("--base", required=True)
    p.add_argument("--kr0", type=float, required=True)
    p.add_argument("--sizes", type=int, nargs="+", default=None)
    p.add_argument("--starts", type=int, nargs="+", default=[0])
    _add_common(p)
    p.set_defaults(func=cmd_converge)

    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
