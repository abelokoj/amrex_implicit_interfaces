# !/usr/bin/env python3
"""
Generate matched AMReX input decks for the Rayleigh-Plateau NS-MFP study.

Every deck is expanded from inputs/inputs_rayleigh_template.txt so that a perturbed run and its unperturbed twin differ in exactly one token (`ns.radblob`).  That is the whole point: the twin subtraction stands in for the NS-MFP body force B_f, and it is only valid if nothing else differs. Hand-editing decks is how that invariant gets broken silently, so generate them instead.

Typical use
-----------
    # one wavenumber, one amplitude, plus its twin
    python make_inputs.py single --kr0 0.7 --eps 1e-3 --outdir runs

    # amplitude sweep for the proportionality test (shares one twin)
    python make_inputs.py sweep-eps --kr0 0.7 --eps 1e-4 1e-3 1e-2 \
        --outdir runs

    # wavenumber sweep to trace the whole dispersion curve
    python make_inputs.py sweep-k --kr0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 \
        --eps 1e-3 --outdir runs
"""

from __future__ import annotations

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import dispersion  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATE = os.path.join(HERE, "..", "inputs", "inputs_rayleigh_template.txt")

# Default physical parameters.  Chosen to match inputs.growthrate.LSA so the results are comparable with the deck that ships with the solver.
DEFAULTS = dict(
    r0=1.0,
    sigma=1.0,
    rho_l=1.0,
    rho_g=0.001225,
    mu_l=0.02,
    mu_g=0.00026,
    rmax=4.0,        # domain radius; 4 r0 keeps the outer wall far from the
                     # interface without wasting cells on quiescent gas
    cells_per_r0=32, # radial resolution
    safety=0.25,     # fraction of the capillary-wave dt limit to use
    n_periods=5.0,   # e-folding times to march; 5.0 clears the modal
                     # transient, 2.6 does not and biases the rate low
    snapshots=300,   # target number of plotfiles
)


def capillary_dt(dx, sigma, rho_l, rho_g):
    """Brackbill capillary-wave time-step limit.

        dt < sqrt( (rho_l + rho_g) dx^3 / (2 pi sigma) )

    Surface tension is treated implicitly in this solver (ns.implicit_surface_tension), which relaxes the constraint, but we stay under the explicit limit because the perturbation we are measuring is tiny and we do not want temporal error competing with it.
    """
    return math.sqrt((rho_l + rho_g) * dx**3 / (2.0 * math.pi * sigma))


def plan_run(kr0, eps, role, **kw):
    """Work out every derived quantity for one deck."""
    p = dict(DEFAULTS)
    p.update({k: v for k, v in kw.items() if v is not None})

    r0 = p["r0"]
    k = kr0 / r0
    lam = 2.0 * math.pi / k

    dr = r0 / p["cells_per_r0"]
    nr = _round_to_block(p["rmax"] / dr)
    dr = p["rmax"] / nr
    # Aim for near-square cells, then round nz to a block-friendly value.
    nz = _round_to_block(lam / dr)
    dz = lam / nz

    dt_cap = capillary_dt(min(dr, dz), p["sigma"], p["rho_l"], p["rho_g"])
    fixed_dt = p["safety"] * dt_cap

    sigma_exact = float(dispersion.growth_rate(
        k, r0=r0, sigma=p["sigma"], rho=p["rho_l"]))
    mu_eff = p["mu_l"]
    sigma_visc = float(dispersion.growth_rate_viscous(
        k, r0=r0, sigma=p["sigma"], rho=p["rho_l"], mu=mu_eff))

    # March for n_periods e-folding times of the *viscous* estimate (the slower, more conservative one).  In the stable band there is no e-folding time, so fall back on the capillary time scale.
    t_c = math.sqrt(p["rho_l"] * r0**3 / p["sigma"])
    if sigma_visc > 1e-12:
        stop_time = p["n_periods"] / sigma_visc
    else:
        stop_time = p["n_periods"] * t_c

    max_step = int(math.ceil(stop_time / fixed_dt))
    plot_int = max(1, max_step // p["snapshots"])
    # Land max_step on an exact multiple of plot_int so both twins produce the same final plotfile.
    max_step = plot_int * (max_step // plot_int)
    check_int = max(plot_int * 50, plot_int)

    return dict(
        role=role, kr0=kr0, k=k, eps=eps, lam=lam,
        nr=nr, nz=nz, dr=dr, dz=dz,
        fixed_dt=fixed_dt, dt_cap=dt_cap,
        stop_time=max_step * fixed_dt, max_step=max_step,
        plot_int=plot_int, check_int=check_int,
        n_snapshots=max_step // plot_int,
        sigma_exact=sigma_exact, sigma_visc=sigma_visc,
        oh=dispersion.ohnesorge(mu_eff, p["rho_l"], p["sigma"], r0),
        **p,
    )


def _round_to_block(x, block=8):
    """Round a cell count up to a multiple of `block` (>= block).

    amr.blocking_factor is 2 here, but keeping cell counts divisible by 8 leaves headroom for coarsening in the multigrid solve.
    """
    n = int(math.ceil(x))
    return max(block, ((n + block - 1) // block) * block)


def render(plan, template=TEMPLATE):
    with open(template) as f:
        text = f.read()

    subs = {
        "@ROLE@": plan["role"],
        "@KR0@": f"{plan['kr0']:.6g}",
        "@LAMBDA@": f"{plan['lam']:.10g}",
        "@EPS@": f"{plan['eps']:.10g}",
        "@SIGMA_EXACT@": f"{plan['sigma_exact']:.10g}",
        "@RMAX@": f"{plan['rmax']:.10g}",
        "@NR@": str(plan["nr"]),
        "@NZ@": str(plan["nz"]),
        "@DR@": f"{plan['dr']:.6g}",
        "@DT_CAP@": f"{plan['dt_cap']:.6g}",
        "@FIXED_DT@": f"{plan['fixed_dt']:.10g}",
        "@MAX_STEP@": str(plan["max_step"]),
        "@STOP_TIME@": f"{plan['stop_time']:.10g}",
        "@PLOT_INT@": str(plan["plot_int"]),
        "@CHECK_INT@": str(plan["check_int"]),
        "@SIGMA@": f"{plan['sigma']:.10g}",
        "@RHO_L@": f"{plan['rho_l']:.10g}",
        "@RHO_G@": f"{plan['rho_g']:.10g}",
        "@MU_L@": f"{plan['mu_l']:.10g}",
        "@MU_G@": f"{plan['mu_g']:.10g}",
    }
    for key, val in subs.items():
        text = text.replace(key, val)

    leftover = [t for t in ("@" + s.split("@")[1] + "@"
                            for s in text.split("@")[1::2]) if t in text]
    if leftover:
        raise RuntimeError(f"unsubstituted template tokens: {sorted(set(leftover))}")
    return text


def write_deck(plan, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, name)
    with open(path, "w") as f:
        f.write(render(plan))
    return path


def describe(plan):
    return (
        f"  k r0={plan['kr0']:<5g} lambda={plan['lam']:<9.4f} eps={plan['eps']:<9g}\n"
        f"    grid {plan['nr']}x{plan['nz']}  dr={plan['dr']:.4g} dz={plan['dz']:.4g}\n"
        f"    dt={plan['fixed_dt']:.4g} (cap limit {plan['dt_cap']:.4g})  "
        f"steps={plan['max_step']}  snapshots={plan['n_snapshots']}\n"
        f"    sigma_exact(inviscid)={plan['sigma_exact']:.6f}  "
        f"sigma(viscous est)={plan['sigma_visc']:.6f}  Oh={plan['oh']:.4g}"
    )


def _common_physics_args(ap):
    ap.add_argument("--rmax", type=float, help="domain radius (default 4 r0)")
    ap.add_argument("--cells-per-r0", type=int, dest="cells_per_r0",
                    help="radial cells per column radius (default 32)")
    ap.add_argument("--sigma", type=float, help="surface tension")
    ap.add_argument("--rho-l", type=float, dest="rho_l")
    ap.add_argument("--rho-g", type=float, dest="rho_g")
    ap.add_argument("--mu-l", type=float, dest="mu_l")
    ap.add_argument("--mu-g", type=float, dest="mu_g")
    ap.add_argument("--n-periods", type=float, dest="n_periods",
                    help="e-folding times to march (default 6)")
    ap.add_argument("--snapshots", type=int, help="target plotfile count")
    ap.add_argument("--safety", type=float,
                    help="fraction of capillary dt limit (default 0.25)")
    ap.add_argument("--outdir", default="runs", help="output directory")


def _phys(args):
    return {k: getattr(args, k, None) for k in
            ("rmax", "cells_per_r0", "sigma", "rho_l", "rho_g", "mu_l",
             "mu_g", "n_periods", "snapshots", "safety")}


def cmd_single(args):
    made = []
    twin = plan_run(args.kr0, 0.0, "unperturbed twin (base state)", **_phys(args))
    made.append((twin, write_deck(twin, args.outdir,
                                  f"inputs.k{args.kr0:g}.base")))
    pert = plan_run(args.kr0, args.eps, "perturbed", **_phys(args))
    made.append((pert, write_deck(pert, args.outdir,
                                  f"inputs.k{args.kr0:g}.eps{args.eps:g}")))
    _report(made)


def cmd_sweep_eps(args):
    made = []
    twin = plan_run(args.kr0, 0.0, "unperturbed twin (base state)", **_phys(args))
    made.append((twin, write_deck(twin, args.outdir,
                                  f"inputs.k{args.kr0:g}.base")))
    for eps in args.eps:
        p = plan_run(args.kr0, eps, f"perturbed eps={eps:g}", **_phys(args))
        made.append((p, write_deck(p, args.outdir,
                                   f"inputs.k{args.kr0:g}.eps{eps:g}")))
    _report(made)
    print("\nAll amplitudes share the single base deck above: the twin is\n"
          "independent of eps, so it is run once and reused.")


def cmd_sweep_k(args):
    made = []
    eps = args.eps[0] if isinstance(args.eps, list) else args.eps
    for kr0 in args.kr0:
        t = plan_run(kr0, 0.0, "unperturbed twin (base state)", **_phys(args))
        made.append((t, write_deck(t, args.outdir, f"inputs.k{kr0:g}.base")))
        p = plan_run(kr0, eps, "perturbed", **_phys(args))
        made.append((p, write_deck(p, args.outdir,
                                   f"inputs.k{kr0:g}.eps{eps:g}")))
    _report(made)
    print("\nEach wavenumber needs its own twin: the box length changes with\n"
          "k, so the base state is a different discrete problem each time.")


def _report(made):
    print(f"wrote {len(made)} deck(s)\n")
    for plan, path in made:
        print(os.path.basename(path))
        print(describe(plan))
        print()
    total = sum(p["max_step"] for p, _ in made)
    print(f"total time steps across all decks: {total:,}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("single", help="one wavenumber, one amplitude + twin")
    s.add_argument("--kr0", type=float, default=0.7)
    s.add_argument("--eps", type=float, default=1e-3)
    _common_physics_args(s)
    s.set_defaults(func=cmd_single)

    s = sub.add_parser("sweep-eps", help="amplitude sweep (proportionality test)")
    s.add_argument("--kr0", type=float, default=0.7)
    s.add_argument("--eps", type=float, nargs="+",
                   default=[1e-4, 1e-3, 1e-2])
    _common_physics_args(s)
    s.set_defaults(func=cmd_sweep_eps)

    s = sub.add_parser("sweep-k", help="wavenumber sweep (dispersion curve)")
    s.add_argument("--kr0", type=float, nargs="+",
                   default=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    s.add_argument("--eps", type=float, nargs="+", default=[1e-3])
    _common_physics_args(s)
    s.set_defaults(func=cmd_sweep_k)

    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
