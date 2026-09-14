# !/usr/bin/env python3
"""Animations of the growing mode: a GIF to embed, and a slider to explore.

A static figure shows three or four instants and asks the reader to interpolate between them. An animation shows the whole record, which is worth more for a talk than any single panel, and a slider is worth more still when someone asks what happens around a particular time.

Two outputs are written for each animation. The GIF is self-contained, embeds in slides and needs no player. The HTML carries the same frames with a slider and a play control, opens in any browser with nothing installed, and is the one to use when exploring or when sending a result to a supervisor.

Two animations are available and they differ sharply in cost.

The interface animation is free. `profiles.dat` already holds the interface at every snapshot of the run, typically a hundred and sixty of them, so the animation is made from data already exported and needs no new work.

The field animation is limited by what was exported. A field slice is the whole plane and runs to about a megabyte of text at the finest resolution, so `dump_profiles` writes only a handful by default. The animation can show no more instants than were exported, and a smooth one needs the export repeated with a larger slice count. The script says which it found rather than silently producing something jerky.

Run from `amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Python`:

    python3 python/animate.py --bundle runs/bundles/meas32 --outdir figures
    python3 python/animate.py --bundle runs/bundles/meas32 --field y_velocity
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import bundle as bundle_io
import transient
from pub_style import apply_style, PALETTE

GUIDE_LW = 0.8


def write_outputs(anim, stem, outdir, fps, dpi, make_gif=True,
                  make_html=True, html_dpi=None):
    """Write the animation as a GIF and as a self-contained HTML slider.

    The HTML embeds every frame, so it opens from a local file with no server and no network. That matters on a cluster, where the alternative is forwarding a display.
    """
    written = []
    if make_gif:
        path = os.path.join(outdir, f"{stem}.gif")
        # Pillow ships with matplotlib's dependencies, so a GIF needs nothing that a plotting environment lacks. An MP4 would need ffmpeg, which a cluster may not have.
        anim.save(path, writer=animation.PillowWriter(fps=fps), dpi=dpi)
        written.append(path)
    if make_html:
        # Every frame is embedded as a PNG, so the file grows linearly with frame count and resolution. Matplotlib refuses past twenty megabytes and drops the remaining frames with a warning that is easy to miss, leaving an animation that simply stops part way. The limit is raised and the resolution reduced so that the whole record is present.
        matplotlib.rcParams["animation.embed_limit"] = 512.0
        path = os.path.join(outdir, f"{stem}.html")
        # Written through HTMLWriter rather than to_jshtml, because the latter takes no resolution argument and renders at the figure's own, which for a hundred frames of embedded PNG runs to tens of megabytes. The writer accepts dpi directly.
        anim.save(path, writer=animation.HTMLWriter(fps=fps,
                                                    embed_frames=True,
                                                    default_mode="once"),
                  dpi=html_dpi if html_dpi is not None else dpi)
        written.append(path)
    return written


def animate_interface(profiles, t_hi, outdir, r_window=(0.825, 1.175),
                      fps=20, dpi=110, stride=1):
    """The interface shape through the linear regime.

    The ordinate is fixed at the same limits as the static figure, so that the two are read on the same scale and the animation is not a differently scaled version of a figure the reader has already seen.
    """
    z = profiles["z"]
    keep = profiles["times"] <= t_hi
    times = profiles["times"][keep][::stride]
    radii = profiles["radius"][keep][::stride]
    wavelength = profiles["wavelength"]

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    line, = ax.plot(z, radii[0], color=PALETTE[0])
    ax.axhline(1.0, ls=":", lw=GUIDE_LW, color="0.45")
    ax.set_xlabel("$z$")
    ax.set_ylabel(r"$r_{\mathrm{interface}}$")
    ax.set_xlim(0.0, float(wavelength))
    ax.set_ylim(*r_window)
    title = ax.set_title("")

    def draw(frame):
        line.set_ydata(radii[frame])
        amp = 0.5 * float(np.nanmax(radii[frame]) - np.nanmin(radii[frame]))
        title.set_text(f"Interface shape    $t={times[frame]:.2f}$    "
                       f"$a/r_0={amp:.4f}$")
        return line, title

    anim = animation.FuncAnimation(fig, draw, frames=len(times),
                                   interval=1000 / fps, blit=False)
    fig.tight_layout()
    out = write_outputs(anim, "anim_interface", outdir, fps, dpi,
                        html_dpi=max(dpi // 2, 60))
    plt.close(fig)
    return out, len(times)


def animate_field(data, profiles, name, outdir, fps=8, dpi=110, r_limit=2.0):
    """A plotfile variable on the plane, with the interface over it.

    The colour scale is fixed across every frame, from the largest magnitude in the record. Rescaling per frame would make the early frames look as vigorous as the late ones and hide the growth the animation exists to show.
    """
    times = data["times"]
    planes = data["data"]
    r, z = data["r"], data["z"]

    vmax = float(np.abs(planes).max()) or 1.0
    levels = np.linspace(-vmax, vmax, 21)
    picked = [int(np.argmin(np.abs(profiles["times"] - t))) for t in times]
    radii = profiles["radius"][picked]
    label = name.replace("_", " ")

    fig, ax = plt.subplots(figsize=(3.6, 5.0))
    ax.set_xlabel("$r$")
    ax.set_ylabel("$z$")
    ax.set_xlim(0.0, min(r_limit, float(r.max())))
    ax.set_ylim(float(z.min()), float(z.max()))
    ax.grid(False)
    contour = ax.contourf(r, z, planes[0].T, levels=levels, cmap="RdBu_r",
                          extend="both")
    interface, = ax.plot(radii[0], profiles["z"], "-", color="k")
    fig.colorbar(contour, ax=ax, pad=0.03, label=label)
    title = ax.set_title("")

    def draw(frame):
        # A filled contour set cannot be updated in place, so the old one is removed and redrawn. With a few dozen frames this is fast enough and avoids the artefacts of reusing collections.
        for art in ax.collections[:]:
            art.remove()
        ax.contourf(r, z, planes[frame].T, levels=levels, cmap="RdBu_r",
                    extend="both")
        interface.set_xdata(radii[frame])
        ax.add_line(interface)
        title.set_text(f"{label}    $t={times[frame]:.2f}$")
        return ()

    anim = animation.FuncAnimation(fig, draw, frames=len(times),
                                   interval=1000 / fps, blit=False)
    fig.tight_layout()
    out = write_outputs(anim, f"anim_field_{name}", outdir, fps, dpi,
                        html_dpi=max(dpi // 2, 60))
    plt.close(fig)
    return out, len(times)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--field", default=None,
                    help="also animate this field; default every one found")
    ap.add_argument("--no-field", action="store_true", dest="no_field",
                    help="interface only, which needs no extra export")
    ap.add_argument("--r0", type=float, default=1.0)
    ap.add_argument("--linear-bound", type=float, default=0.1,
                    dest="linear_bound")
    ap.add_argument("--fps", type=int, default=20)
    ap.add_argument("--dpi", type=int, default=110)
    ap.add_argument("--stride", type=int, default=2,
                    help="keep every nth interface snapshot, for a lighter file")
    ap.add_argument("--r-min", type=float, default=0.825, dest="r_min")
    ap.add_argument("--r-max", type=float, default=1.175, dest="r_max")
    ap.add_argument("--outdir", default="figures")
    args = ap.parse_args(argv)

    apply_style()
    os.makedirs(args.outdir, exist_ok=True)

    profiles = bundle_io.read_profiles(
        os.path.join(args.bundle, "profiles.dat"))

    # The animation stops where the linear description does, matching the static figures. An animation running on into break-up would be a different figure with a different purpose.
    record = os.path.join(args.bundle, "amplitude.dat")
    t_hi = float(profiles["times"].max())
    if os.path.isfile(record):
        times, amps, _ = bundle_io.read_amplitude(record)
        t_hi, reached = transient.linear_end(times, amps, r0=args.r0,
                                             bound=args.linear_bound)
        if reached:
            print(f"linear range ends at t = {t_hi:.3f}; "
                  f"later snapshots are excluded")

    out, n = animate_interface(profiles, t_hi, args.outdir,
                               r_window=(args.r_min, args.r_max),
                               fps=args.fps, dpi=args.dpi, stride=args.stride)
    print(f"interface animation, {n} frames:")
    for path in out:
        print(f"    {path}")

    if args.no_field:
        return 0

    found = bundle_io.find_field_bundles(args.bundle)
    wanted = [args.field] if args.field else list(found)
    for name in wanted:
        if name not in found:
            print(f"no field bundle for {name} in {args.bundle}")
            continue
        data = bundle_io.read_field(found[name])
        n_slice = data["times"].size
        if n_slice < 8:
            print(f"\n{name}: only {n_slice} slices exported, which is too "
                  f"few for a smooth animation. Re-export with more, for "
                  f"example\n    ./bin/dump_profiles <rundir> <wavelength> "
                  f"{args.bundle} {name} 60")
        out, n = animate_field(data, profiles, name, args.outdir,
                               fps=max(args.fps // 2, 4), dpi=args.dpi)
        print(f"{name} animation, {n} frames:")
        for path in out:
            print(f"    {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
