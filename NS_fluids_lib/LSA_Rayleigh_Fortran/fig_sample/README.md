# Sample figures

These files show what each figure of the analysis should look like. They are a reference for appearance, not a set of results: figures generated from a sweep are written to `figures/`, which this directory does not replace.

Everything at the top level was drawn by the current code from `validation_data/`, which holds the reference calculation at $kr_0 = 0.7$. Regenerate the whole set at any time with

```bash
python3 plotting/make_report_figures.py --datadir validation_data --archive validation_data/validation_kr0_0.7.json --outdir fig_sample
python3 plotting/plot_results.py --json validation_data/validation_kr0_0.7.json --kr0 0.7 --t-fit 4.5 --outdir fig_sample
```

in the Fortran tree, or the same two commands with `python/` in place of `plotting/` in the Python tree. Both draw from stored records rather than from plotfiles, so neither needs `yt` and neither needs the run directories.

| Figure | Shows |
|---|---|
| `fig1_validation` | the fitted growth rate beside the local rate, which together justify the fit window |
| `fig2_replication` | the two solver builds on the same configuration, with the distribution of their disagreement |
| `fig3_proportionality` | the amplitude sweep collapsed by $\varepsilon$, and the departure from the mean |
| `fig4_dispersion` | the predicted dispersion relation with the measured point on it |
| `growth_kr0_0p7` | interface mode amplitude against time, with the fit and both theoretical rates |
| `local_sigma_kr0_0p7` | the local growth rate; here still rising, which is how a run stopped inside the transient shows itself |
| `ohnesorge_surface` | growth rate over the wavenumber and Ohnesorge plane, as a surface |
| `ohnesorge_contour` | the same as contours, with the locus of the most unstable wavenumber |
| `ohnesorge_peak_drift` | drift of that locus away from the inviscid value of 0.697 |
| `schematic_twin_run` | the twin-run architecture; the one figure carrying no data |

The four Ohnesorge and schematic figures need no calculation at all: the first three are evaluated from the dispersion relation and the last is drawn, so both can be regenerated on any machine with matplotlib. `fig4_dispersion` carries the single measured point that `validation_data` supports. Once a sweep exists, `plot_dispersion_sweep.py` draws the whole measured curve from `runs/records/` and should be preferred.

`local_sigma_kr0_0p7` and the right panel of `fig1_validation` present the same quantity, so a document wants one of the two rather than both.

These samples come from a calculation run for 2.6 e-folding times, which is now known to stop inside the modal transient and to bias the fitted rate low by about three per cent. They are kept because they show correctly what each figure looks like, and because `local_sigma` shows the diagnostic symptom plainly: a curve still climbing at the right-hand edge. Runs are now advanced for 5.0 e-folding times, and these samples will be replaced once those complete.

## `legacy/`

Four figures from an earlier calculation, kept for one reason: they show the intended appearance of the two figures that the archived data cannot reproduce.

`growth.png` and `local_sigma.png` have direct equivalents above, drawn from the archived record. `interface.png` and `field_y_velocity.png` do not, because they need the interface radius profile $r(z)$ and a plotfile variable on the $(r,z)$ plane, and `validation_data` carries neither; it holds the amplitude record only. Both are produced by the current code as soon as a run exists, through the figure-data bundle:

```bash
./bin/dump_amplitude runs/k0.7_eps0.02 8.976 runs/bundles/k0.7/amplitude.dat
./bin/dump_profiles runs/k0.7_eps0.02 8.976 runs/bundles/k0.7 y_velocity
python3 plotting/plot_results.py --bundle runs/bundles/k0.7 --kr0 0.7 --outdir figures
```

`scripts/run_all.sh` and `scripts/run_nsmfp.sh` perform the first two steps for every case automatically, so in normal use only the last command is issued.

These four predate the present house style and differ from current output in typography, line width and colour. They are retained for content, not for style, and should not be cited in a document or used as a style reference.
