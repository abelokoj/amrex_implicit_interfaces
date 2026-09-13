# LSA_Rayleigh_Python

Linear stability analysis of the Rayleigh–Plateau capillary instability by the Navier–Stokes-based mean flow perturbation method of Ranjan, Unnikrishnan and Gaitonde (*J. Comput. Phys.* **403**, 109076, 2020), applied to the two-phase solver in `amrex_implicit_interfaces/NS_fluids_lib`.

The method and its adaptation to an interfacial instability are described in `analysis_report.typ`, which also reports the measurements. This file covers installation and use. A Fortran implementation of the same analysis accompanies this one; the two are independent and serve as mutual checks.

## Placement and invocation

This directory is intended to sit inside the solver tree:

```
amrex_implicit_interfaces/
    amrex-master/                  vendored AMReX, supplied with the repository
    NS_fluids_lib/
        amr2d.gnu.MPI.ex           the solver executable
        run2d/                     decks distributed with the solver
        LSA_Rayleigh_Python/       this directory
        LSA_Rayleigh_Fortran/      the Fortran implementation
        runs/                      simulation output, created on first use
```

All commands below are issued from `amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Python`, and every file produced is written inside that directory: simulation output to `runs`, figures to `figures`. Nothing is written above it, so the whole study may be moved, archived or deleted by acting on this one directory. The only path referred to outside it is the solver executable, which is read and never written.

## Requirements

```bash
pip install numpy scipy matplotlib yt pytest
```

Only `numpy` and `scipy` are needed for the dispersion relation and the decomposition. Reading AMReX plotfiles requires `yt`; `matplotlib` is needed only for figures. The test suite reports 83 assertions when `yt` is present, and skips the tests needing it otherwise.

The Fortran implementation requires no Python at all, and is preferable where installing this stack on a cluster is unwelcome.

## Procedure

Four documents describe the procedure at different depths. `Running_the_Analysis_Python.md` is the step-by-step guide, running from a machine on which nothing has been installed through to a measured growth rate and the full sweep, and is the place to start. `Running_the_Analysis_HPC_Python.md` covers what differs on the FSU Research Computing Center cluster, where a job array advances the cases of the sweep concurrently, and states what must hold when part of the study is run on a workstation and part on the cluster. `further_analysis.md` is the runbook for the remaining calculations, with the expected result and an acceptance criterion for each stage. The outline below is the shortest form of the same sequence, and presumes that the solver has already been built.

```bash
# 1. matched decks: a perturbed case and its unperturbed twin
python3 python/make_inputs.py single --kr0 0.7 --eps 2e-2 \
    --outdir runs/decks

# 2. advance both cases, then extract the growth rate
python3 python/run_lsa.py analyse \
    --perturbed runs/k0.7_eps0.02 --base runs/k0.7_base \
    --kr0 0.7 --skip 60 --t-min 4.0

# 3. proportionality of the response, which must precede any measurement
python3 python/run_lsa.py linearity --base runs/k0.7_base --kr0 0.7 \
    --case 1e-2=runs/k0.7_eps0.01 \
    --case 2e-2=runs/k0.7_eps0.02 \
    --case 4e-2=runs/k0.7_eps0.04 --t-min 4.0
```

The decks are generated rather than edited by hand, because a perturbed deck and its twin must differ in exactly one functional line, `ns.radblob`; that invariant is enforced by a test.

The argument `--skip` is not optional for velocity-based subspaces. Because the perturbation is interfacial, both calculations begin from rest and the perturbation field vanishes identically at the initial instant, which makes the amplitude ranking degenerate. Section 2.4 of `analysis_report.typ` explains the consequence and the remedy.

## Running the whole study in one sweep

Every calculation of the study is listed in `full_analysis.txt`, one per line, with comments describing what each group tests and what result to expect. Lines beginning with `#` are ignored, so a block may be commented out to run a subset without altering the field values.

```bash
./scripts/run_all.sh                          # full_analysis.txt, one rank
./scripts/run_all.sh full_analysis.txt 16   # sixteen ranks
./scripts/run_all.sh my_subset.txt 8        # a subset of the cases
```

The script generates the deck for each case, advances it to completion by restarting from the most recent checkpoint as often as necessary, and extracts the interface mode amplitude to `runs/records/<label>.dat`. Progress is appended to `runs/run_all.log` and a summary table is written to `runs/summary.tsv`.

The sweep is intended to be started once and left to run, and it is safe to interrupt: a case whose record already exists is skipped, and a partly advanced case resumes from its checkpoint, so the script may simply be run again. A failure in one case is recorded and does not stop the others, so that a single unstable configuration cannot cost the whole sweep.

Cost scales roughly as the cube of the resolution, since the capillary time-step limit varies as $\Delta x^{3/2}$. Order the file so that the cases most needed finish first, because the runner proceeds from top to bottom.

## Figures

Figures follow a fixed house style defined in `python/pub_style.py`: Computer Modern typography, inward ticks on all four sides, no axis padding, and PNG, PDF and SVG output at 600 dpi. Settings are not restated in the plotting scripts.

```bash
# figures for the report, from the stored amplitude records
python3 python/make_report_figures.py --datadir validation_data \
    --archive validation_data/validation_kr0_0.7.json --outdir figures

# all four figures for a case, from the exported figure data
python3 python/plot_results.py --bundle runs/bundles/k0.7_eps0.02 \
    --kr0 0.7 --outdir figures
```

Four figures describe a case: the growth curve, the local growth rate, the interface shape and a field slice with the interface over it. The last two need the interface profile and a plotfile variable on the plane, which no amplitude record carries, so the sweep scripts export a figure-data bundle for every case while the plotfiles are still present. A bundle is a few hundred kilobytes of plain text against tens of gigabytes of plotfiles, and `--bundle` draws all four figures from it without `yt` and without the plotfiles. Passing `--rundir` instead reads the plotfiles directly and produces identical output; `--json` accepts a stored record and draws the two amplitude figures only.

## Contents

```
LSA_Rayleigh_Python/
    README.md
    Running_the_Analysis_Python.md         step-by-step procedure on a workstation
    Running_the_Analysis_HPC_Python.md     the same study on the FSU RCC cluster
    analysis_report.typ           method, measurements and discussion, Typst source
    references.bib                bibliography, DOIs for peer-reviewed entries
    further_analysis.md           runbook for the remaining calculations
    FINDINGS_FROM_REAL_RUNS.md    what the solver runs established
    python/
        dispersion.py             exact Rayleigh-Plateau dispersion relation
        dmd.py                    SVD-based dynamic mode decomposition
        nsmfp.py                  plotfile reading and twin subtraction
        interface_mode.py         interface-amplitude observable
        proportionality.py        linearity test
        make_inputs.py            deck generation
        run_lsa.py                analysis driver
        dump_amplitude.py         amplitude record as plain text
        dump_profiles.py          interface profiles and field slices
        bundle.py                 the portable figure-data format
        plot_dmd.py               spectrum, amplitudes and convergence
        plot_mode_shape.py        eigenvector against the analytic form
        plot_ohnesorge.py         growth rate over the parameter plane
        plot_pinchoff.py          break-up, and the linear extrapolation
        benchmark_jvp.py          cost, with the machine recorded
        plot_schematic.py         the twin-run architecture, drawn
        machine_info.py           machine and run provenance
        plot_results.py           figures for a single calculation
        make_report_figures.py    figures reproduced in the report
        pub_style.py              figure house style
    tests/                        74 assertions under pytest
    full_analysis.txt             every calculation of the study, one per line
    scripts/                      sweep runner and single-case driver
    inputs/                       deck template
    validation_data/              measured amplitude records
    fig_sample/                   sample figures, for reference rather than results
        legacy/                   earlier figures, kept to show intended appearance
    figures/                      figures generated from a sweep, created on first use
```

## Relation to the Fortran implementation

The two implementations differ in every major dependency: NumPy and SciPy against LAPACK, `yt` against a hand-written plotfile parser, and `scipy.special.iv` against a truncated Bessel series. Three checks establish their equivalence. Both reproduce the utility distributed with the solver, `run2d/rayleigh_capillary_growth.F90`, to machine precision. Both deck generators produce identical parameters for identical input, verified at a wavenumber not used during development. Their plotfile readers and writers are mutually compatible, and a plotfile written by either is read correctly by `yt`.

The two are interchangeable in use. The Fortran tools write the same quantities in plain text and JSON, which the plotting scripts here consume directly; this was exercised during the analysis reported in `analysis_report.typ`, where the amplitude records were produced by the Fortran reader on a machine without `yt` installed.

The analysis report, `analysis_report.typ`, and its bibliography, `references.bib`, sit one level above this directory, alongside the sibling implementation. There is one report because there is one set of results: the two trees implement the same method in different languages and are checked against each other, so a copy in each would be two documents free to drift apart.
