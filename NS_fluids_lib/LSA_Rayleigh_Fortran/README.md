# LSA_Rayleigh_Fortran

Linear stability analysis of the Rayleigh–Plateau capillary instability by the Navier–Stokes-based mean flow perturbation method of Ranjan, Unnikrishnan and Gaitonde (*J. Comput. Phys.* **403**, 109076, 2020), applied to the two-phase solver in `amrex_implicit_interfaces/NS_fluids_lib`.

The method, its adaptation to an interfacial instability, and the measurements obtained are described in `analysis_report.typ`. This file covers installation and use. A Python implementation of the same analysis accompanies this one; the two are independent and serve as mutual checks.

The Fortran implementation requires no Python, and depends on nothing beyond a Fortran 2008 compiler with LAPACK and BLAS. It is therefore the preferable choice where installing a scientific Python stack on a cluster is unwelcome, or where a single language across solver and analysis is convenient.

## Placement and invocation

This directory is intended to sit inside the solver tree:

```
amrex_implicit_interfaces/
    amrex-master/                  vendored AMReX, supplied with the repository
    NS_fluids_lib/
        amr2d.gnu.MPI.ex           the solver executable
        run2d/                     decks distributed with the solver
        LSA_Rayleigh_Fortran/      this directory, holding runs/ and figures/
        LSA_Rayleigh_Python/       the Python implementation
```

All commands below are issued from `amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Fortran`, and every file produced is written inside that directory: simulation output to `runs`, figures to `figures`. Nothing is written above it, so the whole study may be moved, archived or deleted by acting on this one directory. The only path referred to outside it is the solver executable, which is read and never written.

## Building

```bash
sudo apt-get install gfortran liblapack-dev libblas-dev   # Debian, Ubuntu
brew install gcc                                          # macOS
```

On macOS, LAPACK is provided through the Accelerate framework:

```bash
make LAPACK_LIBS="-framework Accelerate"
```

Then, from this directory:

```bash
make          # produces bin/make_inputs, bin/run_lsa, bin/dump_amplitude,
              # bin/dump_profiles
make test     # 169 assertions across five test programs
```

The test suite must be run from this directory, since the deck-generation and end-to-end tests read `inputs/inputs_rayleigh_template.txt` by relative path.

The solver itself is built from `NS_fluids_lib`:

```bash
cd .. && make -j16
```

AMReX is vendored in the repository at `../amrex-master`; no separate clone is required. Should the link fail with an undefined reference to `amrex::FillPatchTower`, an interrupted build has left an object file of zero length, which `make` treats as current:

```bash
find tmp_build_dir -name '*.o' -size 0 -delete && make -j16
```

## Procedure

Four documents describe the procedure at different depths. `Running_the_Analysis_Fortran.md` is the step-by-step guide, running from a machine on which nothing has been installed through to a measured growth rate and the full sweep, and is the place to start. `Running_the_Analysis_HPC_Fortran.md` covers what differs on the FSU Research Computing Center cluster, where a job array advances the cases of the sweep concurrently, and states what must hold when part of the study is run on a workstation and part on the cluster. `further_analysis.md` is the runbook for the remaining calculations, with the expected result and an acceptance criterion at each stage. The outline below is the shortest form of the same sequence, and presumes that the solver and the tools of this directory have already been built.

```bash
# 1. matched decks: a perturbed case and its unperturbed twin
./bin/make_inputs single --kr0 0.7 --eps 2e-2 --outdir runs/decks

# 2. advance both cases, then extract the growth rate
./bin/run_lsa analyse --perturbed runs/k0.7_eps0.02 \
    --base runs/k0.7_base --kr0 0.7 --skip 60 --t-min 4.0

# 3. proportionality of the response, which must precede any measurement
./bin/run_lsa linearity --base runs/k0.7_base --kr0 0.7 \
    --case 1e-2=runs/k0.7_eps0.01 \
    --case 2e-2=runs/k0.7_eps0.02 \
    --case 4e-2=runs/k0.7_eps0.04 --t-min 4.0

# 4. interface mode amplitude as plain text, for plotting elsewhere
./bin/dump_amplitude runs/k0.7_eps0.02 8.976 runs/bundles/k0.7/amplitude.dat
./bin/dump_profiles  runs/k0.7_eps0.02 8.976 runs/bundles/k0.7 y_velocity
```

Decks are generated rather than edited by hand, because a perturbed deck and its twin must differ in exactly one functional line, `ns.radblob`; that invariant is enforced by a test rather than left to inspection.

The argument `--skip` is not optional for velocity-based subspaces. Because the perturbation is interfacial, both calculations begin from rest and the perturbation field vanishes identically at the initial instant, which makes the amplitude ranking degenerate. Section 2.4 of `analysis_report.typ` explains the consequence and the remedy.

Two further points bear on how runs are managed. The solver does not honour `amr.plot_file`; with `ns.visual_nddata_format = 1` it writes directories named `nddataPLT00000023`, which is the default prefix assumed throughout and may be overridden with `--prefix`. Checkpoint restart accepts command-line overrides, so a calculation may be advanced beyond the end of its deck without editing it:

```bash
mpirun -n 16 ../amr2d.gnu.MPI.ex inputs amr.check_int=200
mpirun -n 16 ../amr2d.gnu.MPI.ex inputs amr.check_int=200 \
    amr.restart=chk00400 max_step=1600 stop_time=20.0
```

The checkpoint interval must be smaller than the number of steps achievable within one scheduler allocation, or no checkpoint is written and each restart silently repeats from the beginning.

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

The analysis tools write plain text and JSON; figures are produced by the Python scripts in `plotting/`, which are copies of those in the Python tree and are included so that this directory is self-contained.

```bash
python3 plotting/make_report_figures.py --datadir validation_data \
    --archive validation_data/validation_kr0_0.7.json --outdir figures
```

These follow a fixed house style defined in `plotting/pub_style.py`: Computer Modern typography, inward ticks on all four sides, no axis padding, a data line width of 1.0 point, and PNG, PDF and SVG output at 600 dpi. They require `numpy`, `scipy` and `matplotlib`, and read the plain text written by `dump_amplitude` and `dump_profiles`, so `yt` is not needed.

Four figures describe a case: the growth curve, the local growth rate, the interface shape and a field slice with the interface over it. The last two need the interface profile and a plotfile variable on the plane, which no amplitude record carries, so the sweep scripts export a figure-data bundle for every case while the plotfiles are still present. A bundle is a few hundred kilobytes against tens of gigabytes of plotfiles, and `plot_results.py --bundle` draws all four figures from it. Where a Python installation is unavailable altogether, a bundle may be copied to another machine and plotted there.

## Contents

```
LSA_Rayleigh_Fortran/
    Makefile
    README.md
    Running_the_Analysis_Fortran.md         step-by-step procedure on a workstation
    Running_the_Analysis_HPC_Fortran.md     the same study on the FSU RCC cluster
    analysis_report.typ           method, measurements and discussion, Typst source
    references.bib                bibliography, DOIs for peer-reviewed entries
    further_analysis.md           runbook for the remaining calculations
    FINDINGS_FROM_REAL_RUNS.md    what the solver runs established
    src/
        dispersion_mod.f90        exact Rayleigh-Plateau dispersion relation
        dmd_mod.f90               SVD-based decomposition through LAPACK
        plotfile_mod.f90          native AMReX plotfile reader
        nsmfp_mod.f90             twin subtraction and subspace assembly
        interface_mode_mod.f90    interface-amplitude observable
        proportionality_mod.f90   linearity test
        make_inputs_mod.f90       deck generation
        make_inputs.f90           deck generator, command line
        run_lsa.f90               analysis driver, command line
        dump_amplitude.f90        amplitude record as plain text
        dump_profiles.f90         interface profiles and field slices
        bundle_mod.f90            the portable figure-data format
    plotting/                     the figure scripts, mirroring the Python tree
        plot_dmd.py               spectrum, amplitudes and convergence
        plot_mode_shape.py        eigenvector against the analytic form
        plot_ohnesorge.py         growth rate over the parameter plane
        plot_pinchoff.py          break-up, and the linear extrapolation
        benchmark_jvp.py          cost, with the machine recorded
        plot_schematic.py         the twin-run architecture, drawn
        machine_info.py           machine and run provenance
    tests/                        169 assertions, five programs
    plotting/                     figure scripts and house style, Python
    full_analysis.txt             every calculation of the study, one per line
    scripts/                      sweep runner and single-case driver
    inputs/                       deck template
    validation_data/              measured amplitude records
    fig_sample/                   sample figures, for reference rather than results
        legacy/                   earlier figures, kept to show intended appearance
    figures/                      figures generated from a sweep, created on first use
```

## Verification

The test suite comprises 169 assertions and completes in a few seconds.

The dispersion relation is checked in 43 assertions. The modified Bessel functions are compared against trapezoidal quadrature of the integral representation rather than against the same power series, which would be circular, and the assembled growth rate is compared against the utility distributed with the solver.

The decomposition is checked in 20 assertions, recovering planted growth rates and frequencies to approximately $10^{-8}$ for mixed stationary and oscillatory spectra, together with rank truncation and input validation.

The proportionality test is checked in 38 assertions, including the separation of nonlinear contamination from round-off and the twin invariant on generated decks.

The end-to-end path is checked in 37 assertions. Synthetic plotfiles written in the genuine on-disk format are read back through the production parser, differenced and decomposed, recovering the exact Rayleigh growth rate to $10^{-6}$. The planted base-state drift is two orders of magnitude larger than the perturbation, so the test fails if the twin subtraction is omitted or the snapshots are mispaired.

The plotfile parser carries the greatest risk of any component here, since no library equivalent to `yt` exists for Fortran. Confidence rests on three-way interoperability: a plotfile written by the Fortran writer is read correctly by `yt`; a plotfile written by the Python implementation is read correctly by this parser; and the writer and parser round-trip in both byte orders. During the analysis reported in `analysis_report.typ` this parser produced the amplitude records on a machine where `yt` was unavailable, which exercised the dependency-free path in earnest rather than only in test.

The analysis report, `analysis_report.typ`, and its bibliography, `references.bib`, sit one level above this directory, alongside the sibling implementation. There is one report because there is one set of results: the two trees implement the same method in different languages and are checked against each other, so a copy in each would be two documents free to drift apart.
