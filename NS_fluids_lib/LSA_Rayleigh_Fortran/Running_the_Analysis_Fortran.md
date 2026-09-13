# Running the analysis, step by step

This document carries the procedure from a machine on which nothing has yet been installed to a measured growth rate, and thereafter to the full sweep. It assumes no prior familiarity with either the solver or this directory. `README.md` describes what the tools are; `analysis_report.typ` describes what the method does and what it has established; the present file describes only the order in which commands are issued and what each is expected to produce.

Every command below is written out in full relative to the solver tree. Where a step must be issued from a particular directory, that directory is named at the head of the step, and the prompt shown in the listing states it.

## 0. Notation and the one assumption

A single environment variable stands for the top of the solver tree, so that no path in this document is ambiguous. Set it once per shell session, adjusting the location to the machine in use.

```bash
export AII=$HOME/amrex_implicit_interfaces
```

Three directories are referred to throughout:

| Symbol | Path | Contents |
|---|---|---|
| solver tree | `$AII` | the repository, including the vendored AMReX |
| solver directory | `$AII/NS_fluids_lib` | the solver source, its executable, and `run2d` |
| this directory | `$AII/NS_fluids_lib/LSA_Rayleigh_Fortran` | the analysis code, and all of its output |

Every file produced by this work is written inside this directory. Simulation output goes to `./runs`, created on first use, and figures to `./figures`. Nothing is written above this directory, to the home directory, or anywhere else on the system, so the whole study may be moved, archived or deleted by acting on this one directory. The only path referred to outside it is the solver executable itself, which is read and never written.

No program used here resides outside `$AII`. The executables `bin/make_inputs`, `bin/run_lsa` and `bin/dump_amplitude` do not exist until Step 5 builds them; if they are absent, the build has not yet been performed, and no installed copy elsewhere on the system is expected or required.

## 1. Obtain the solver

```bash
cd $(dirname $AII)
git clone https://github.com/msussman42/amrex_implicit_interfaces.git
cd $AII
git checkout 0b8df42
```

Should the clone fail on a file whose name contains quotation marks, `"u_vector".xml`, the option `git -c core.protectNTFS=false clone` succeeds on a Linux filesystem.

AMReX is vendored in the repository at `$AII/amrex-master` and requires no separate clone.

The checkout of `0b8df42` is deliberate and is the reason this step is not simply a clone. That commit corrects the defect recorded in Section 2 of `FINDINGS_FROM_REAL_RUNS.md`, under which every fluid-only calculation terminates during initialisation with `expecting num_FSI_outer_sweeps>=2`. The defect is present again at `b014c1f`, which is the current head, so a plain clone reproduces it. Working from `0b8df42` has been confirmed to run correctly and requires no modification of the solver source, which is why it is preferred to correcting a later commit by hand: a pinned commit is stated in one line of a thesis, whereas a local edit must be described, justified and shown not to have altered the physics.

The commit in use should be recorded with every result, since the solver is not versioned by release. Before moving to a later head, repeat Step 4 to confirm that the assertion has been corrected there.

## 2. Install the compiler and the libraries

```bash
sudo apt-get install gfortran g++ libopenmpi-dev openmpi-bin \
    liblapack-dev libblas-dev            # Debian, Ubuntu

brew install gcc open-mpi lapack         # macOS
```

The analysis tools require a Fortran 2008 compiler with LAPACK and BLAS, and nothing further. Python is not required for the analysis itself; it is required only for the figures of Step 9, and those may equally be produced on another machine from the plain-text records.

On macOS, LAPACK is supplied through the Accelerate framework, so the build of Step 5 takes the form `make LAPACK_LIBS="-framework Accelerate"`.

## 3. Build the solver

Issued from the solver directory.

```bash
cd $AII/NS_fluids_lib
make -j16
```

The expected product is `$AII/NS_fluids_lib/amr2d.gnu.MPI.ex`, of approximately 109 MB. The build takes tens of minutes.

Confirm that this is what was built, because the name depends on settings in `GNUmakefile` and the analysis assumes the MPI build:

```bash
cd $AII/NS_fluids_lib
ls -la amr2d.gnu.*.ex
```

A report of `amr2d.gnu.DEBUG.ex is up to date` means `DEBUG = TRUE` is set, in which case `make` leaves any existing `amr2d.gnu.MPI.ex` untouched and possibly stale, while appearing to have succeeded. Set `DEBUG = FALSE` and `USE_MPI = TRUE` in `GNUmakefile` and rebuild, or pass the debug executable explicitly as the fourth argument of the driver script in Step 6. A debug build is correct but far slower, so it is unsuitable for the sweep.

Should the link fail with an undefined reference to `amrex::FillPatchTower`, an interrupted build has left an object file of zero length, which `make` treats as current. The remedy is to remove such files and rebuild:

```bash
cd $AII/NS_fluids_lib
find tmp_build_dir -name '*.o' -size 0 -delete && make -j16
```

## 4. Confirm that the solver runs before going further

Issued from the distributed deck directory. This step exists so that a failure of the solver is distinguished from a failure of the analysis, which is otherwise a costly confusion.

```bash
cd $AII/NS_fluids_lib/run2d
mpirun -n 1 $AII/NS_fluids_lib/amr2d.gnu.MPI.ex inputs.growthrate.LSA
```

Time-stepping should begin within a few seconds, with lines beginning `STEP` appearing on standard output. The calculation may be interrupted once that is seen. Termination with `expecting num_FSI_outer_sweeps>=2` means the working tree is not at `0b8df42`; return to Step 1. The assertion prevents every fluid-only calculation from starting, so nothing further is worth attempting until it is gone.

The launcher itself must also be checked, because a failure of this kind does not announce itself as an error:

```bash
mpirun -n 2 hostname
```

Two lines should appear. A message reading `No PMIx server was reachable, but a PMI1/2 was detected. If srun is being used to launch application, 2 singletons will be started` means the launcher cannot reach a process manager and starts independent single-rank processes instead of one job of the requested size. Each would then advance the entire problem in the same directory and overwrite the output of the others, so the run appears to succeed while producing corrupted results. The driver scripts of Steps 6 and 8 now test for this and refuse to start, but it is better to establish it here. The usual cause is that `mpirun` comes from a different installation than the libraries the solver was linked against; `which mpirun` and `ldd $AII/NS_fluids_lib/amr2d.gnu.MPI.ex | grep -i mpi` identify the mismatch. Until it is resolved, run on a single rank, which is correct if slow.

## 5. Build and test the analysis tools

Issued from this directory. The test suite reads `inputs/inputs_rayleigh_template.txt` by a relative path and must therefore be run from here rather than from `scripts/`.

```bash
cd $AII/NS_fluids_lib/LSA_Rayleigh_Fortran
chmod +x scripts/*.sh
make
make test
```

The `chmod` is required because an archive does not always preserve the executable bit through extraction, and without it the driver scripts of Steps 6 and 8 fail with `Permission denied`. Running them as `bash scripts/run_nsmfp.sh` works regardless and is the alternative if the permissions cannot be changed.

`make` produces `bin/make_inputs`, `bin/run_lsa` and `bin/dump_amplitude`. `make test` builds four test programs and reports 138 assertions, completing in a few seconds. Both commands must succeed before any calculation is attempted, since the test suite exercises the dispersion relation, the decomposition, the proportionality test and the complete plotfile path.

## 6. A first calculation

Issued from this directory. This is the reference case: wavenumber $kr_0 = 0.7$, amplitude $\varepsilon = 2\times10^{-2}$, 8 cells per column radius. The single script performs the whole of Steps 6.1 to 6.3, and is the shortest route to a result.

```bash
cd $AII/NS_fluids_lib/LSA_Rayleigh_Fortran
./scripts/run_nsmfp.sh 0.7 2e-2 8
```

The arguments are the wavenumber, the amplitude, the number of MPI ranks and, optionally, the cells per column radius, which defaults to the reference value of 8; the solver is located automatically at `$AII/NS_fluids_lib/amr2d.gnu.MPI.ex`. At the reference resolution this is 1800 steps on a 32 by 72 grid, which takes some hours on a single core and some minutes on sixteen ranks. Raising the resolution is expensive and Section 12 gives the scaling before it is attempted.

The steps the script performs are set out below, since they are also the form in which an individual case is run by hand.

### 6.1 Generate the matched decks

```bash
cd $AII/NS_fluids_lib/LSA_Rayleigh_Fortran
./bin/make_inputs single --kr0 0.7 --eps 2e-2 \
    --cells-per-r0 8 --rmax 4.0 --n-periods 5.0 --snapshots 120 \
    --outdir runs/decks
```

Two decks are written, `runs/decks/inputs.k0.7.eps0.02` and `runs/decks/inputs.k0.7.base`. Decks are generated rather than edited by hand because a perturbed deck and its twin must differ in exactly one functional line, `ns.radblob`; the twin subtraction stands in for the NS-MFP body force and is valid only under that condition, which a test enforces.

### 6.2 Advance the two calculations

Each calculation is advanced from its own directory, since the solver writes its plotfiles and checkpoints into the working directory.

```bash
mkdir -p $AII/NS_fluids_lib/LSA_Rayleigh_Fortran/runs/k0.7_eps0.02
cd $AII/NS_fluids_lib/LSA_Rayleigh_Fortran/runs/k0.7_eps0.02
cp ../decks/inputs.k0.7.eps0.02 inputs
mpirun -n 8 $AII/NS_fluids_lib/amr2d.gnu.MPI.ex inputs \
    mac.mac_abs_tol=1.0e-10 mac.visc_abs_tol=1.0e-10 amr.check_int=200
```

The twin is advanced identically, from `$AII/NS_fluids_lib/LSA_Rayleigh_Fortran/runs/k0.7_base` with the deck `inputs.k0.7.base`, and on the same number of ranks. Altering the domain decomposition between the two perturbs round-off-level behaviour, which is the scale at which the measurement is made.

The twin is required only for the velocity subspace of Step 10 and for the $\lVert Q' \rVert$ form of the proportionality test. The growth rate itself is obtained from the perturbed calculation alone, so a study confined to growth rates may omit every twin and halve its cost.

Two points govern how a long calculation is managed. The solver does not honour `amr.plot_file`: with `ns.visual_nddata_format = 1` it writes directories named `nddataPLT00000023` rather than `plt00023`, which is the prefix the analysis tools assume. And the checkpoint interval must be smaller than the number of steps achievable within one scheduler allocation, or no checkpoint is written and each restart silently repeats from the beginning while accumulating plotfiles give every appearance of progress. A calculation is resumed with

```bash
cd $AII/NS_fluids_lib/LSA_Rayleigh_Fortran/runs/k0.7_eps0.02
mpirun -n 8 $AII/NS_fluids_lib/amr2d.gnu.MPI.ex inputs \
    amr.check_int=200 amr.restart=chk00400 max_step=1600 stop_time=20.0
```

### 6.3 Extract the growth rate

Issued from this directory. The wavelength argument is `ns.yblob` of the deck, which for $kr_0 = 0.7$ is 8.976; it may be read back from the deck rather than recalled, with `grep '^ns.yblob' runs/k0.7_eps0.02/inputs`.

```bash
cd $AII/NS_fluids_lib/LSA_Rayleigh_Fortran
./bin/dump_amplitude runs/k0.7_eps0.02 8.976 runs/records/k0.7.dat

./bin/run_lsa analyse \
    --perturbed runs/k0.7_eps0.02 --base runs/k0.7_base \
    --kr0 0.7 --skip 60 --t-min 4.0 \
    --save runs/result_k0.7.json
```

The expected result at this resolution is a fitted rate near 0.319 with a coefficient of determination above 0.999, lying some 4 per cent below the viscous estimate of 0.3289 and below the inviscid value of 0.3433. A value above the inviscid figure is physically impossible and indicates an error in the deck, most often a disagreement between `ns.yblob` and the axial domain height.

The argument `--skip` is not optional for velocity-based subspaces. Because the perturbation is interfacial, both calculations begin from rest and the perturbation field vanishes identically at the initial instant, which renders the amplitude ranking degenerate; Section 2.4 of `analysis_report.typ` explains the consequence. The argument `--t-min` excludes the startup transient, without which the fitted rate is a substantial underestimate: a fit over the whole record returns 0.173 rather than 0.319.

## 7. The proportionality test, which precedes any quoted result

Issued from this directory. The measurement presumes that the response is linear in the imposed amplitude. Since the snapshots here are constructed as actions of the linear operator rather than drawn from nonlinear simulation data, that presumption is load-bearing and must be verified rather than assumed.

Three further calculations at $\varepsilon \in \{5\times10^{-3}, 1\times10^{-2}, 4\times10^{-2}\}$ are advanced exactly as in Step 6.2, sharing the single twin already computed. The decks are generated together:

```bash
cd $AII/NS_fluids_lib/LSA_Rayleigh_Fortran
./bin/make_inputs sweep-eps --kr0 0.7 --eps 5e-3,1e-2,2e-2,4e-2 \
    --cells-per-r0 8 --n-periods 1.2 --snapshots 60 \
    --outdir runs/decks_eps
```

The verdict is then obtained from the assembled set:

```bash
./bin/run_lsa linearity --base runs/k0.7_base --kr0 0.7 \
    --case 5e-3=runs/k0.7_eps0.005 \
    --case 1e-2=runs/k0.7_eps0.01 \
    --case 2e-2=runs/k0.7_eps0.02 \
    --case 4e-2=runs/k0.7_eps0.04 \
    --t-min 4.0
```

The expected verdict is `LINEAR`, with a scaled-norm spread below 2 per cent and fitted rates agreeing to within 2 per cent of their mean. The tool returns a non-zero exit status when the sweep is unusable, so that it may gate a scripted workflow. The two failure modes call for opposite corrections and are distinguished by the time at which the collapse degrades: degradation at late time indicates the neglected second-order term and calls for a smaller amplitude, whereas degradation at early time indicates round-off and calls for a larger one, or for tighter solver tolerances.

## 8. The full sweep

Issued from this directory. Every calculation of the study is listed in `full_analysis.txt`, one per line, with comments describing what each group tests and what result to expect. Lines beginning with `#` are ignored, so a block may be commented out to run a subset without altering the field values.

```bash
cd $AII/NS_fluids_lib/LSA_Rayleigh_Fortran
./scripts/run_all.sh                        # full_analysis.txt, one rank
./scripts/run_all.sh full_analysis.txt 16   # sixteen ranks
./scripts/run_all.sh my_subset.txt 8        # a subset of the cases
```

The runner builds the analysis tools if they are absent, generates the deck for each case, advances it to completion by restarting from the most recent checkpoint as often as necessary, and extracts the interface mode amplitude to `runs/records/<label>.dat`. Progress is appended to `runs/run_all.log` and a summary table is written to `runs/summary.tsv` as each case finishes.

The sweep is intended to be started once and left to run, and it is safe to interrupt: a case whose record already exists is skipped, and a partly advanced case resumes from its checkpoint, so the script may simply be run again. A failure in one case is recorded and does not stop the others.

Cost scales roughly as the cube of the resolution, since the capillary time-step limit varies as $\Delta x^{3/2}$. The file should be ordered so that the cases most needed finish first, because the runner proceeds from top to bottom.

`further_analysis.md` sets out the same programme in narrative form, with the acceptance criterion for each stage and the interpretation of each failure mode.

## Run length and the transient

Every calculation is advanced for 5.0 e-folding times of the reference growth rate, not the 2.6 used earlier. The reason is worth stating, because a shorter run does not fail, it merely returns a wrong answer that looks right.

The imposed initial condition is a single Fourier mode of the interface, but the velocity field that accompanies it is not the matching eigenvector. The initial state is therefore a mixture of the growing mode and several decaying ones, and the amplitude grows at the eigenvalue only after the decaying part has died away. An exponential fitted before that point measures the mixture, and since a decaying contribution subtracts from the apparent rate, the fit comes out low. At 2.6 e-folding times roughly three per cent of transient remains, which was the whole of a discrepancy previously attributed to spatial discretisation.

The amplitude curve gives no warning of this: on a logarithmic axis it looks like a clean exponential long before it is one. The local growth rate does give warning, by still rising at the end of the record. That check is now performed automatically by `plotting/transient.py` for every measurement, and it prints the residual transient, the extrapolated asymptote and the run length required. A rate should not be quoted when it reports that the run is too short.


## 9. Figures

Issued from this directory. Figures follow the house style defined in `plotting/pub_style.py`: Computer Modern typography, inward ticks on all four sides, no axis padding, a data line width of 1.0 point, and PNG, PDF and SVG output at 600 dpi. Nothing in the figure scripts restates those settings, so a change to the house style takes effect everywhere at once.

Four figures describe a single case. `growth` shows the interface mode amplitude against time with the fitted exponential and both theoretical rates; `local_sigma` shows the local growth rate, whose plateau is what demonstrates that modal growth has been reached; `interface` shows the interface radius at several instants; and `field_<name>` shows a plotfile variable on the (r,z) plane with the interface laid over it.

The last two need more than an amplitude record. They need the interface profile r(z) at each instant and a plotfile variable on the plane, and neither can be recovered from a fitted growth rate. Reading them from the plotfiles requires yt and requires the plotfiles, which is impractical when a sweep leaves tens of gigabytes on a cluster carrying no matplotlib. The figure-data bundle exists for that reason: `dump_amplitude` and `dump_profiles` extract exactly what the four figures need into a few hundred kilobytes of plain text, once, on the machine holding the plotfiles.

```bash
cd $AII/NS_fluids_lib/LSA_Rayleigh_Fortran
./bin/dump_amplitude runs/k0.7_eps0.02 8.976 runs/bundles/k0.7/amplitude.dat
./bin/dump_profiles runs/k0.7_eps0.02 8.976 runs/bundles/k0.7 y_velocity

# all four figures, from the bundle alone; no yt, no plotfiles
python3 plotting/plot_results.py --bundle runs/bundles/k0.7 --kr0 0.7 --outdir figures
```

Both `scripts/run_nsmfp.sh` and `scripts/run_all.sh` write a bundle for every case automatically, immediately after the amplitude record and while the plotfiles are still present, so in normal use only the plotting command is issued. Set `FIELDS=` in the environment to skip the field slices, which are the larger part of a bundle, or `FIELDS="y_velocity pressure"` to export more than one.

Where the plotfiles are at hand, `--rundir` draws the same four figures directly, and the two paths have been checked to produce identical output:

```bash
python3 plotting/plot_results.py --rundir runs/k0.7_eps0.02 --wavelength 8.976 --kr0 0.7 --field y_velocity --outdir figures
```

A stored JSON record may also be passed with `--json`, which draws the two amplitude figures only, since a record carries no interface profile. The archived reference calculation is in that form.

The report figures are separate, and likewise read stored records rather than plotfiles:

```bash
python3 plotting/make_report_figures.py --datadir validation_data --archive validation_data/validation_kr0_0.7.json --outdir figures
```

`fig_sample/` holds the same set drawn from `validation_data/`, kept as a reference for what each figure should look like; see `fig_sample/README.md`.

## 9a. Growth rates and the dispersion figure from the records

The records written by the sweep are plain text, so the two tools below run on a login node or a workstation without the plotfiles and without yt.

```bash
python3 plotting/fit_records.py                    # every record, fitted and tabulated
python3 plotting/fit_records.py --t-min 4.0        # an explicit fit window
python3 plotting/plot_dispersion_sweep.py          # the measured dispersion figure
```

`fit_records.py` reports, for each case, the fitted growth rate, the coefficient of determination, the inviscid prediction and the departure from it, and flags a rate above the inviscid value, a poor fit, or a mean-radius drift beyond one per cent. Cases above $kr_0 = 1$ are reported as oscillatory with their predicted period rather than given a fitted rate, since the squared growth rate is negative there and a log-linear fit of an oscillation is meaningless. The unperturbed twin is likewise not fitted; its amplitude is reported as the noise floor against which every perturbed case is measured.

`plot_dispersion_sweep.py` draws the measured curve against the inviscid and viscous predictions, with a lower panel showing the departure from the viscous prediction as a percentage. That panel is where the character of the error becomes legible: a discretisation error appears as a roughly constant offset across the band, whereas a mistake in a deck or a fit window appears as scatter. Wavenumbers above the cutoff are excluded from the curve and the stable band is shaded instead. The sample `fig4_dispersion` carries the single measured point that `validation_data` supports; once a sweep exists this tool supersedes it.

## 9c. The extended figure set

Seven further scripts produce the figures that support the method rather than the single measurement. All follow the house style and none requires a new calculation beyond those already described.

```bash
# the modal extraction: spectrum, amplitudes, truncation and convergence
python3 plotting/plot_dmd.py --perturbed runs/k0.7_eps0.02 --base runs/k0.7_base --outdir figures
python3 plotting/plot_dmd.py --replot --outdir figures        # redraw without decomposing again

# the extracted eigenvector against the analytic eigenfunction
python3 plotting/plot_mode_shape.py --perturbed runs/k0.7_eps0.02 --base runs/k0.7_base --outdir figures

# growth rate over the wavenumber and Ohnesorge plane, surface and contour
python3 plotting/plot_ohnesorge.py --measured 0.7 0.02 --outdir figures

# break-up, from a run carried past the linear regime to pinch-off
python3 plotting/plot_pinchoff.py --bundle runs/bundles/pinch32 --kr0 0.7 --cells-per-r0 32 --outdir figures

# cost of the Jacobian-free path, with the machine recorded alongside
python3 plotting/benchmark_jvp.py --perturbed runs/k0.7_eps0.02 --base runs/k0.7_base --outdir figures

# the twin-run architecture, drawn rather than measured
python3 plotting/plot_schematic.py --outdir figures

# machine and run provenance on its own
python3 plotting/machine_info.py $AII
```

`plot_dmd.py` and `plot_mode_shape.py` read the plotfiles and therefore need the reader; both store the numbers behind their figures so that `--replot` redraws without repeating the decomposition. `plot_ohnesorge.py` and `plot_schematic.py` need no calculation at all. `plot_pinchoff.py` reads a bundle, but the run must be advanced until the neck closes, which the reference case is not: raise `--n-periods` until the minimum radius approaches zero, and the script reports plainly when the record stops short.

The break-up figures need a calculation the reference configuration does not provide, and two settings change. The run must reach pinch-off: the stop time is measured in e-folding times of the viscous rate, and from an amplitude of 2e-2 the neck closes after a growth by a factor of fifty, which with the startup transient puts break-up near $t = 14$ at the measured rate; `n_periods` of 6.0 gives a stop time near 18 and a margin of roughly a quarter for the nonlinear slowdown. And the mesh must resolve the neck: at 8 cells per radius a cell is an eighth of the radius, so the neck falls below one cell long before it closes and a break-up time measured there describes the mesh rather than the flow.

Both are set in the `pinch16` and `pinch32` cases of `full_analysis.txt`, which are run like any other case. There are two of them deliberately. Break-up is a singularity, no mesh resolves it to the end, and the break-up time therefore depends on the mesh; quoting one value without the other invites the first question a committee will ask. Run `pinch16` first, at roughly an eighth of the cost, and confirm that it breaks before committing to `pinch32`.

The threshold at which the neck counts as closed defaults to two cell widths, taken from `--cells-per-r0`, since a neck thinner than a cell is not represented by the interface reconstruction. It is printed with every result so that the break-up time and the criterion that produced it are read together. Where the record stops short of break-up the script says so and reports how far the neck actually fell, rather than returning a number.

`benchmark_jvp.py` writes `jvp_cost.json` beside its figure, holding every timing together with the processor, the memory, the compiler, the scheduler allocation where one exists, the solver commit, and both the processors present and the processors the process may actually use. That last distinction matters on a cluster, where a node advertising many cores may grant a handful inside an allocation, and a parallel efficiency computed against the wrong denominator is wrong by the ratio between them.

## 9b. The written report

`../analysis_report.typ` is the Typst source of the analysis report, with its bibliography in `../references.bib`; both sit one level above this directory because the report discusses one set of results, obtained by two implementations of the same method, and duplicating it in each tree would allow the two copies to drift. Compile it with

```bash
pip install typst
python3 -c "import typst; typst.compile('../analysis_report.typ', output='../analysis_report.pdf')"
```

or with the standalone binary, `typst compile ../analysis_report.typ`. Page geometry follows the thesis convention of 1.7 cm left and right and 1.9 cm top and bottom at 12 pt, mathematics is left in New Computer Modern Math, and the bibliography is rendered in APA style. Every peer-reviewed entry in `references.bib` carries a DOI and the two entries without one, a monograph and a laboratory report, carry a URL, so that each reference resolves independently of the document.

Figures cited in the report are those in `fig_sample/`, drawn by the current code from `validation_data/` and kept as a reference for what each should look like. Figures generated from a sweep are written to `figures/`; see `fig_sample/README.md`.

## 10. What follows

Steps 6 to 8 measure the growth rate through the interface-amplitude observable, which validates the physics but bypasses the NS-MFP machinery proper. Stage 5 of `further_analysis.md` exercises the decomposition and its mode shapes, Stage 6 compares against the Krylov-subspace stability module already present in the solver, and Stage 7 closes the loop against the nonlinear pinch-off calculation.

## 11. Diagnosis of common failures

Every entry below has been observed in practice rather than anticipated.

| Symptom | Cause | Remedy |
|---|---|---|
| `expecting num_FSI_outer_sweeps>=2` | the working tree is not at `0b8df42`; the defect is present at `cbcee7c` and again at `b014c1f` | `git checkout 0b8df42`, per Step 1 |
| `Permission denied` on a script in `scripts/` | extraction did not preserve the executable bit | `chmod +x scripts/*.sh`, or invoke as `bash scripts/run_nsmfp.sh` |
| `N singletons will be started`, and every rank reports `NProcs()= 1` | the launcher cannot reach a process manager, so N independent copies start instead of one N-rank job | run on one rank, or resolve the mismatch between `mpirun` and the libraries the solver was linked against; do not trust output produced this way |
| a calculation appears not to finish | the resolution defaults higher than intended, or the build is a debug build, or the ranks are singletons | read the grid and step count the deck generator reports, and see Section 12 |
| `make` reports `amr2d.gnu.DEBUG.ex is up to date` | `DEBUG = TRUE` in `GNUmakefile`, so the MPI executable was neither built nor refreshed | set `DEBUG = FALSE` and `USE_MPI = TRUE` and rebuild, or pass the debug executable explicitly |
| undefined reference to `amrex::FillPatchTower` | interrupted build left a zero-length object file | `find tmp_build_dir -name '*.o' -size 0 -delete` and rebuild |
| `ModuleNotFoundError: No module named 'matplotlib'` | the figure scripts are the one part of this tree that needs Python | `pip install numpy scipy matplotlib`; the analysis itself does not require them, and the records may be plotted elsewhere |
| `no nddataPLT* directories under the perturbed run` | the run never produced plotfiles, usually because the solver aborted at startup | read `runs/<case>/solver.err`; the FSI defect above is the most common cause |
| `Runtime environment uses unsupported PMI version PMIx` | `mpirun` belongs to a different MPI than the solver was linked against | match the launcher to the library `ldd` reports, per Section 15 |
| a restart repeats from the beginning | the checkpoint interval exceeds the steps completed per allocation | reduce `amr.check_int` |
| every mode amplitude returned as zero | the perturbation field vanishes at the first snapshot | increase `--skip`, per Step 6.3 |
| fitted rate far below expectation | the startup transient is included in the fit | set `--t-min` beyond the transient, near $t = 4$ |
| fitted rate above the inviscid value | an error in the deck | check `ns.yblob` against the axial domain height, and `ns.tension` and `ns.denconst` |
| mean radius drifts by more than 1 per cent | failure of mass conservation | the growth rate is invalid irrespective of the quality of the fit |

## 12. Cost, and what to do when a calculation is slow

The reference case at 8 cells per column radius is 1800 steps on a 32 by 72 grid and completes in minutes on a few ranks. Anything markedly slower than that has a cause worth identifying before more time is spent, because each of the causes below is a factor of several and they multiply.

Resolution dominates everything else. Work scales roughly as the cube of the cells per radius, since the cell count varies as the inverse square of the mesh spacing while the capillary stability limit on the time step varies as the three-halves power of it. The driver script of Step 6 takes the resolution as its fourth argument and defaults to 8, the reference value. Passing 32 instead gives a 128 by 288 grid and 33110 steps, which is approximately three hundred times the work of the reference case for a single wavenumber, and is the single most likely explanation for a calculation that appears not to finish. The deck generator reports the grid and the step count before anything is advanced, so the cost is knowable in advance and should be read:

```
grid 32x72  dr= 1.250E-01 dz= 1.247E-01
dt= 4.393E-03 (cap limit  1.757E-02)  steps=1800  snapshots=120
```

A debug build is the next largest factor. If `make` in `NS_fluids_lib` reports `amr2d.gnu.DEBUG.ex`, then `DEBUG = TRUE` is set, bounds checking and assertions are compiled in, and the solver runs several times slower than an optimised build. Step 3 covers how to confirm which executable is in use.

A launcher that starts singletons rather than one job of the requested size gives no speed-up at all, however many ranks are requested, while appearing to use them. Step 4 covers the test for this.

Two properties of the method reduce the total cost and are easy to overlook. The growth rate is obtained from the perturbed calculation alone through the interface amplitude, so a study confined to growth rates may omit every unperturbed twin and halve its cost; the twin is needed only for the velocity subspace and for the norm form of the proportionality test. And within a sweep over amplitudes at fixed wavenumber, a single twin is shared by every amplitude, which is why Step 7 generates those decks together.

The sweep of Step 8 is designed to be interrupted. A case whose record already exists is skipped and a partly advanced case resumes from its checkpoint, so a long sweep may be stopped and restarted freely, and `full_analysis.txt` should be ordered so that the cases most needed appear first.

## 13. Output format

All four reports are rendered as aligned tables with a rule above and below, so that a growth rate, its residual and its departure from theory may be read without counting columns, and so that a result pasted into a thesis appendix or a batch log keeps its alignment. The Fortran implementation renders them through `src/table_mod.f90` and the Python implementation through `python/report.py`, in the same shapes, so that the two may be compared line by line.

The Python renderer uses `tabulate` when it is installed and an equivalent renderer of its own when it is not. Neither the structure of the output nor the ability to run the analysis depends on which path is taken, because a missing formatting dependency is not a reason for a calculation of several hours to fail at the point of printing its answer. Both renderers emit plain ASCII rather than box-drawing characters, since output is read through terminals, log files and editors whose locale settings are not known in advance.

Formatted output is written only to the terminal. Records written to disk, in `runs/records/*.dat` and the JSON files written by `--save`, remain in their plain machine-readable form and are unaffected by any of this, so downstream scripts and the plotting code need not parse a table.

## 14. Running on a cluster

This document describes the workstation case. On a scheduled cluster almost every step differs: the solver is built against modules that the job must then load, the launcher is `srun` rather than `mpirun`, the sweep is submitted rather than started, and the cases may be advanced concurrently as a job array instead of one after another. `Running_the_Analysis_HPC_Fortran.md` covers all of that for the FSU Research Computing Center HPC and assumes the present document has been read.

Nothing here needs undoing before going to the cluster. The build, the decks, the observable and the interpretation of the result are the same in both places; what changes is how a calculation is started and what may interrupt it.

## 15. When mpirun and the solver come from different MPI installations

A machine with more than one MPI installed will usually have one `mpirun` first on the PATH and the solver linked against another. Nothing detects this at build time, and the symptom appears only at launch. The clearest form is

```
Runtime environment uses unsupported PMI version PMIx. Aborting.
```

which means the executable was built against MPICH while the `mpirun` that started it is Open MPI, or the reverse; the two use different process-management protocols and neither will speak the other's. The singleton failure of Section 11, in which each rank reports one rank, is the quieter form of the same disorder.

Establish which is which before changing anything:

```bash
which -a mpirun mpiexec
mpirun --version
ldd $AII/NS_fluids_lib/amr2d.gnu.MPI.ex | grep -iE 'mpi|mpich|pmix'
```

The third command names the library the solver actually loads. If it lists `libmpich`, the launcher must be MPICH's `mpiexec`; if it lists `libmpi.so` from Open MPI, the launcher must be Open MPI's `mpirun`. Whichever it is, the fix is to use the launcher that belongs to it, by absolute path if the PATH cannot be relied upon:

```bash
/usr/lib64/mpich/bin/mpiexec -n 4 $AII/NS_fluids_lib/amr2d.gnu.MPI.ex inputs
```

The alternative is to rebuild the solver against whichever MPI is first on the PATH, which is cleaner where the machine is yours to arrange. On Debian and Ubuntu, `update-alternatives --config mpi` and `--config mpirun` set the default pair consistently; where environment modules are in use, loading the matching module does the same. Either way the rule is that the launcher and the executable must come from one installation.

A single rank is not a way around this. `mpirun -n 1` still initialises MPI and still fails, as the error above was produced with one rank. To test the solver without a launcher at all, run the executable directly:

```bash
cd $AII/NS_fluids_lib/run2d
$AII/NS_fluids_lib/amr2d.gnu.MPI.ex inputs.growthrate.LSA
```

An MPI program started this way runs as a single process and, for most implementations, initialises without a process manager. That is enough to confirm the solver itself works and to separate a launcher problem from a solver problem, and it is enough for the smoke test, but it is not a way to run the study: the driver scripts take a rank count, and the reference case wants four.
