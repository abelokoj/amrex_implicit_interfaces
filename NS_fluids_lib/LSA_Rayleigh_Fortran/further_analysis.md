# Programme of further calculations

The calculations that remain, in order of dependency, with the expected outcome and an acceptance criterion for each. All build upon the protocol validated in Section 5 of `analysis_report.typ`.

Installation, the build of the solver and the first reference calculation are set out step by step in `Running_the_Analysis_Fortran.md`, which this document presumes has been followed. Every calculation described here is also listed in `full_analysis.txt`, which `scripts/run_all.sh` reads and executes in a single sweep. The stages below explain what each group of cases tests and how its output should be interpreted; the case file is the operational form of the same programme.

Commands assume the working directory `amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Fortran`, with all output written inside that directory: simulation output to `runs`, figures to `figures`. The Python implementation is invoked identically from its own directory, substituting `python3 python/<script>.py` for `./bin/<program>`.

## Established position

At $kr_0 = 0.7$, with 8 cells per column radius and $\varepsilon = 2\times10^{-2}$, the measured growth rate reaches a plateau at $\sigma = 0.3189$ with a coefficient of determination of 0.99997, lying 3.7 per cent below the viscous estimate of 0.3289. The response is proportional to the imposed amplitude to better than 0.11 per cent over a fourfold range, and the result is independent of the solver version to 0.006 per cent. The measured record is retained in `validation_data/`.

Three matters remain open: the independence of the plateau from the mesh, the behaviour at wavenumbers other than 0.7, and the extent to which the outer boundary influences the result.

## Reference values

Growth rates for $r_0 = \sigma_s = \rho_l = 1$ and $\mu_l = 0.02$, corresponding to an Ohnesorge number of 0.02, which are the defaults of the generated decks.

| $kr_0$ | $\lambda = 2\pi/k$ | inviscid | viscous estimate | $1/\sigma$ | run to $t \approx$ |
|---|---|---|---|---|---|
| 0.20 | 31.416 | 0.138220 | 0.137025 | 7.30 | 29 |
| 0.30 | 20.944 | 0.201236 | 0.198555 | 5.04 | 20 |
| 0.40 | 15.708 | 0.256692 | 0.251937 | 3.97 | 16 |
| 0.50 | 12.566 | 0.301558 | 0.294151 | 3.40 | 14 |
| 0.60 | 10.472 | 0.332128 | 0.321503 | 3.11 | 12 |
| 0.697 | 9.015 | 0.343339 | 0.329074 | 3.04 | 12 |
| 0.70 | 8.976 | 0.343327 | 0.328941 | 3.04 | 12 |
| 0.80 | 7.854 | 0.326909 | 0.308273 | 3.24 | 13 |
| 0.90 | 6.981 | 0.264730 | 0.241543 | 4.14 | 17 |
| 0.95 | 6.614 | 0.199198 | 0.173955 | 5.75 | 23 |

These may be regenerated from `dispersion_mod` or from `python/dispersion.py`. The required run length scales as the reciprocal of the growth rate, so cases near the ends of the unstable band must be advanced considerably further in time than those near the maximum.

Measured values are expected to fall a few per cent below the viscous column. The viscous estimate is a quadratic approximation rather than the exact viscous eigenvalue, and finite resolution under-resolves the curvature. A value exceeding the inviscid column is physically impossible and indicates an error.

## Stage 0. Building

```bash
git clone https://github.com/msussman42/amrex_implicit_interfaces.git
cd amrex_implicit_interfaces/NS_fluids_lib
make -j16
```

Should the clone fail on the file named `"u_vector".xml`, whose name contains quotation marks, the option `git -c core.protectNTFS=false clone` succeeds on a Linux filesystem.

The expected product is `amr2d.gnu.MPI.ex` of approximately 109 MB. AMReX is vendored at `../amrex-master` and requires no separate clone. Should the link fail with an undefined reference to `amrex::FillPatchTower`, remove any zero-length object files and rebuild, as described in Section 1 of `FINDINGS_FROM_REAL_RUNS.md`.

Confirm that a distributed deck runs before proceeding:

```bash
cd run2d && mpirun -n 1 ../amr2d.gnu.MPI.ex inputs.growthrate.LSA
```

Time-stepping should begin. Termination with `expecting num_FSI_outer_sweeps>=2` indicates a commit carrying the fluid-only assertion defect. The defect was corrected upstream at `0b8df42` but is present again at `b014c1f`, so updating the working tree is not by itself a remedy; Section 2 of `FINDINGS_FROM_REAL_RUNS.md` gives the commands that locate the assertion at the commit in use and states which of them to guard.

Then build and test the analysis tools:

```bash
cd ../LSA_Rayleigh_Fortran && make && make test    # 138 assertions
```

## Stage 1. Independence of the plateau from the mesh

This is the first priority, since no growth rate should be quoted until it is established. The comparison reported in Section 6.3 of `analysis_report.typ` covers only the startup transient and does not bear on the plateau value.

```bash
for N in 8 16 32; do
  ./bin/make_inputs single --kr0 0.7 --eps 2e-2 \
      --cells-per-r0 $N --rmax 4.0 --safety 0.9 \
      --n-periods 5.0 --snapshots 120 --outdir runs/decks_res$N
done
```

A run length of 2.6 e-folding times reaches $t \approx 8$, comfortably beyond the transient and before the amplitude reaches $a/r_0 \approx 0.15$, at which nonlinearity intervenes.

Only the perturbed deck of each pair need be advanced. The interface-amplitude measurement reads $r(z,t)$ directly and requires no twin, which halves the cost; the twin is needed only for the velocity subspace of Stage 4.

```bash
for N in 8 16 32; do
  mkdir -p runs/res$N && cd runs/res$N
  cp ../decks_res$N/inputs.k0.7.eps0.02 inputs
  mpirun -n 16 ../../amr2d.gnu.MPI.ex inputs \
      mac.mac_abs_tol=1.0e-8 mac.visc_abs_tol=1.0e-8 amr.check_int=200
  cd ../../LSA_Rayleigh_Fortran
done

for N in 8 16 32; do
  ./bin/dump_amplitude runs/res$N 8.976 runs/records/res$N.dat
done
```

Expected outcome:

| cells per $r_0$ | grid | $\sigma$ for $t > 5$ | $r^2$ |
|---|---|---|---|
| 8 | 64 by 144 | approximately 0.319 | above 0.999 |
| 16 | 128 by 288 | 0.321 to 0.328 | above 0.999 |
| 32 | 256 by 576 | 0.325 to 0.330 | above 0.999 |

The criterion is a monotone approach to approximately 0.329 with a change below 2 per cent between the two finest meshes, after which the finest value may be quoted and the sequence reported as the evidence of convergence.

Three failure modes should be distinguished. A value exceeding 0.343 is physically impossible and indicates an error in the deck, most often in `ns.tension`, `ns.denconst`, or a disagreement between `ns.yblob` and the axial domain height. A value falling with refinement, or failing to reach a plateau, suggests confinement by the outer boundary, which Stage 4 addresses. A drift in mean radius exceeding roughly 1 per cent indicates a failure of mass conservation, which invalidates the growth rate irrespective of the quality of the fit.

Cost scales approximately as the cube of the resolution, since the capillary time-step limit varies as $\Delta x^{3/2}$. The 32-cell case should be attempted last, and only if the coarser pair has not already converged.

## Stage 2. Proportionality of the response

The proportionality test reported in Section 6.1 of `analysis_report.typ` was conducted at 8 cells per radius. It should be repeated at whichever resolution Stage 1 identifies as converged, since the round-off floor scales with the mesh spacing and the admissible range of amplitude moves with it.

```bash
./bin/make_inputs sweep-eps --kr0 0.7 --eps 5e-3,1e-2,2e-2,4e-2 \
    --cells-per-r0 16 --n-periods 5.0 --snapshots 120 \
    --outdir runs/decks_eps
```

The unperturbed twin is required here, since the test in its $\lVert Q'\rVert$ form compares differenced fields. A single twin serves every amplitude. Solver tolerances should be tightened to $10^{-10}$, so that the residual remains well below the smallest amplitude in the sweep.

```bash
./bin/run_lsa linearity --base runs/lin_base --kr0 0.7 \
    --case 5e-3=runs/lin_eps0.005 \
    --case 1e-2=runs/lin_eps0.01  \
    --case 2e-2=runs/lin_eps0.02  \
    --case 4e-2=runs/lin_eps0.04  \
    --t-min 4.0
```

The argument `--t-min` excludes the transient, without which the fitted rates are contaminated and the verdict is meaningless.

The expected verdict is `LINEAR`, with a scaled-norm spread below 2 per cent and fitted rates agreeing to within 2 per cent of their mean. The tool returns a non-zero exit status when the sweep is unusable, so that it may gate a scripted workflow.

The two failure modes call for opposite corrections and are distinguished by the time at which the collapse degrades. Degradation at late time indicates the neglected second-order term and calls for a smaller amplitude; degradation at early time indicates round-off and calls for a larger one, or for tighter solver tolerances.

## Stage 3. The dispersion relation

This is the principal result and presupposes both preceding stages.

```bash
./bin/make_inputs sweep-k --kr0 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
    --eps 2e-2 --cells-per-r0 16 --n-periods 5.0 --snapshots 120 \
    --outdir runs/decks_k
```

Each wavenumber constitutes a separate discrete problem, since the box length is one wavelength and therefore varies with $k$. The perturbed decks alone are required.

The fit window must scale with the growth rate rather than being fixed, since the transient persists for approximately $1.3/\sigma$. Taking the final 40 per cent of each record achieves this without further intervention.

The expected outcome is a set of points lying 3 to 8 per cent below the viscous curve across the unstable band, with the largest departures at the extremes, where the required box is long and the growth slow. The principal check is that the measured maximum falls at $kr_0 \approx 0.7$, against the theoretical 0.697.

Acceptance requires a monotone rise to a maximum near 0.7, a monotone fall thereafter, every point below the inviscid curve, and coefficients of determination above 0.99. A point above the inviscid curve, or a maximum displaced by more than one sampled wavenumber, indicates an error; the disagreement of `ns.yblob` with the axial domain height is the most common cause, since the imposed sinusoid is then not periodic within the box.

A calculation at $kr_0 = 1.1$ serves as a negative control. That wavenumber lies beyond the cut-off, so the amplitude must decay. Recovering a negative growth rate there is evidence that the procedure does not manufacture growth from discretisation error, and is worth reporting alongside the curve.

## Stage 4. Influence of the outer boundary

The outer boundary is a no-slip wall at $r = 4r_0$, which confines the surrounding gas. Repeating the reference case at $r_{\max} = 6r_0$ tests whether that confinement affects the measured rate.

```bash
./bin/make_inputs single --kr0 0.7 --eps 2e-2 --cells-per-r0 16 \
    --rmax 6.0 --n-periods 5.0 --snapshots 120 --outdir runs/decks_dom
```

Agreement with the corresponding Stage 1 case to within the mesh convergence error establishes that the domain is adequate. A material difference would require the entire study to be repeated at the larger radius, which is why this case is inexpensive insurance rather than an optional extra.

## Stage 5. The subspace and the decomposition

Stages 1 to 4 employ the interface-amplitude scalar, which validates the physics but bypasses the NS-MFP machinery. This stage exercises the method proper and requires twin pairs.

```bash
./bin/run_lsa analyse \
    --perturbed runs/lin_eps0.02 --base runs/lin_base \
    --kr0 0.7 --skip 60 --t-min 4.0 --save runs/dmd_k0.7.json
```

The argument `--skip` is essential rather than optional. Since the perturbation is interfacial, both calculations begin from rest and the perturbation field vanishes identically at the initial instant, which renders the amplitude ranking degenerate; Section 5.1 of `FINDINGS_FROM_REAL_RUNS.md` records the consequence.

Expected: a leading stationary mode whose growth rate matches Stage 1 to within a few per cent, agreement between the decomposition and the norm fit reported as consistent, and strongly damped remaining modes. Amplitudes returned as zero indicate that `--skip` is too small. A leading mode that is oscillatory indicates that the subspace is dominated by acoustics or by startup transients, since the capillary mode is non-oscillatory.

Convergence with respect to the subspace should then be established, following Section 5 of Ranjan et al.:

```bash
./bin/run_lsa converge --perturbed runs/lin_eps0.02 \
    --base runs/lin_base --kr0 0.7 --skip 60
```

The growth rate is expected to be stable to within 1 per cent across subspace sizes, since a single mode dominates.

Comparison of the mode shapes is a stronger test than the growth rate alone and is available at no additional computational cost. The eigenvectors are held in `dmd_result_t%modes`; reshaping the leading one onto the grid and comparing it with the analytic interface eigenfunction tests the spatial structure directly.

## Stage 6. Comparison with the existing stability module

The solver contains a Krylov-subspace stability module of its own, whose comments cite the same decomposition literature. A direct comparison on the same configuration quantifies the difference and forestalls the obvious question.

```bash
mkdir -p runs/builtin && cd runs/builtin
cp ../../run2d/inputs.growthrate.LSA inputs
mpirun -n 16 ../../amr2d.gnu.MPI.ex inputs
```

Comparable growth rates are expected. The argument advanced in Section 1.4 of the README concerns structure rather than accuracy, namely one forward march per case against restart-and-repeat iteration, and drift removed by a computed twin against drift handled within the restart cycle. This calculation renders that argument quantitative in cost and in the number of parameters that must be supplied in advance.

## Stage 7. Nonlinear validation

```bash
mkdir -p runs/nonlinear && cd runs/nonlinear
cp ../../run2d/inputs.growthrate.LSA.Design_of_Experiments inputs
mpirun -n 16 ../../amr2d.gnu.MPI.ex inputs
```

This is a finite-amplitude single-mode calculation at the critical wavelength. The amplitude is expected to follow $e^{\sigma t}$ with the measured rate while $a/r_0 \lesssim 0.15$, then to depart and proceed to pinch-off. Extrapolating the linear phase should predict the time of break-up to within roughly 20 per cent, which closes the loop between the linear analysis and the nonlinear behaviour it is intended to describe.

## Priority

Stage 1 is prerequisite to quoting any growth rate, and Stage 2 to claiming that the operator decomposed is the linearised one. These are independent of one another and may be run concurrently. Stage 3 depends on the resolution and amplitude they determine. Stages 4 to 7 strengthen the case without blocking it.

Independently of the calculations, the defect described in Section 2 of `FINDINGS_FROM_REAL_RUNS.md` should be reported as having affected results obtained during the interval concerned, together with the equivalence established across solver versions. The defect is present at commit `b014c1f`, so a local correction is required there; it must be recorded alongside any result obtained from such a build.
