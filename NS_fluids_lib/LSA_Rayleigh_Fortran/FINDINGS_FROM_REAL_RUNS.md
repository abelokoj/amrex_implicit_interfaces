# Observations from the solver calculations

Every observation recorded here derives from building and running `amrex_implicit_interfaces`, rather than from inspection of the source or from synthetic testing. Several could not have been established by any other means. Each is recorded with the symptom by which it is recognised, so that it can be identified again rather than rediscovered.

**Status.** The defect described in Section 2 was corrected upstream at commit `0b8df42`, but it is present again at commit `b014c1f`, where it was observed on 3 September 2026 to abort both a generated deck and the distributed deck `run2d/inputs.growthrate.LSA`. The assertion must therefore be treated as live rather than closed. The remedy adopted is to work from `0b8df42`, which has been confirmed to run correctly; no local modification of the solver is then required, and the solver version is stated by a single commit identifier rather than by a description of an edit. The reference calculation was subsequently repeated on the corrected solver and agrees with the earlier result to 0.006 per cent, as reported in Section 6.2 of `analysis_report.typ`. Sections 1 and 2 are retained because results obtained during the affected interval must be reported as such, and because the build guidance remains current.

---

## 1. Building the solver

AMReX is vendored within the repository at `../amrex-master`, and is tracked at the same commit as the solver, so no separate clone is required despite the suggestion in the `GNUmakefile` header. The build additionally requires MPI, obtained on Debian and Ubuntu systems through `apt-get install libopenmpi-dev openmpi-bin`. From `NS_fluids_lib`, `make -j8` produces `amr2d.gnu.MPI.ex`, of approximately 109 MB and 187 object files, in roughly ten minutes on eight cores.

One failure mode deserves note because its symptom is misleading. If a build is interrupted, `make` may leave an object file of zero length, which it subsequently treats as current. The consequence appears much later, at the link stage, as an undefined reference to a symbol whose source file had apparently compiled successfully; in the present case the symbol was `amrex::FillPatchTower`. The remedy is to remove such files and rebuild:

```bash
find tmp_build_dir -name '*.o' -size 0 -delete && make -j8
```

## 2. A defect preventing all fluid-only calculations

At commit `cbcee7c`, every calculation lacking an ice, rigid or elastic material terminated during initialisation with the message `expecting num_FSI_outer_sweeps>=2`.

The defect was not specific to the decks used here. The decks distributed with the solver, `run2d/inputs.growthrate.LSA` and `run2d/inputs.capillary`, failed in precisely the same manner. Running a distributed deck is the most direct means of establishing whether a failure originates upstream, and is worth attempting before examining a locally written input file.

The quantity `num_FSI_outer_sweeps` is derived rather than supplied by the user. In `NavierStokes.cpp`, near line 3503, it is initialised to unity and incremented once for each material that is ice, FSI-rigid or FSI-elastic. For fluid-only problems it is therefore necessarily unity, and no input parameter can raise it. A requirement that it be at least two consequently rejects the entire class of fluid-only problems by construction.

Comparison of the two commits shows that six assertions requiring a value of at least two were introduced with the elastic-material work, and that the same change deleted a comment documenting the unit case. Of the three branches carrying the assertion, two lie in branches implying a rigid or elastic primary material and are sound; the third lies in the fall-through branch that fluid-only problems always take. A sibling routine elsewhere in the same file continued to handle the unit case explicitly, so the source was internally inconsistent.

The defect was corrected upstream at commit `0b8df42`, where the assertion is guarded by a test on the primary material, with an explicit no-check branch provided for the fluid-only case. That correction was verified by building from the corrected source without modification and confirming that both a distributed deck and a generated deck advance normally.

It is present again at commit `b014c1f`. The solver reports `num_FSI_outer_sweeps: 1` and `NFSI_LIMIT: 2` during initialisation and then aborts with `expecting num_FSI_outer_sweeps>=2`, exactly as at `cbcee7c`, and it does so for the distributed deck as well as for a generated one, which again establishes that the fault is not specific to this configuration.

The remedy is to check out `0b8df42`, as Step 1 of `Running_the_Analysis_Fortran.md` directs. Should a later head be required for some other reason, the assertion may be located and guarded directly:

```bash
cd $AII/NS_fluids_lib
grep -n "num_FSI_outer_sweeps>=2" *.cpp
git log --oneline -S "num_FSI_outer_sweeps>=2" -- NavierStokes.cpp | head
```

The second command lists every commit that added or removed the string, which distinguishes a revert from a fresh occurrence. Of the assertions it reports, the one to correct is the one in the fall-through branch that a fluid-only problem reaches, not those in branches that presuppose a rigid or elastic primary material. Guarding it on the primary material, in the manner of `0b8df42`, restores fluid-only operation; disabling the check outright also runs but removes a test that is sound for the elastic case, so the guard is preferable.

Any result obtained from a locally corrected build must be reported as such, and the equivalence across solver versions in Section 6.2 of `analysis_report.typ` is what licenses comparing it with results from an uncorrected one.

## 3. Naming of the output directories

The solver does not honour `amr.plot_file`. With `ns.visual_nddata_format = 1` it writes directories named `nddataPLT00000023`, together with `MOF_PLT*` directories and Tecplot files, rather than the `plt00023` that the input parameter would suggest. A tool expecting the documented name locates no output at all, and reports an empty run rather than a naming mismatch, which is why the symptom is worth stating explicitly.

The prefix `nddataPLT` is now the default in both, and may be overridden where a build behaves differently. The plotfiles themselves are in the standard AMReX format and carry 46 components, among them `x_velocity`, `y_velocity`, `PRES_MG`, the volume fractions `F01` and `F02`, and the level sets `L0101` and `L0202`.

## 4. Correctness of the plotfile readers

Applied to genuine solver output, the Fortran parser and the Python path through `yt` agree to the last digit recorded. For the calculation examined, both return a grid of 32 by 72 cells matching the deck, a time at step 23 of 0.101034047364, consistent with 23 steps of 0.0043928, a domain of 4.0 by 8.9760, and a maximum absolute field value of 5.977060209964134e-15.

That final quantity is physically informative rather than merely a check on parsing. With the perturbation amplitude set to zero, the base state remains at machine-precision velocity, which confirms that the quiescent column is an exact discrete equilibrium. This is the premise underlying the replacement of the constraining body force by a twin calculation, described in Section 2.2 of `analysis_report.typ`.

## 5. Choice of observable

This is the most consequential of the observations recorded here, and the one least likely to be anticipated from the source alone.

### 5.1 Velocity fields vanish identically at the initial instant

The perturbation is imposed on the interface, so both calculations of a twin pair begin from rest and the first perturbation snapshot is identically zero. The amplitude ranking of the decomposition projects each mode onto that first snapshot, so every amplitude was returned as zero and the ordering became meaningless. The growth rate reported under these conditions was negative and of order ten, which bears no relation to the physics.

Velocity fields remain a valid state vector, but only once the startup transient has been excluded.

### 5.2 The raw level set dilutes the signal

The level set is a signed distance function, so the difference between the perturbed and unperturbed calculations is approximately $\varepsilon\sin(kz)$ throughout the domain, whereas only the neighbourhood of the interface participates in the instability. In one calculation the level-set norm remained constant to 0.6 per cent while the interface amplitude grew by 37 per cent.

A stationary perturbation norm therefore does not by itself demonstrate a stable interface. It may instead indicate that the observable is unsuitable.

### 5.3 Resolution

The measurements reported use the quantity the dispersion relation predicts, namely the amplitude of the $k$-Fourier component of the interface radius, obtained from the zero contour of the level set. For $r = r_0 + a(t)\cos kz$ linear theory gives $a(t) = a(0)e^{\sigma t}$, so a log-linear fit returns the growth rate directly, unaffected by either the vanishing initial condition or the dilution described above.

The quantity is used alongside the decomposition rather than in place of it. The decomposition supplies the spectrum and the mode shapes; the interface amplitude supplies a robust scalar against which those may be checked.

## 6. Scaling of the perturbation amplitude

An amplitude of $10^{-3}$ on a grid of spacing 0.125 amounts to 0.8 per cent of a cell width, which lies far below the scale of the interface reconstruction. The resulting calculation returned noise, which is the outcome the proportionality test is designed to identify as round-off dominance. An amplitude of $2\times10^{-2}$, or 16 per cent of a cell, behaved as expected.

The value $\varepsilon_0 = 10^{-5}$ adopted by Ranjan et al. does not transfer. Their amplitude is a fraction of a free-stream state variable, whereas the present one is an interface displacement that must exceed the reconstruction error. Amplitude should therefore be scaled against the mesh spacing rather than against the state.

## 7. Practical notes on running

Speed was approximately 2 s per step at 32 by 72 cells on a single core with tight tolerances. Relaxing `mac.mac_abs_tol` and `mac.visc_abs_tol` from $10^{-10}$ to $10^{-8}$ approximately halves this, and both remain far below the perturbation amplitude, which is the governing requirement.

Checkpoint restart operates as expected and is the practical means of advancing a calculation beyond a scheduler limit. AMReX accepts the relevant settings as command-line overrides, so the deck need not be edited:

```bash
mpirun -n 16 ./amr2d.gnu.MPI.ex inputs amr.check_int=200
mpirun -n 16 ./amr2d.gnu.MPI.ex inputs amr.check_int=200 amr.restart=chk00400
mpirun -n 16 ./amr2d.gnu.MPI.ex inputs amr.restart=chk00400 max_step=1600 stop_time=20.0
```

The checkpoint interval must be smaller than the number of steps achievable within one allocation. If it is not, no checkpoint is written and each restart repeats from the beginning, while the accumulating plotfiles give every appearance of progress.

Restarts write additional plotfiles when `amr.plotfile_on_restart` is set, so the number of output directories exceeds `max_step / plot_int`. This is harmless, since pairing is by step number, but the file count should not be taken as a measure of progress.

Sorting requires care, and the appropriate key differs by purpose. Twin calculations are paired by step number, since matching steps must be differenced. A time series must instead be ordered by time, because a restart may emit a plotfile whose step number does not follow the time ordering of its neighbours, so that a step-ordered series is not monotonic in time. This distinction was encountered in practice: the Fortran and Python amplitude series agreed on every fitted quantity, which are selected by time value, but disagreed on the final element of the record. The last path returned by a glob should not be assumed to be the latest.

## 8. What the amplitude record does not carry

The amplitude record of a case, three columns of time, mode amplitude and mean radius, is sufficient to fit a growth rate and to draw the growth curve and the local rate. It is not sufficient to draw the interface shape or a field slice, because it carries no interface profile r(z) and no plotfile variable on the plane. That distinction is easy to miss, since a record is described as holding the result of a case, and the result in question is the growth rate.

The consequence is practical. A sweep leaves its plotfiles on the machine that ran it, typically a cluster carrying no plotting stack, and a record copied to a workstation permits two of the four figures and not the other two. Recovering the other two afterwards means returning to the plotfiles, which is possible only for as long as they are retained.

Both implementations therefore export a figure-data bundle for every case, at the moment the amplitude record is written and while the plotfiles are still present. A bundle holds the interface profile at every snapshot and a field slice at a few instants, in plain text, and runs to a few hundred kilobytes against tens of gigabytes of plotfiles. `plot_results.py --bundle` draws all four figures from it with neither yt nor the plotfiles present.

The two paths have been checked against each other rather than assumed equivalent: drawing the four figures from a bundle and from the plotfiles of the same calculation gives byte-identical output for three of them, the fourth differing only in sub-pixel antialiasing, with the underlying records agreeing to 7.3e-16 relative.

One hazard emerged in making the format portable, and is recorded because it is invisible from within Fortran. Written with the default two-digit exponent descriptor, `es24.16`, a real whose exponent needs three digits is emitted with the exponent letter dropped, so 4.66e-310 appears as `4.6598606362008991-310`. Fortran reads that form back correctly, so a round trip within one language succeeds; no other language accepts it, and the reader that draws every figure rejects the file. The writers use `es25.16e3`, which makes the letter unconditional, and the reader additionally tolerates the two-digit form so that a bundle written by an earlier build remains readable.

## 9. Run length, and a bias that looked like discretisation error

A departure of about four per cent between measurement and theory was attributed to spatial discretisation until three candidate causes were tested against each other.

Replacing the quadratic viscous approximation with the exact viscous eigenvalue moves the reference from 0.328941 to 0.329367 and therefore widens the departure. Adding the ambient phase, at a density ratio of 1.225e-3, lowers the theoretical rate by 0.013 per cent, three orders of magnitude too little to matter. Neither is the cause.

The cause is the duration of the runs. The initial condition perturbs the interface but not the velocity field, so the state begins as a mixture of the growing mode and decaying ones, and the modal rate is reached only once the latter have gone. At 2.6 e-folding times roughly three per cent of transient remains, biasing every fitted rate low by that amount.

Three observations follow, and each was previously read as evidence of something else. The local growth rate was still rising at the end of every record, which the figure caption described as a plateau. The deficit was nearly constant across the wavenumber band, which was read as a discretisation signature, but follows directly from specifying run length in e-folding times, since every case is then truncated at the same point in its own development. And mesh refinement appeared to worsen agreement, because the 8- and 16-cell cases were both truncated at the same physical time and converged correctly on the answer for a truncated record.

The run length is now 5.0 e-folding times, at which the residual transient falls below one part in a thousand. The check is automatic: `transient.py` assesses every measurement and states whether a rate may be quoted.
