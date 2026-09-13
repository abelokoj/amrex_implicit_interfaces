# Running the analysis on the FSU RCC HPC

This document covers what is specific to the Florida State University Research Computing Center cluster: how the work is submitted, how it is made to run in parallel, and which failures are properties of a scheduled machine rather than of the analysis. `Running_the_Analysis_Fortran.md` describes the method itself, what each command does and how to read the result, and is assumed here rather than repeated. Nothing in it needs undoing to work on the cluster; the build, the decks, the observable and the interpretation of a growth rate are the same in both places. What changes is how a calculation is started and what may interrupt it.

The short version, for a reader who wants the sequence and not the reasoning.

Connect, then set the environment. This must be repeated in every new shell, since neither the modules nor the variable persist between sessions, and it must be done before building or submitting anything:

```bash
ssh $USER@hpc-login.rcc.fsu.edu

module purge
module load gnu/13
module load openmpi
export AII=$HOME/amrex_implicit_interfaces
```

The exact module names vary as the cluster is updated, so confirm them with `module avail gnu` and `module avail openmpi` rather than assuming; `module load gnu openmpi` takes the current defaults. Whatever is loaded here must also be loaded inside every job, which is why the same lines appear in the submit scripts, and Section 3 explains what goes wrong when the two disagree.

Then, once:

```bash
cd $AII/NS_fluids_lib && git checkout 0b8df42 && make -j16
cd LSA_Rayleigh_Fortran && chmod +x scripts/*.sh && make && make test
```

And thereafter:

```bash
srun --pty -t 00:10:00 -n 4 -A quicktest bash                  # confirm it starts
scripts/submit_hpc_array.sh full_analysis.txt 8                # the sweep
scripts/merge_summaries.sh                                     # when it finishes
```

Adding `export AII=$HOME/amrex_implicit_interfaces` and the `module load` lines to `~/.bashrc` saves repeating them, at the cost of making the environment implicit; recording which modules a result was produced under matters more here than the keystrokes saved, so the explicit form is preferred while the study is in progress.

## 0. Working on a workstation and the cluster at once

Running part of the study in each place is reasonable and saves time, but three things have to hold or the results cannot be combined.

Both builds must be at commit `0b8df42`. This is the commit at which the fluid-only assertion defect is corrected; later commits, including the current head, reintroduce it. A study assembled from two different solver commits is not a study of one solver.

A twin pair must never be split across the two machines. The perturbed calculation and its unperturbed twin are differenced snapshot by snapshot, and the subtraction stands in for the NS-MFP body force only because both experience the same drift and therefore cancel it. Across two machines the compiler, its optimisation flags, the BLAS and LAPACK implementations, the MPI library and the domain decomposition all differ, and each of those reorders floating-point operations. The two base states then drift differently, the cancellation fails, and what emerges is not an obvious error but a plausible growth rate that is wrong by an unknown amount. Base-state velocities sit at $6\times10^{-15}$ in a matched pair, so the tolerance here is far tighter than in an ordinary simulation. Divide the work by case, not within a case: run $kr_0 = 0.7$ entirely in one place and $kr_0 = 0.5$ entirely in the other.

The machine and toolchain must be recorded with each result. Assembling a dispersion relation means comparing growth rates obtained in different places, so a systematic offset between two builds would appear as scatter with no physical cause. Guard against it once, cheaply: run the reference case at $kr_0 = 0.7$ on both machines and compare the fitted rates. Agreement to a few hundredths of a per cent is what has been observed across solver versions and settles the question; a larger discrepancy should be understood before the sweep is trusted. Note also that for growth rates alone the twin may be omitted entirely, since the interface-amplitude observable reads the perturbed run by itself, which removes this concern along with half the compute.

## 1. The distinction that matters most

The cluster offers two kinds of job, and choosing wrongly costs days. Batch jobs are submitted to the scheduler, wait in queue, run without further input and leave their output behind. Interactive jobs reserve a compute node for a period so that work may be done directly at the command line.

The sweep is batch work. An interactive job ties the calculation to a live terminal, so a dropped connection, a closed laptop or an expired session ends it, and interactive allocations are short. Use an interactive job for exactly two things: the smoke test that confirms the solver starts, and inspecting output afterwards. Everything that takes hours goes through `sbatch`.

Jobs do not start immediately. They queue, typically for a few minutes, and the wait grows with the resources requested.

## 2. Accounts and their time limits

The account determines how long a job may run, and a sweep submitted to the wrong one is killed partway through.

| Account | Limit | Suited to |
|---|---|---|
| `genacc_q` | 14 days | the default choice for this work |
| `condor` | 90 days | long serial work; not suitable, as it does not take MPI jobs |
| `backfill` | 4 hours | short jobs |
| `backfill2` | 4 hours | short jobs, and may be pre-empted |
| `quicktest` | 10 minutes | confirming that a submit script is correct; the smoke test of Section 4 |

An owner-based account, where the research group holds one, is preferable to any of these and carries priority. Use `quicktest` first, every time a submit script changes, because ten minutes spent finding a typo in a module name is better than discovering it after a day in queue.

The limits in this table are ceilings on what may be requested, not defaults that a longer request is trimmed to. Asking for more than the account allows is rejected at submission with `Requested time limit is invalid (missing or exceeds some limit)`, which refers to the time limit alone and implies nothing about the account or the resources requested alongside it.

Preemption on `backfill2` deserves a specific warning here. A pre-empted case is not merely slower; it stops wherever it happens to be, and unless a checkpoint has been written the restart repeats from the beginning. Section 6 covers the checkpoint interval that makes this survivable.

## 3. Building on the cluster

The solver is built once, from the pinned commit, on a login node:

```bash
cd $AII
git checkout 0b8df42
cd NS_fluids_lib
make -j16
ls -la amr2d.gnu.*.ex
module list
```

Three things are being checked. That the working tree is at `0b8df42`, since later commits reintroduce the assertion that aborts every fluid-only calculation. That the product is `amr2d.gnu.MPI.ex` and not `amr2d.gnu.DEBUG.ex`, since a debug build carries bounds checking and no optimisation and runs several times slower; if the debug name appears, set `DEBUG = FALSE` and `USE_MPI = TRUE` in `GNUmakefile` and rebuild cleanly. And what `module list` reports, because those same modules must be loaded inside every job.

That last point is the most common cause of a job that fails immediately after waiting hours in queue. A solver built against one MPI implementation and launched under another does not fail at link time; it fails at run time, on the compute node, after the allocation has been granted. The submit scripts load `gnu/13` and `openmpi`; if `module list` on the build node reports something else, edit them to match before submitting anything.

The analysis tools are built inside each job rather than on the login node, so nothing needs doing here for them.

## 4. Confirming that the solver runs, interactively

This is the one place an interactive job is the right instrument. Request a short allocation, run the distributed deck, and stop as soon as time-stepping begins.

```bash
srun --pty -t 00:10:00 -n 4 -A quicktest bash
cd $AII/NS_fluids_lib/run2d
srun -n 4 $AII/NS_fluids_lib/amr2d.gnu.MPI.ex inputs.growthrate.LSA
```

The time limit requested must not exceed the limit of the account, or Slurm refuses the allocation outright with `Requested time limit is invalid (missing or exceeds some limit)` rather than reducing the request to what the account permits. Ten minutes is the ceiling on `quicktest` and is ample here, since the calculation is to be interrupted as soon as it starts stepping. Where longer is wanted, use an account from the table in Section 2 that permits it:

```bash
srun --pty -t 00:30:00 -n 4 -A backfill bash     # up to 4 hours
srun --pty -t 00:30:00 -n 4 -A genacc_q bash     # up to 14 days
```

For a session this short `backfill` usually starts sooner than `genacc_q`, since brief jobs fill gaps in the schedule. Avoid `backfill2`, where the session may be pre-empted while in use.

Lines beginning `STEP` are the whole of what is wanted; interrupt once they appear. This deck is a smoke test and there is no reason to let it finish.

Termination with `expecting num_FSI_outer_sweeps>=2` means the working tree is not at `0b8df42`; return to Section 3.

Note that the launch inside the allocation uses `srun` rather than `mpirun`. That is not incidental, and Section 5 explains why.

## 5. Why every launch uses srun

The RCC documentation recommends `srun` over `mpirun` inside a submit script, because `srun` integrates with the scheduler. The practical consequence here is larger than a matter of preference.

A launcher that cannot reach the process manager does not stop. It starts independent single-rank copies of the program instead of one job of the requested size, announcing this only through a message about no PMIx server being reachable and a report from each process that it sees one rank. Every copy then advances the whole problem in the same directory and overwrites the plotfiles of the others. The calculation appears to succeed and the results are meaningless. This has already been observed with `mpirun` on a workstation and is a standing risk wherever a resource manager is present.

Under `srun` the rank count and the placement come from the allocation, so there is no second number to disagree with the first. Both runners here read a `LAUNCHER` variable and the cluster submit scripts set it to `srun`; the runners also refuse to start if they detect the singleton condition. Outside a job script the same override is available:

```bash
LAUNCHER=srun ./scripts/run_all.sh full_analysis.txt 16
```

## 6. Two settings that decide whether a job survives its time limit

The sweep is designed to be interrupted: a case whose record already exists is skipped and a partly advanced case resumes from its most recent checkpoint. Resubmitting after a job hits its wall-clock limit therefore continues rather than repeats. This holds only under two conditions.

The checkpoint interval must be small enough that a checkpoint is written within a single allocation. If it is not, no checkpoint is ever written, each restart begins again from the start, and accumulating plotfiles give every appearance of progress. The runner sets `amr.check_int` to 50 by default and it is adjustable through the `CHECK_INT` environment variable. On a queue with a four-hour limit, or on `backfill2` where the job may be pre-empted at any moment, it should be smaller.

The wall-clock limit requested must be an honest overestimate. Underestimating kills the job; overestimating by a wide margin lengthens the wait in queue. Better to overestimate than to underestimate.

## 7. The sequential route

One job, every case in order:

```bash
cd $AII/NS_fluids_lib/LSA_Rayleigh_Fortran
sbatch scripts/submit_hpc.sh                  # full_analysis.txt
sbatch scripts/submit_hpc.sh my_subset.txt    # a subset
squeue --me
tail -f runs/run_all.log
scancel JOBID
```

The rank count is read from `SLURM_NTASKS` rather than written a second time, so the value requested from the scheduler and the value passed to the solver cannot drift apart. Output goes to `runs/`, and the summary table to `runs/summary.tsv`.

This route is simple and suits a handful of cases. Its cost is the sum of the case times, which for the full file is long.

## 8. The array route, which is what makes it fast

The cases of the sweep are independent calculations, so they should be advanced at the same time rather than one after another. A Slurm job array does exactly this: one task per case, each with its own allocation.

```bash
cd $AII/NS_fluids_lib/LSA_Rayleigh_Fortran
scripts/submit_hpc_array.sh                        # every case, 8 ranks each
scripts/submit_hpc_array.sh full_analysis.txt 16   # 16 ranks each
scripts/submit_hpc_array.sh my_subset.txt 8
```

The wrapper counts the cases, since the array range must be known at submission, and then submits `1-N`. Where the sequential route takes the sum of the case times, the array takes the longest single case plus the queueing. For the sixteen cases of `full_analysis.txt` that is the difference between a week and a day, and the improvement costs nothing in compute: the same work is done, merely not in single file.

Three environment variables tune it:

```bash
THROTTLE=4 WALLTIME=12:00:00 ACCOUNT=genacc_q scripts/submit_hpc_array.sh
```

`THROTTLE` caps how many tasks run at once and defaults to 8. Submitting a hundred at once is discourteous on a shared machine and does not finish sooner, since the scheduler grants what is free either way. `WALLTIME` is per task, not for the whole array, which is why it may be much shorter than a sequential job would need. `ACCOUNT` selects from the table in Section 2.

Each task writes its own log and its own summary, so that concurrent tasks do not interleave lines in a shared file. Merge them when the array finishes:

```bash
scripts/merge_summaries.sh
```

This produces `runs/summary.tsv` in the same format the sequential runner writes, so the plotting and reporting code need not know which route produced the results.

A property worth stating explicitly, because it is the concern that usually arises with parallel submission. Each case in the file is generated, advanced and reduced entirely within one task. A perturbed calculation and its unperturbed twin are therefore never split across tasks, and never across machines or builds. The twin subtraction stands in for the NS-MFP body force and relies on both twins experiencing the same discretisation and the same drift, which would not survive being run under different compilers or different rank counts. The array respects that; splitting a pair by hand across two machines does not.

## 9. How many ranks a case can use

More ranks is not always faster, and at the reference resolution it stops helping early.

AMReX decomposes the domain into boxes and distributes those. The reference case is 32 by 72 cells, which at the usual `amr.max_grid_size` yields only a handful of boxes; ranks beyond the box count sit idle, and well before that the surface-to-volume ratio makes communication and the implicit pressure solve dominate. Four ranks is realistic there, eight shows diminishing returns, and sixteen is unlikely to beat eight.

The picture changes with resolution. At 16 and 32 cells per radius the domain is 64 by 144 and 128 by 288, there are many more boxes, and parallel efficiency is much better. Since those are also the expensive cases, this works in the right direction: give the coarse cases few ranks and the fine cases many. Where scaling disappoints at high resolution, reducing `amr.max_grid_size` to 16 or 32 produces more and smaller boxes to distribute.

This interacts with the array. Sixteen tasks at sixteen ranks each is 256 cores, which will queue for a long time. Eight ranks per task with a throttle of eight is 64 cores and will start far sooner, and for this problem will very likely finish first.

## 10. Requesting resources without waiting longer than necessary

The scheduler documentation is emphatic that fewer constraints start sooner, and the advice applies directly here.

Do not specify the number of nodes. Requesting a node count, or a large core count on one node, is a common reason for a job that sits pending indefinitely, since the scheduler must find that exact shape. Asking for a number of tasks and letting the scheduler place them across an arbitrary number of nodes is far more likely to start. The submit scripts here use `--ntasks` alone for this reason.

Do not request memory unless it is genuinely needed. The default allocation is about 3.9 GB per CPU, which is ample for these grids; the reference case ran on a single core in 3 GB. Where more is wanted, raising the task count is the recommended route in preference to `--mem` or `--mem-per-cpu`.

The RCC provides a submit script generator at `https://manage.rcc.fsu.edu/ssg`, which is a convenient check on the syntax of anything written by hand.

## 11. Watching and diagnosing a running sweep

```bash
squeue --me                                  # queued and running jobs
squeue --me -t PENDING                       # what is still waiting, and why
sstat -j JOBID                               # a running job
scancel JOBID                                # a whole job or array
scancel JOBID_3                              # one task of an array
sacct -j JOBID --format=JobID,State,Elapsed,MaxRSS,ExitCode
tail -f runs/run_all.log                     # the sequential route
tail -f runs/array/run_*_3.log               # task 3 of an array
```

`squeue --me` reports one line per job, of which the state column is the one to read:

```
JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
123456 genacc_q      trap   abc12a PD       0:00      1 (Priority)
```

`PD` is pending, with the reason in the final column, and `R` is running, in which case the node names appear there instead. A job absent from the listing has either finished or failed, and `squeue` alone does not distinguish the two; `sacct -j JOBID` does, and after the fact is the quickest way to tell a case that failed from one that ran out of time. The two call for different responses: the first wants investigation, the second merely resubmission.

A job stuck in pending is usually asking for a shape the cluster cannot conveniently supply. Reduce the task count, drop any node constraint, and check the account against the table in Section 2.

### Reading what a job wrote

The scheduler sends an electronic message when a job ends, whether it succeeded or failed, and it also leaves an output file behind. Unless the submit script directs it elsewhere, that file is named `slurm-JOBID.out` and appears in the directory the job was submitted from:

```bash
more slurm-123456.out
less runs/slurm-123456.out                   # where submit_hpc.sh puts it
less runs/array/task-123456_3.out            # task 3 of an array
```

By default Slurm combines the standard output and error streams into that one file. The scripts here separate them, sending errors to a matching `.err`, because a solver abort is far easier to find in a file of its own than buried in several thousand lines of ordinary progress output.

Three files answer different questions and it is worth knowing which to open. The `slurm-*.out` file holds what the job script itself printed, including the module and build steps. `runs/run_all.log`, or its per-task equivalent under `runs/array/`, holds the sweep runner's own record of which case started, finished or failed. The solver's own output for an individual case is in `runs/<label>/solver.out` and `solver.err`. A case reported as failed in the sweep log is diagnosed from the last of these, not the first.

### When most tasks finish in seconds

An array in which one task runs normally and the rest exit after two or three seconds has not hit a resource problem; the environment evidently works, since one task is using it. Read `runs/array/run_JOBID_N.log` first, because it distinguishes the two cases that look alike from outside. A task that never generated its deck was starved of the analysis tools; a task that generated its deck and then reported `step 0` reached the solver and the solver never started, in which case the reason is in `runs/LABEL/solver.err` and nowhere else.

The most common cause of the second is submitting the array from inside an existing allocation. The submit scripts pass `--export=ALL`, so the outer job's Slurm environment travels with the new one, and a nested `srun` that inherits another job's identifiers attempts to create its step in that job and is cancelled at once. Nothing the task script prints reveals this, since the task itself started perfectly well. The submit scripts now refuse to run when `SLURM_JOB_ID` is already set, and the task script clears the inherited step-level variables; submit from a login node.

The first has two possible collisions, both now prevented.

Two collisions are possible and both are now prevented. Building the analysis tools inside a task would have every task deleting and rewriting `bin/` and `build/` while its siblings were executing those same binaries, so `scripts/submit_hpc_array.sh` builds once before submitting and the task script only checks the result. And generated deck names encode the wavenumber and the amplitude but not the resolution or the domain size, so the three resolution cases and the outer-boundary case would all write `inputs.k0.7.eps0.02` into one directory; each task is now given a private deck directory under `runs/array/`.

The second of these is the more dangerous, because it need not crash. A task can advance a deck that another task wrote moments earlier, producing a complete and plausible record at the wrong resolution.

The symptom was obscured by the exit status. A sweep runner that logs a failed case and carries on returns success, so `sacct` reported `COMPLETED` with an exit code of `0:0` for tasks in which nothing ran. Both the runner and the task script now return a non-zero status when any case fails, so a broken sweep is recorded as failed. The internal job steps of a task, shown by `sacct` as `JOBID_N.0` and `.1`, are the individual solver launches; those being `CANCELLED` while the task itself reports `COMPLETED` is the signature of this class of fault.

Diagnosis, in order: `runs/array/run_JOBID_N.log` says which case was attempted and whether the deck was generated, `runs/array/task-JOBID_N.err` holds anything the task script itself printed, and `runs/LABEL/solver.err` holds the solver's own output for that case.

### Unable to satisfy cpu bind request

The most common way for a task to reach the solver and produce nothing is this, recorded in `runs/LABEL/solver.err`:

```
srun: error: CPU binding outside of job step allocation, allocated CPUs are: 0x140000.
srun: error: Task launch for StepId=N.0 failed on node ...: Unable to satisfy cpu bind request
srun: error: Application launch failed: Unable to satisfy cpu bind request
```

Slurm satisfies a request for eight tasks from whatever cores happen to be free, which on a busy cluster is a few cores on each of several nodes rather than eight together. The masks in the message say as much: `0x140000` is two CPUs and `0x000001` is one. The default binding cannot be mapped onto a set that ragged, so every launch is refused and the sweep reports zero steps.

The characteristic signature is that one task of an array succeeds while the rest fail, and the successful one is the one `squeue` shows on a single node. Placement, not the case, decides which.

Both runners now pass `--cpu-bind=none` when the launcher is `srun`, which leaves placement to the operating system and removes the failure. Where a run must be launched by hand inside an allocation, the same applies:

```bash
srun --cpu-bind=none -n 8 $AII/NS_fluids_lib/amr2d.gnu.MPI.ex inputs
```

Scattered placement is also slower than packed placement for a job this small, since the pressure solve communicates several times per step and those messages then cross the network. `PACK=1 scripts/submit_hpc_array.sh` asks for all ranks on one node, which runs faster and queues longer. That is a reasonable trade at eight ranks and a poor one at thirty-two.

## 12. Failures specific to this machine

| Symptom | Cause | Remedy |
|---|---|---|
| `Requested time limit is invalid (missing or exceeds some limit)` | the wall-clock time asked for is above the ceiling of the account | lower `-t`, or choose an account from Section 2 that permits it; the message concerns the time limit only |
| `Runtime environment uses unsupported PMI version PMIx` | `mpirun` belongs to a different MPI than the solver was linked against | use the launcher from that MPI, or rebuild the solver against the one on the PATH; see Section 12 of the workstation guide |
| job fails seconds after starting, with an MPI or shared-library error | the modules in the job differ from those used to build the solver | run `module list` on the build node and edit the `module load` lines in the submit scripts to match |
| `expecting num_FSI_outer_sweeps>=2` | the working tree is not at `0b8df42` | `git checkout 0b8df42` and rebuild |
| every rank reports one rank, or a message about no PMIx server | `mpirun` used in place of `srun` | the submit scripts set `LAUNCHER=srun`; do not override it on the cluster |
| job killed at the wall-clock limit with the sweep unfinished | the limit was underestimated | resubmit; completed cases are skipped and partial ones resume, provided checkpoints were written |
| a restart repeats from the beginning | no checkpoint was written within the allocation | reduce `CHECK_INT`, per Section 6 |
| a job sits pending indefinitely | a node count was specified, or too many cores were asked for on one node | request `--ntasks` only and let the scheduler place them |
| an array task fails while the others succeed | that case alone | inspect `runs/array/task-*_N.err`, then resubmit that case as a one-line case file |
| most array tasks end in seconds while one runs normally | the array was submitted from inside an allocation, or the tasks collided in the shared working directory | submit from a login node; see the end of Section 11 |
| a case logs `the solver produced no time steps` | the solver never started, as distinct from failing partway | read `runs/LABEL/solver.err`; the launcher is the usual culprit |
| `Unable to satisfy cpu bind request` | the ranks were scattered across nodes in fragmented core masks that srun's default binding cannot match | the runners now pass `--cpu-bind=none`; by hand, add it to `srun`. See Section 11 |
| far slower than expected | a debug build, or too high a resolution | check for `amr2d.gnu.DEBUG.ex`, and read the grid and step count the deck generator reports |

## 13. Afterwards

The analysis reads plain text and JSON and needs no allocation, so it belongs on a login node or a workstation rather than in a job:

```bash
cd $AII/NS_fluids_lib/LSA_Rayleigh_Fortran
scripts/merge_summaries.sh                   # after an array run
./bin/run_lsa dispersion --root runs --kr0 0.2,0.3,0.5,0.7,0.9
```

Two directories carry everything the figures need, and both are plain text. `runs/records/*.dat` holds the amplitude record of each case, which is what the growth rate is fitted from. `runs/bundles/<case>/` holds the interface profiles and field slices, written by the sweep while the plotfiles were still present, which is what the interface and field figures are drawn from. Together they are a few hundred kilobytes per case against tens of gigabytes of plotfiles.

Copying those two directories to a workstation is therefore enough to produce every figure, and the plotting stack need not be installed on the cluster. Note that the records alone are not sufficient: they carry no interface profile, so a growth curve can be drawn from them but an interface or field figure cannot. If the bundles are missing because a sweep predates them, they can be written after the fact from the plotfiles that remain, with `./bin/dump_profiles` on the machine holding those plotfiles.

```bash
# on a workstation, once both directories have been copied across
python3 plotting/plot_results.py --bundle runs/bundles/k0.7 --kr0 0.7 --outdir figures
```
