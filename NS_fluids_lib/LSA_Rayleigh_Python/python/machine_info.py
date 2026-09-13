"""Machine and run provenance, recorded alongside every timing measurement.

A cost figure without the machine it was measured on is not reproducible and not comparable. A reader cannot tell whether a Jacobian-vector product taking 40 seconds indicates an efficient method or a fast processor, and a second measurement taken a year later on different hardware cannot be placed beside the first. This module collects what is needed to make a timing meaningful and writes it into the record the timing script saves, so that provenance travels with the number rather than being reconstructed from memory.

Everything here is read from the operating system and degrades gracefully: a field that cannot be determined is reported as unknown rather than raising, since a missing processor model is not a reason for a benchmark of several hours to fail at the point of writing its answer. Nothing outside the standard library is required, so the module runs on a login node with no scientific stack installed.

The distinction between processors present and processors used matters more than it appears. A cluster node commonly reports many cores while a Slurm allocation grants a few, and timing that assumes the former overstates the parallel efficiency of everything measured on it. Both are recorded, together with the CPU affinity mask actually in force, which is what the process can use rather than what the scheduler was asked for.
"""

from __future__ import annotations

import json
import os
import platform
import re
import subprocess
import sys
import time

UNKNOWN = "unknown"


def _read(path):
    """Contents of a file, or None if it cannot be read."""
    try:
        with open(path) as fh:
            return fh.read()
    except OSError:
        return None


def _run(command):
    """Standard output of a command, or None if it cannot be run."""
    try:
        out = subprocess.run(command, capture_output=True, text=True,
                             timeout=10, check=False)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout if out.returncode == 0 else None


def processor_model():
    """Human-readable processor model."""
    info = _read("/proc/cpuinfo")
    if info:
        match = re.search(r"^model name\s*:\s*(.+)$", info, re.MULTILINE)
        if match:
            return match.group(1).strip()
        # Some architectures, notably POWER and certain ARM parts, use other keys.
        for key in ("cpu", "Processor", "Hardware"):
            match = re.search(rf"^{key}\s*:\s*(.+)$", info, re.MULTILINE)
            if match:
                return match.group(1).strip()
    out = _run(["sysctl", "-n", "machdep.cpu.brand_string"])
    if out:
        return out.strip()
    return platform.processor() or UNKNOWN


def physical_cores():
    """Physical cores present, distinguished from hardware threads.

    Hyper-threading inflates the logical count by a factor of two on most Intel parts, and a scaling study reported against logical processors shows an apparent collapse in efficiency past the physical core count that has nothing to do with the algorithm. The pairing of physical identifier and core identifier counts each physical core once.
    """
    info = _read("/proc/cpuinfo")
    if info:
        pairs = set()
        physical = core = None
        for line in info.splitlines():
            if line.startswith("physical id"):
                physical = line.split(":")[1].strip()
            elif line.startswith("core id"):
                core = line.split(":")[1].strip()
                if physical is not None:
                    pairs.add((physical, core))
        if pairs:
            return len(pairs)
    out = _run(["sysctl", "-n", "hw.physicalcpu"])
    if out and out.strip().isdigit():
        return int(out.strip())
    return None


def available_cores():
    """Processors this process may actually use.

    `sched_getaffinity` reports the mask in force, which under a scheduler is the allocation rather than the machine. `os.cpu_count` reports the machine. Where the two disagree, the affinity mask is the honest denominator for a parallel efficiency.
    """
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count()


def memory_gb():
    """Total physical memory in gibibytes."""
    info = _read("/proc/meminfo")
    if info:
        match = re.search(r"^MemTotal:\s*(\d+)\s*kB$", info, re.MULTILINE)
        if match:
            return round(int(match.group(1)) / (1024.0 * 1024.0), 1)
    out = _run(["sysctl", "-n", "hw.memsize"])
    if out and out.strip().isdigit():
        return round(int(out.strip()) / (1024.0 ** 3), 1)
    return None


def scheduler_allocation():
    """Slurm allocation details, when running inside a job.

    Recorded because a timing taken inside an allocation is bounded by what the allocation granted, and because the job identifier is what ties a figure back to a specific run in the cluster accounting records.
    """
    keys = {
        "job_id": "SLURM_JOB_ID",
        "job_name": "SLURM_JOB_NAME",
        "partition": "SLURM_JOB_PARTITION",
        "nodes": "SLURM_JOB_NUM_NODES",
        "ntasks": "SLURM_NTASKS",
        "cpus_per_task": "SLURM_CPUS_PER_TASK",
        "nodelist": "SLURM_JOB_NODELIST",
    }
    found = {name: os.environ[var] for name, var in keys.items()
             if os.environ.get(var)}
    return found or None


def mpi_ranks():
    """MPI world size, taken from whichever launcher variable is set.

    Open MPI, MPICH and Slurm each advertise the rank count under a different name, and a benchmark launched under one must not silently report the default of one because it looked only for another.
    """
    for var in ("OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "MPI_LOCALNRANKS",
                "SLURM_NTASKS"):
        value = os.environ.get(var)
        if value and value.isdigit():
            return int(value)
    return 1


def compiler_version(command="gfortran"):
    """First line of the compiler's version banner."""
    out = _run([command, "--version"])
    return out.splitlines()[0].strip() if out else None


def solver_commit(path):
    """Commit of the solver working tree, which pins the numerics.

    A timing is only comparable against another taken from the same source. The commit is recorded rather than described, since a description of a local edit cannot be checked.
    """
    if not path or not os.path.isdir(path):
        return None
    out = _run(["git", "-C", path, "rev-parse", "--short", "HEAD"])
    if not out:
        return None
    commit = out.strip()
    dirty = _run(["git", "-C", path, "status", "--porcelain"])
    if dirty and dirty.strip():
        # A modified tree is not the commit it claims to be, and a timing taken from one must say so.
        commit += " (modified)"
    return commit


def collect(solver_path=None, extra=None):
    """Everything above, as one dictionary ready to be stored."""
    physical = physical_cores()
    logical = os.cpu_count()
    usable = available_cores()

    record = {
        "recorded": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
        "hostname": platform.node() or UNKNOWN,
        "system": f"{platform.system()} {platform.release()}",
        "machine": platform.machine(),
        "processor": processor_model(),
        "cores_physical": physical,
        "cores_logical": logical,
        "cores_available": usable,
        "memory_gb": memory_gb(),
        "mpi_ranks": mpi_ranks(),
        "python": sys.version.split()[0],
        "compiler": compiler_version(),
        "scheduler": scheduler_allocation(),
        "solver_commit": solver_commit(solver_path),
    }
    if extra:
        record.update(extra)
    return record


def format_table(record, width=22):
    """The record as aligned plain text, in the house report style."""
    lines = ["-" * 66,
             f"{'machine and run provenance':<{width}}",
             "-" * 66]
    for key, value in record.items():
        if value is None:
            value = UNKNOWN
        if isinstance(value, dict):
            value = ", ".join(f"{k}={v}" for k, v in value.items())
        lines.append(f"{key:<{width}} {value}")
    lines.append("-" * 66)
    return "\n".join(lines)


def caption(record):
    """A one-sentence provenance note for a figure caption.

    A cost figure carries this so that the machine is stated where the number is read, rather than in a file the reader does not have.
    """
    cores = record.get("cores_available") or UNKNOWN
    ranks = record.get("mpi_ranks") or 1
    return (f"Measured on {record.get('hostname', UNKNOWN)}, "
            f"{record.get('processor', UNKNOWN)}, "
            f"{cores} processors available, {ranks} MPI rank(s).")


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    solver = argv[0] if argv else None
    record = collect(solver_path=solver)
    print(format_table(record))
    print()
    print(json.dumps(record, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
