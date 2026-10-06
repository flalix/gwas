"""Shared configuration and shell wrapper for the PennCATH GWAS pipeline.

Every external command in this project goes through `run()`, so nothing is
executed by a bare shell string and every invocation is logged with its exit
status and duration.
"""
from __future__ import annotations

import logging
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

# --------------------------------------------------------------------------
# Configuration.  Override any of these with environment variables so the same
# scripts run unchanged on a different machine.
# --------------------------------------------------------------------------
# Layout, matching the existing project:
#     <GWAS_ROOT>/data/<GWAS_DATASET>/{raw,interim,results,figures}
#     <GWAS_ROOT>/data/reference
#     <GWAS_ROOT>/logs
# The root is found by walking up from this file until a directory containing
# 'data/' appears, so the scripts can be moved to any depth (src/scripts,
# scripts, notebooks/...) without editing this line.  GWAS_ROOT overrides it.
def _find_project_root(start: Path) -> Path:
    for candidate in (start, *start.parents):
        if (candidate / "data").is_dir():
            return candidate
    raise FileNotFoundError(
        f"no directory containing 'data/' found at or above {start}. "
        f"Set GWAS_ROOT to the project root."
    )


_env_root = os.environ.get("GWAS_ROOT")
PROJECT_ROOT = (Path(_env_root) if _env_root
                else _find_project_root(Path(__file__).resolve().parent))
DATASET = os.environ.get("GWAS_DATASET", "PennCATH")

ROOT_DATA = PROJECT_ROOT / "data"
ROOT_RAW = ROOT_DATA / DATASET / "raw"
ROOT_INTERIM = ROOT_DATA / DATASET / "interim"
ROOT_RESULTS = ROOT_DATA / DATASET / "results"
ROOT_FIGURES = ROOT_DATA / DATASET / "figures"
ROOT_REFERENCE = ROOT_DATA / "reference"
ROOT_LOGS = PROJECT_ROOT / "logs"

ALL_DIRS = [ROOT_RAW, ROOT_REFERENCE, ROOT_INTERIM, ROOT_RESULTS, ROOT_FIGURES, ROOT_LOGS]

# Tool executables.  Point these at your conda env if they are not on PATH,
# e.g. PLINK=/home/flavio/miniforge3/envs/genome_env/bin/plink
PLINK = os.environ.get("PLINK", "plink")
PLINK2 = os.environ.get("PLINK2", "plink2")

# Optional: run every command inside an activated conda env.
#
#     GWAS_CONDA_PREFIX=/home/flavio/miniforge3/envs/genome_env
#
# Leave this unset if you have already run `conda activate genome_env` -- then
# PATH and JAVA_HOME are inherited and nothing extra is needed.  Set it when
# driving the pipeline from a different environment (a notebook, cron, another
# conda env), because some tools need more than PATH:  FastQC is a Perl wrapper
# that locates the JVM through JAVA_HOME, which only the env's activate.d
# scripts set.  Without it the JVM cannot find java.security and FastQC dies
# with 'java.lang.InternalError: Error loading java.security file'.
CONDA_PREFIX_ENV = os.environ.get("GWAS_CONDA_PREFIX")

# Sources the env's activate.d hooks, then execs the command.  The prefix and
# the command arrive as positional arguments, never interpolated into the
# script text.
_ACTIVATE_SHIM = (
    'prefix="$1"; shift; '
    'export CONDA_PREFIX="$prefix"; export PATH="$prefix/bin:$PATH"; '
    'for s in "$prefix"/etc/conda/activate.d/*.sh; do [ -f "$s" ] && . "$s"; done; '
    'exec "$@"'
)

# QC thresholds.  These follow the GWASTutorial / Reed et al. defaults.
QC = {
    "geno": 0.10,        # max per-variant missing rate
    "mind": 0.10,        # max per-sample missing rate
    "maf": 0.01,         # min minor allele frequency
    "hwe": 1e-10,        # HWE exact-test p-value floor (controls only)
    "het_sd": 3.0,       # heterozygosity outlier cutoff, in SD of F
    "king_cutoff": 0.0884,   # ~2nd-degree relatives
    "ld_window": 50,     # --indep-pairwise window, in variants
    "ld_step": 5,
    "ld_r2": 0.2,
    "n_pcs": 10,
}

LOG = logging.getLogger("gwas")


# --------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------
def setup_logging(name: str) -> Path:
    """Log to both the console and logs/<name>.log.  Returns the log path."""
    ROOT_LOGS.mkdir(parents=True, exist_ok=True)
    logfile = ROOT_LOGS / f"{name}.log"
    LOG.handlers.clear()
    LOG.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%H:%M:%S")
    for handler in (logging.StreamHandler(sys.stdout), logging.FileHandler(logfile)):
        handler.setFormatter(fmt)
        LOG.addHandler(handler)
    LOG.info("logging to %s", logfile)
    return logfile


# --------------------------------------------------------------------------
# Shell wrapper
# --------------------------------------------------------------------------
class CommandFailed(RuntimeError):
    pass


def run(cmd: list[str], *, cwd: Path | None = None, check: bool = True,
        quiet: bool = True) -> subprocess.CompletedProcess:
    """Run an external command given as a list of arguments.

    No shell is involved, so filenames containing spaces are safe and nothing
    is interpolated into a command line.
    """
    cmd = [str(c) for c in cmd]
    LOG.info("$ %s", " ".join(cmd))
    if CONDA_PREFIX_ENV:
        cmd = ["bash", "-c", _ACTIVATE_SHIM, "gwas-run", CONDA_PREFIX_ENV, *cmd]
    start = time.time()
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    elapsed = time.time() - start

    if proc.returncode != 0:
        LOG.error("exit %d after %.1fs", proc.returncode, elapsed)
        LOG.error("stdout tail:\n%s", "\n".join(proc.stdout.splitlines()[-25:]))
        LOG.error("stderr tail:\n%s", "\n".join(proc.stderr.splitlines()[-25:]))
        if check:
            raise CommandFailed(f"{cmd[0]} exited {proc.returncode}")
    else:
        LOG.info("ok (%.1fs)", elapsed)
        if not quiet:
            LOG.info(proc.stdout)
    return proc


def pipe(stages: list[list[str]], *, stdout: Path | None = None,
         check: bool = True) -> int:
    """Run `a | b | c`, returning the exit status of the LAST stage.

    Used where writing the intermediate to disk would be wasteful -- aligning
    (bwa mem | samtools sort) and CRAM-to-FASTQ (samtools collate | samtools
    fastq) both move tens of GB that never need to land.

    Every stage's status is checked, not just the last, so a failure in an
    upstream stage cannot be masked by a downstream one exiting 0.
    """
    rendered = " | ".join(" ".join(str(a) for a in s) for s in stages)
    LOG.info("$ %s%s", rendered, f" > {stdout}" if stdout else "")
    start = time.time()

    if CONDA_PREFIX_ENV:
        stages = [["bash", "-c", _ACTIVATE_SHIM, "gwas-run", CONDA_PREFIX_ENV,
                   *[str(a) for a in s]] for s in stages]

    procs: list[subprocess.Popen] = []
    try:
        previous = None
        for index, stage in enumerate(stages):
            last = index == len(stages) - 1
            out = (open(stdout, "wb") if (last and stdout)
                   else (None if last else subprocess.PIPE))
            proc = subprocess.Popen([str(a) for a in stage],
                                    stdin=previous.stdout if previous else None,
                                    stdout=out, stderr=subprocess.PIPE)
            if previous:
                # let the upstream stage see SIGPIPE if this one exits early
                previous.stdout.close()
            procs.append(proc)
            previous = proc
        statuses = [(p.wait(), p) for p in procs]
    finally:
        for proc in procs:
            if proc.poll() is None:
                proc.kill()

    elapsed = time.time() - start
    failed = [(code, p) for code, p in statuses if code != 0]
    if failed:
        for code, proc in failed:
            tail = (proc.stderr.read().decode(errors="replace").splitlines()[-15:]
                    if proc.stderr else [])
            LOG.error("stage %s exited %d: %s", proc.args[0], code, "\n".join(tail))
        if check:
            raise CommandFailed(f"pipeline failed after {elapsed:.1f}s: {rendered}")
    else:
        LOG.info("ok (%.1fs)", elapsed)
    return statuses[-1][0]


def which(tool: str) -> str | None:
    """Locate an executable, looking inside GWAS_CONDA_PREFIX/bin first."""
    if CONDA_PREFIX_ENV:
        candidate = Path(CONDA_PREFIX_ENV) / "bin" / tool
        if candidate.exists():
            return str(candidate)
    return shutil.which(tool)


def require(tool: str) -> str:
    """Abort early with a useful message if an executable is missing."""
    path = which(tool)
    if path is None:
        raise FileNotFoundError(
            f"'{tool}' not found on PATH. Activate the conda env that provides it, "
            f"set GWAS_CONDA_PREFIX to that env, or set the matching variable "
            f"(e.g. PLINK=/path/to/plink)."
        )
    return path


def tool_version(tool: str, version_flag: str = "--version") -> str:
    """First line of `<tool> --version`, for the provenance record.

    Routed through run() so the conda env is activated if configured -- FastQC
    in particular reports a JVM error otherwise.
    """
    try:
        proc = run([tool, version_flag], check=False, quiet=True)
        out = (proc.stdout or proc.stderr).strip().splitlines()
        return out[0] if out else "unknown"
    except Exception as exc:                                   # noqa: BLE001
        return f"unavailable ({exc})"


# --------------------------------------------------------------------------
# PLINK convenience
# --------------------------------------------------------------------------
def plink(args: list[str], *, out: Path, v2: bool = False) -> Path:
    """Run PLINK with --out <out> and return the output stem.

    `args` must NOT contain --out; it is appended here so every call writes a
    .log next to its outputs.
    """
    exe = PLINK2 if v2 else PLINK
    out.parent.mkdir(parents=True, exist_ok=True)
    run([exe, *args, "--out", out])
    return out


def count_bfile(stem: Path) -> tuple[int, int]:
    """(variants, samples) in a PLINK binary fileset, read from .bim/.fam."""
    n_var = sum(1 for _ in open(f"{stem}.bim"))
    n_sam = sum(1 for _ in open(f"{stem}.fam"))
    return n_var, n_sam
