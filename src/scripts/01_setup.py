#!/usr/bin/env python3
"""Step 1 - create the directory tree and record the toolchain.

Run this first.  It fails loudly if a required executable is missing, so the
pipeline never dies halfway through for a reason as dull as a missing binary.

    python 01_setup.py
"""
from __future__ import annotations

import json
import os
import platform
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import gwaslib as G

# Tools needed for the PLINK track (step 3).  The FASTQ -> VCF track adds
# fastqc, fastp, bwa-mem2, gatk, samtools, bcftools, multiqc and mosdepth.
REQUIRED = {"plink": G.PLINK, "plink2": G.PLINK2}
# tool -> flag that prints its version.  bwa-mem2 uses `version`, not `--version`.
OPTIONAL = {"fastqc": "--version", "fastp": "--version", "samtools": "--version",
            "bcftools": "--version", "bwa-mem2": "version", "gatk": "--version",
            "multiqc": "--version", "mosdepth": "--version"}


def main() -> int:
    G.setup_logging("01_setup")
    G.LOG.info("project root: %s", G.PROJECT_ROOT)
    G.LOG.info("conda env: %s", G.CONDA_PREFIX_ENV or
               f"not set (using PATH: {os.environ.get('CONDA_PREFIX', 'no conda env active')})")

    for directory in G.ALL_DIRS:
        directory.mkdir(parents=True, exist_ok=True)
    G.LOG.info("created/verified %d directories", len(G.ALL_DIRS))

    missing = []
    versions = {"python": sys.version.split()[0], "platform": platform.platform()}

    for name, exe in REQUIRED.items():
        try:
            G.require(exe)
        except FileNotFoundError as exc:
            missing.append(str(exc))
            continue
        versions[name] = G.tool_version(exe)
        G.LOG.info("%-10s %s", name, versions[name])

    for name, flag in OPTIONAL.items():
        if G.which(name):
            versions[name] = G.tool_version(name, flag)
            G.LOG.info("%-10s %s", name, versions[name])
        else:
            G.LOG.info("%-10s not installed (only needed for the FASTQ track)", name)

    if missing:
        for message in missing:
            G.LOG.error(message)
        return 1

    record = G.RESULTS / "tool_versions.json"
    record.write_text(json.dumps(versions, indent=2))
    G.LOG.info("wrote %s", record)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
