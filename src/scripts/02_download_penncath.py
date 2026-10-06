#!/usr/bin/env python3
"""Step 2 - fetch and verify the PennCATH dataset.

The archive is 143 MB, so by default this script only PRINTS the command and
exits.  Add --run to let it download, or run the printed wget yourself.

    python 02_download_penncath.py            # show the command, download nothing
    python 02_download_penncath.py --run      # download, checksum and extract

Expected contents (verified against the live archive on 2026-10-02):
    penncath.bed  302,377,026 B   penncath.bim  24,053,xxx B
    penncath.fam       29,xxx B   penncath.csv       38,xxx B
    1,401 samples x 861,473 autosomal variants, build GRCh37/hg19
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import tarfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import gwaslib as G

URL = "https://d1ypx1ckp5bo16.cloudfront.net/penncath/penncath.tar.gz"
ARCHIVE_BYTES = 142_819_504
EXPECTED = {"penncath.bed": 302_377_026, "penncath.bim": None,
            "penncath.fam": None, "penncath.csv": None}
EXPECTED_SAMPLES = 1401
EXPECTED_VARIANTS = 861_473


def sha256(path, chunk=1 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true",
                        help="actually download (143 MB) instead of printing the command")
    args = parser.parse_args()

    G.setup_logging("02_download")
    G.RAW.mkdir(parents=True, exist_ok=True)
    archive = G.RAW / "penncath.tar.gz"

    if not args.run and not archive.exists():
        print("\nRun this yourself, then re-run with --run to verify and extract:\n")
        print(f"  cd {G.RAW}")
        print(f"  wget -c {URL}\n")
        print(f"  # expected size: {ARCHIVE_BYTES:,} bytes\n")
        return 0

    if not archive.exists():
        # -c resumes a partial download; the server supports range requests.
        G.run(["wget", "-c", "-O", archive, URL])

    size = archive.stat().st_size
    if size != ARCHIVE_BYTES:
        G.LOG.error("size mismatch: got %d bytes, expected %d. Delete the file and retry.",
                    size, ARCHIVE_BYTES)
        return 1
    G.LOG.info("archive size ok (%d bytes)", size)

    checksum = sha256(archive)
    G.LOG.info("sha256 %s", checksum)

    # filter="data" refuses absolute paths and symlinks escaping the target dir.
    with tarfile.open(archive) as tar:
        tar.extractall(G.RAW, filter="data")
    G.LOG.info("extracted into %s", G.RAW)

    # The tarball stores files under data/, so flatten into data/raw/.
    nested = G.RAW / "data"
    if nested.is_dir():
        for item in nested.iterdir():
            item.rename(G.RAW / item.name)
        nested.rmdir()

    for name, expected_size in EXPECTED.items():
        path = G.RAW / name
        if not path.exists():
            G.LOG.error("missing after extraction: %s", name)
            return 1
        actual = path.stat().st_size
        if expected_size is not None and actual != expected_size:
            G.LOG.error("%s is %d bytes, expected %d", name, actual, expected_size)
            return 1
        G.LOG.info("%-13s %12d bytes", name, actual)

    n_var, n_sam = G.count_bfile(G.RAW / "penncath")
    G.LOG.info("fileset: %d variants x %d samples", n_var, n_sam)
    if (n_var, n_sam) != (EXPECTED_VARIANTS, EXPECTED_SAMPLES):
        G.LOG.error("dimension mismatch; expected %d x %d",
                    EXPECTED_VARIANTS, EXPECTED_SAMPLES)
        return 1

    (G.RESULTS / "penncath_download.json").write_text(json.dumps(
        {"url": URL, "bytes": size, "sha256": checksum,
         "variants": n_var, "samples": n_sam, "build": "GRCh37/hg19"}, indent=2))
    G.LOG.info("verified. provenance written to results/penncath_download.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
