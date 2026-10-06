#!/usr/bin/env python3
"""Step 5 - fetch everything needed for a whole-genome FASTQ -> VCF run.

Sample is NA12878 (HG001), deliberately: it is the Genome in a Bottle reference
sample, so the final VCF can be scored against a published truth set instead of
merely inspected.

    python 05_fetch_wgs_data.py                 # print the budget, fetch nothing
    python 05_fetch_wgs_data.py --run           # fetch everything (hours)
    python 05_fetch_wgs_data.py --run --only reference,knownsites
    python 05_fetch_wgs_data.py --run --threads 16

Every stage is skipped if its output already exists and is non-empty, so an
interrupted run is resumed by re-issuing the same command.  Downloads use
`wget -c`, which resumes mid-file.

Stages
    reference   GRCh38 analysis-set FASTA (3.26 GB) + .fai + .dict + bwa index.
                bwa, not bwa-mem2: measured on chr1, bwa-mem2's index needs
                ~24 bytes/base => ~75 GB RAM for the genome, over this machine's
                63 GB.  bwa needs 1.5 bytes/base => ~4.7 GB, and ~31 min.
    cram        NA12878 30x CRAM (15.8 GB), then name-collated and converted to
                paired FASTQ (~45 GB gz).  Collation needs ~30 GB of scratch in
                $TMPDIR: samtools fastq cannot pair reads from a
                coordinate-sorted file.
    knownsites  Ensembl per-chromosome variation VCFs for chr1-22,X,Y (12.0 GB).
                Ensembl names contigs '1'..'22','X','Y'; the reference uses
                'chr1'..'chrY', so each file is renamed before being concatenated
                into one BQSR known-sites VCF.  Its header also lacks ##contig
                lines, so the published .csi must be fetched alongside or
                bcftools rejects every record.
    giab        HG001 truth VCF + high-confidence BED, for scoring (0.14 GB).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import gwaslib

SAMPLE = "NA12878"
CRAM_URL = "https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram"
CRAM_BYTES = 15_797_182_294

REF_URL = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/"
           "GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa")
REF_BYTES = 3_263_683_042

ENSEMBL = "https://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens"
CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]

GIAB = ("https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/"
        "NISTv4.2.1/GRCh38")
GIAB_FILES = {
    "HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz": 125_932_193,
    "HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi": 1_589_208,
    "HG001_GRCh38_1_22_v4.2.1_benchmark.bed": 15_479_939,
}

BUDGET = """
Downloads
  GRCh38 analysis-set FASTA                3.26 GB
  NA12878 30x CRAM                        15.80 GB
  Ensembl known sites, chr1-22,X,Y        11.99 GB   (24 files + .csi each)
  GIAB HG001 truth set                     0.14 GB
  ----------------------------------------------
  total                                   31.20 GB

Disk after this script finishes
  FASTQ R1+R2 gz                          45.0 GB
  reference + .fai + .dict                 3.3 GB
  bwa index                                5.4 GB
  known sites, concatenated               12.1 GB
  CRAM (kept; FASTQ is re-derivable)      15.8 GB
  ----------------------------------------------
  ~82 GB, plus ~30 GB transient scratch in $TMPDIR during collation

Time, dominated by transfer rate
  downloads        depends on your link; 31 GB
  bwa index        ~31 min   (single-threaded, ~4.7 GB RAM)
  CRAM -> FASTQ    ~1-2 h    (collate + convert)
  knownsites       ~30-60 min (rename + concatenate 12 GB)
"""


def usable(path: Path) -> bool:
    if path.exists() and path.stat().st_size == 0:
        gwaslib.LOG.info("discarding empty %s", path.name)
        path.unlink()
    return path.exists()


def download(url: str, dest: Path, expect: int | None = None) -> Path:
    if usable(dest) and (expect is None or dest.stat().st_size == expect):
        gwaslib.LOG.info("%s already complete", dest.name)
        return dest
    gwaslib.run(["wget", "-c", "-O", dest, url])
    if expect is not None and dest.stat().st_size != expect:
        raise ValueError(f"{dest.name}: {dest.stat().st_size} bytes, expected {expect}")
    return dest


# --------------------------------------------------------------------------
def stage_reference(threads: int) -> Path:
    fasta = gwaslib.ROOT_REFERENCE / "GRCh38_full_analysis_set_plus_decoy_hla.fa"
    download(REF_URL, fasta, REF_BYTES)
    if not usable(Path(f"{fasta}.fai")):
        gwaslib.run(["samtools", "faidx", fasta])
    dict_path = fasta.with_suffix(".dict")
    if not usable(dict_path):
        gwaslib.run(["gatk", "CreateSequenceDictionary", "-R", fasta, "-O", dict_path])
    if not usable(Path(f"{fasta}.bwt")):
        gwaslib.LOG.info("building bwa index: ~31 min, ~4.7 GB RAM, 5.4 GB output")
        gwaslib.run(["bwa", "index", fasta])
    return fasta


def stage_cram(fasta: Path, threads: int) -> tuple[Path, Path]:
    cram = gwaslib.ROOT_RAW / f"{SAMPLE}.final.cram"
    download(CRAM_URL, cram, CRAM_BYTES)

    r1 = gwaslib.ROOT_RAW / f"{SAMPLE}_R1.fastq.gz"
    r2 = gwaslib.ROOT_RAW / f"{SAMPLE}_R2.fastq.gz"
    if usable(r1) and usable(r2):
        gwaslib.LOG.info("FASTQ already present")
        return r1, r2

    # samtools fastq needs mates adjacent; the CRAM is coordinate-sorted, so
    # collate first.  -u keeps the intermediate uncompressed for speed; the
    # temp prefix must sit somewhere with ~30 GB free.
    scratch = Path(os.environ.get("TMPDIR", "/tmp")) / f"collate_{SAMPLE}"
    gwaslib.LOG.info("collating + converting to FASTQ (scratch prefix %s)", scratch)
    collate = ["samtools", "collate", "-@", str(threads), "-u", "-O",
               "--reference", fasta, cram, str(scratch)]
    fastq = ["samtools", "fastq", "-@", str(threads), "-1", str(r1), "-2", str(r2),
             "-0", "/dev/null", "-s", "/dev/null", "-n", "-"]
    gwaslib.pipe([collate, fastq])
    return r1, r2


def stage_knownsites(fasta: Path, threads: int) -> Path:
    merged = gwaslib.ROOT_REFERENCE / "ensembl_known_sites_GRCh38.vcf.gz"
    if usable(merged):
        gwaslib.LOG.info("%s already present", merged.name)
        return merged

    rename = gwaslib.ROOT_REFERENCE / "rename_chrs_all.txt"
    rename.write_text("".join(f"{c}\tchr{c}\n" for c in CHROMS))

    parts = []
    for chrom in CHROMS:
        raw = gwaslib.ROOT_REFERENCE / f"homo_sapiens-chr{chrom}.vcf.gz"
        out = gwaslib.ROOT_REFERENCE / f"known_chr{chrom}.vcf.gz"
        if not usable(out):
            download(f"{ENSEMBL}/homo_sapiens-chr{chrom}.vcf.gz", raw)
            download(f"{ENSEMBL}/homo_sapiens-chr{chrom}.vcf.gz.csi", Path(f"{raw}.csi"))
            gwaslib.run(["bcftools", "annotate", "--rename-chrs", rename, "--threads", str(threads), "-Oz", "-o", out, raw])
            gwaslib.run(["bcftools", "index", "-t", "--threads", str(threads), out])
            raw.unlink(missing_ok=True)
            Path(f"{raw}.csi").unlink(missing_ok=True)
        parts.append(out)

    gwaslib.LOG.info("concatenating %d per-chromosome files", len(parts))
    gwaslib.run(["bcftools", "concat", "--threads", str(threads), "-Oz", "-o", merged, *parts])
    gwaslib.run(["bcftools", "index", "-t", "--threads", str(threads), merged])
    for part in parts:
        part.unlink(missing_ok=True)
        Path(f"{part}.tbi").unlink(missing_ok=True)
    return merged


def stage_giab() -> dict[str, str]:
    out = {}
    for name, expect in GIAB_FILES.items():
        out[name] = str(download(f"{GIAB}/{name}", gwaslib.ROOT_REFERENCE / name, expect))
    return out


# --------------------------------------------------------------------------
STAGES = ["reference", "cram", "knownsites", "giab"]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--threads", type=int, default=min(16, os.cpu_count() or 4))
    parser.add_argument("--only", default=",".join(STAGES),
                        help=f"comma-separated subset of {STAGES}")
    args = parser.parse_args()

    gwaslib.setup_logging("05_fetch_wgs_data")
    for directory in gwaslib.ALL_DIRS:
        directory.mkdir(parents=True, exist_ok=True)

    if not args.run:
        print(BUDGET)
        print("Re-run with --run.  Everything is resumable; re-issue the same command\n"
              "after an interruption and completed stages are skipped.\n")
        return 0

    wanted = [s.strip() for s in args.only.split(",")]
    unknown = set(wanted) - set(STAGES)
    if unknown:
        raise SystemExit(f"unknown stage(s): {sorted(unknown)}; choose from {STAGES}")

    for tool in ("wget", "samtools", "bcftools", "bwa", "gatk"):
        gwaslib.require(tool)

    manifest: dict[str, object] = {"sample": SAMPLE, "build": "GRCh38",
                                   "threads": args.threads}
    fasta = gwaslib.ROOT_REFERENCE / "GRCh38_full_analysis_set_plus_decoy_hla.fa"

    if "reference" in wanted:
        fasta = stage_reference(args.threads)
        manifest["reference"] = str(fasta)
    if "cram" in wanted:
        r1, r2 = stage_cram(fasta, args.threads)
        manifest["fastq"] = [str(r1), str(r2)]
        manifest["fastq_bytes"] = [r1.stat().st_size, r2.stat().st_size]
    if "knownsites" in wanted:
        manifest["known_sites"] = str(stage_knownsites(fasta, args.threads))
    if "giab" in wanted:
        manifest["giab"] = stage_giab()

    path = gwaslib.ROOT_RESULTS / "wgs_data_manifest.json"
    existing = json.loads(path.read_text()) if path.exists() else {}
    existing.update(manifest)
    path.write_text(json.dumps(existing, indent=2))
    gwaslib.LOG.info("manifest: %s", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
