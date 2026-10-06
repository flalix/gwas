#!/usr/bin/env python3
"""Step 4 - assemble a small, real human FASTQ test set on GRCh38 chr20.

Nothing here downloads a whole-genome file.  The strategy is to take byte
ranges and region slices out of public files instead:

  reference   chr20 is cut out of the 1000 Genomes GRCh38 analysis-set FASTA
              with an HTTP range request (65 MB instead of 3.26 GB), using the
              byte offset in the published .fai.  The result is verified
              against the M5 checksum in the reference's own chr20 header.
  reads       a chr20 region is sliced straight out of a 30x 1000 Genomes CRAM
              over https (samtools fetches only the blocks it needs via the
              .crai) and converted back to paired FASTQ.  A 15.4 GB CRAM yields
              ~320 MB of transfer for all of chr20.
  known sites Ensembl ships per-chromosome variation VCFs, so BQSR needs the
              chr20 file (264 MB) rather than the 10 GB genome-wide dbSNP.

    python 04_fetch_fastq_testdata.py                  # show the plan, fetch nothing
    python 04_fetch_fastq_testdata.py --run            # all of chr20
    python 04_fetch_fastq_testdata.py --run --region chr20:1-10000000   # quick pass

Why chr20 and not the whole genome: a `bwa-mem2 index` of chr20 peaks at
1.56 GB of RAM and takes 13 s.  That scales roughly linearly with reference
length, so the full 3.1 Gb genome would need on the order of 75 GB (RAM) -- 
more than this machine's RAM 63 GB.  
Restricting to one chromosome keeps the index buildable
locally and the whole pipeline runnable in minutes.

Caveat worth stating in any write-up: the reads are extracted from an existing
alignment, so they are the reads that the 1000 Genomes pipeline already mapped
to chr20.  They are real sequencer reads with real quality strings, which is
what the FASTQ QC step needs, but they are not a raw unfiltered run, and
aligning them against a chr20-only reference will mismap the fraction whose
true origin is a homologous region elsewhere in the genome.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import gwaslib

REF_BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome"
REF_FASTA = f"{REF_BASE}/GRCh38_full_analysis_set_plus_decoy_hla.fa"
REF_FAI = f"{REF_FASTA}.fai"

SEQ_INDEX = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index"

ENSEMBL_VCF = "https://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/homo_sapiens-chr20.vcf.gz"

CONTIG = "chr20"
CONTIG_M5 = "b18e6c531b0bd70e949a7fc20859cb01"   # from the reference's chr20 header
DEFAULT_SAMPLE = "NA12718"


def usable(path: Path) -> bool:
    """A file left behind by a failed step is present but empty; treat as absent."""
    if path.exists() and path.stat().st_size == 0:
        gwaslib.LOG.info("discarding empty %s from an earlier failed run", path.name)
        path.unlink()
    return path.exists()


def http_get(url: str, byte_range: tuple[int, int] | None = None) -> bytes:
    request = urllib.request.Request(url)
    if byte_range:
        request.add_header("Range", f"bytes={byte_range[0]}-{byte_range[1]}")
    with urllib.request.urlopen(request, timeout=300) as response:
        return response.read()


def content_length(url: str) -> int | None:
    request = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            value = response.headers.get("content-length")
            return int(value) if value else None
    except Exception:                                          # noqa: BLE001
        return None


# --------------------------------------------------------------------------
def fetch_contig_fasta(out_fasta: Path) -> Path:
    """Cut one contig out of the remote FASTA using the published .fai offsets."""
    if usable(out_fasta):
        gwaslib.LOG.info("%s already present", out_fasta.name)
        return out_fasta

    fai_text = http_get(REF_FAI).decode()
    row = next((line.split("\t") for line in fai_text.splitlines()
                if line.split("\t")[0] == CONTIG), None)
    if row is None:
        raise ValueError(f"{CONTIG} not found in {REF_FAI}")

    # .fai columns: name, length, byte offset of first base, bases/line, bytes/line
    length, offset, bases_per_line = int(row[1]), int(row[2]), int(row[3])
    newlines = -(-length // bases_per_line)          # ceil division
    n_bytes = length + newlines
    gwaslib.LOG.info("%s: %d bases, fetching %d bytes from offset %d",
               CONTIG, length, n_bytes, offset)

    sequence = http_get(REF_FASTA, (offset, offset + n_bytes - 1))
    out_fasta.write_bytes(f">{CONTIG}\n".encode() + sequence)

    digest = hashlib.md5(sequence.replace(b"\n", b"")).hexdigest()
    if digest != CONTIG_M5:
        out_fasta.unlink()
        raise ValueError(f"checksum mismatch: got {digest}, expected {CONTIG_M5}")
    gwaslib.LOG.info("checksum ok (M5 %s)", digest)
    return out_fasta


def resolve_cram(sample: str) -> str:
    """Look the sample's CRAM up in the 1000 Genomes sequence index."""
    gwaslib.LOG.info("resolving CRAM for %s", sample)
    index_text = http_get(SEQ_INDEX).decode()
    for line in index_text.splitlines():
        if line.startswith("#"):
            continue
        path = line.split("\t")[0]
        if path.endswith(".cram") and f"/{sample}." in path:
            # the index lists ftp:// URLs; the same paths serve over https,
            # which supports the range requests samtools needs
            return path.replace("ftp://ftp.sra.ebi.ac.uk", "https://ftp.sra.ebi.ac.uk")
    raise ValueError(f"no CRAM for sample {sample} in the 1000G index")


def slice_to_fastq(cram_url: str, fasta: Path, region: str, sample: str) -> tuple[Path, Path]:
    """Region-slice a remote CRAM and write collated paired FASTQ."""
    r1 = gwaslib.ROOT_RAW / f"{sample}_{region.replace(':', '_')}_R1.fastq.gz"
    r2 = gwaslib.ROOT_RAW / f"{sample}_{region.replace(':', '_')}_R2.fastq.gz"
    if usable(r1) and usable(r2):
        gwaslib.LOG.info("%s / %s already present", r1.name, r2.name)
        return r1, r2

    bam = gwaslib.ROOT_INTERIM / f"{sample}_{region.replace(':', '_')}.bam"
    # -M keeps only reads overlapping the region; collate groups mates together,
    # which samtools fastq needs in order to split R1/R2 correctly.
    gwaslib.run(["samtools", "view", "-T", fasta, "-b", "-o", bam, cram_url, region])
    collated = gwaslib.ROOT_INTERIM / f"{bam.stem}.collated.bam"
    gwaslib.run(["samtools", "collate", "-@", "4", "-o", collated, bam])
    gwaslib.run(["samtools", "fastq", "-@", "4", "-1", r1, "-2", r2,
           "-0", "/dev/null", "-s", "/dev/null", "-n", collated])
    for path in (bam, collated):
        path.unlink(missing_ok=True)
    return r1, r2


def fetch_known_sites(fasta: Path) -> Path:
    """Ensembl chr20 variation VCF, with contigs renamed to match the reference."""
    out = gwaslib.ROOT_REFERENCE / f"ensembl_{CONTIG}_known_sites.vcf.gz"
    if usable(out):
        gwaslib.LOG.info("%s already present", out.name)
        return out

    raw = gwaslib.ROOT_REFERENCE / "homo_sapiens-chr20.vcf.gz"
    if not usable(raw):
        gwaslib.run(["wget", "-c", "-O", raw, ENSEMBL_VCF])

    # The Ensembl VCF header declares no ##contig lines, so bcftools rejects
    # every record as "Contig '20' is not defined in the header".  Its index
    # carries the contig list, which is enough to resolve that -- so fetch the
    # published .csi rather than spending minutes rebuilding one.
    csi = Path(f"{raw}.csi")
    if not usable(csi):
        gwaslib.run(["wget", "-c", "-O", csi, f"{ENSEMBL_VCF}.csi"])

    # Ensembl names the contig '20'; the reference and the CRAM use 'chr20'.
    rename = gwaslib.ROOT_REFERENCE / "rename_chrs.txt"
    rename.write_text(f"20\t{CONTIG}\n")
    gwaslib.run(["bcftools", "annotate", "--rename-chrs", rename, "-Oz", "-o", out, raw])
    gwaslib.run(["bcftools", "index", "-t", out])
    return out


def build_indexes(fasta: Path) -> None:
    if not Path(f"{fasta}.fai").exists():
        gwaslib.run(["samtools", "faidx", fasta])
    dict_path = fasta.with_suffix(".dict")
    if not dict_path.exists():
        gwaslib.run(["gatk", "CreateSequenceDictionary", "-R", fasta, "-O", dict_path])
    if not Path(f"{fasta}.bwt.2bit.64").exists():
        gwaslib.LOG.info("building bwa-mem2 index (~1.6 GB RAM, ~15 s for chr20)")
        gwaslib.run(["bwa-mem2", "index", fasta])


# --------------------------------------------------------------------------
def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true", help="fetch, instead of printing the plan")
    parser.add_argument("--region", default=CONTIG, help="e.g. chr20 or chr20:1-10000000")
    parser.add_argument("--sample", default=DEFAULT_SAMPLE, help="1000 Genomes sample ID")
    args = parser.parse_args()

    gwaslib.setup_logging("04_fetch_fastq_testdata")
    for directory in gwaslib.ALL_DIRS:
        directory.mkdir(parents=True, exist_ok=True)

    if not args.run:
        vcf_bytes = content_length(ENSEMBL_VCF)
        print(f"""
Planned transfers for region {args.region}, sample {args.sample}:

  {CONTIG} reference slice     ~65 MB   (byte range out of a 3.26 GB FASTA)
  CRAM region slice           ~320 MB   for all of chr20; ~50 MB for 10 Mb
  Ensembl chr20 known sites   {(vcf_bytes or 0) / 1e6:>6.0f} MB   (+ ~47 KB index)
  ------------------------------------------------------------------
  total download              ~650 MB
  generated locally           ~1.3 GB paired FASTQ + ~354 MB bwa-mem2 index

Re-run with --run to fetch.  Add --region chr20:1-10000000 for a faster first
pass (~115 MB download, ~200 MB FASTQ).
""")
        return 0

    gwaslib.require("samtools")
    gwaslib.require("bcftools")
    gwaslib.require("bwa-mem2")

    fasta = fetch_contig_fasta(gwaslib.ROOT_REFERENCE / f"GRCh38_{CONTIG}.fa")
    build_indexes(fasta)
    cram_url = resolve_cram(args.sample)
    gwaslib.LOG.info("CRAM: %s", cram_url)
    r1, r2 = slice_to_fastq(cram_url, fasta, args.region, args.sample)
    known = fetch_known_sites(fasta)

    manifest = {
        "region": args.region, "sample": args.sample, "build": "GRCh38",
        "contig": CONTIG, "contig_m5": CONTIG_M5, "cram_url": cram_url,
        "reference": str(fasta), "known_sites": str(known),
        "fastq": [str(r1), str(r2)],
        "fastq_bytes": [r1.stat().st_size, r2.stat().st_size],
    }
    (gwaslib.ROOT_RESULTS / "fastq_testdata.json").write_text(json.dumps(manifest, indent=2))
    gwaslib.LOG.info("R1 %.1f MB, R2 %.1f MB", r1.stat().st_size / 1e6, r2.stat().st_size / 1e6)
    gwaslib.LOG.info("manifest written to results/fastq_testdata.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
