# Project brief: two human GWAS (pancreatic cancer and Alzheimer's disease)

## Goal
Set up two human genome-wide association studies (GWAS), one for pancreatic cancer and one for Alzheimer's disease (AD). Run every association and genotype-QC step in PLINK (1.9 and/or 2.0), following the GWASTutorial (https://cloufield.github.io/GWASTutorial).

## Fixed conventions
- Reference genome: GRCh38 (use the no-ALT analysis set, e.g. `GCA_000001405.15_GRCh38_no_alt_analysis_set`). If a dataset is only available on GRCh37/hg19, flag this and propose a liftover (CrossMap or `bcftools +liftover`) rather than mixing builds.
- Environment: uv with Python 3.12. Wrap every shell command in Python (`subprocess.run(..., check=True)`, with logging). Command-line tools that pip/uv cannot install (FastQC, fastp, BWA-MEM2, samtools, bcftools, GATK, PLINK) go into one separate conda/mamba env (bioconda) or are installed as pinned static binaries. Record every tool version.
- Hardware: local Linux machine, 20 cores, 63 GB RAM, no GPU in the sandbox. Size threads and memory to fit.

Installation using conda

conda create -n genome_env \
    --override-channels \
    -c conda-forge -c bioconda \
    --strict-channel-priority \
    fastqc fastp samtools bcftools bwa-mem2 gatk4 plink


conda activate genome_env

(genome_env) flavio@flavio:~$ conda list -n genome_env '^(fastqc|fastp|samtools|bcftools|bwa-mem2|gatk4.*|plink2?)$'
# packages in environment at /home/flavio/miniforge3/envs/genome_env:
#
# Name                     Version          Build            Channel
bcftools                   1.24             h118bc1c_2       bioconda
bwa-mem2                   2.3              he70b90d_0       bioconda
fastp                      1.3.7            h43da1c4_0       bioconda
fastqc                     0.12.1           hdfd78af_0       bioconda
gatk4                      4.6.2.0          py310hdfd78af_1  bioconda
plink                      1.90b7.7         h18e278d_1       bioconda
samtools                   1.24             h9dcdb79_1       bioconda


## Phase 1: Find data (search only, no downloads)
Note: a GWAS needs a cohort of cases and controls (usually hundreds to thousands of people). "One dataset" here means one cohort per disease.

For each disease, find candidate datasets and report them in a table with these columns: accession, repository, data type (genotyping array / WES / WGS FASTQ / VCF / summary statistics), number of cases and controls, ancestry, genome build, access level (open vs. controlled, e.g. dbGaP/EGA/ADNI/NIAGADS), approximate size, and link.

1. Pancreatic cancer. Prefer low-grade disease and state which definition each dataset uses: tumour grade G1 (well-differentiated PDAC), low-grade precursor lesions (PanIN, low-grade IPMN), or low-grade pancreatic neuroendocrine tumour (PanNET G1). If no grade-stratified data exist, fall back to all-grade PDAC and say so.
2. Alzheimer's disease. Note whether cases are clinical, autopsy-confirmed, or proxy (family history), and whether the dataset is late-onset or early-onset AD.

Rank sources in this order:
1. Open individual-level raw data (FASTQ/BAM/CRAM).
2. Open individual-level genotypes (VCF or PLINK files).
3. Controlled-access individual-level data. List the application steps; do not apply.
4. Open GWAS summary statistics (e.g. GWAS Catalog). These can support clumping, look-ups and comparison, but not a de novo PLINK association test.

Stop condition: if neither disease has openly downloadable individual-level data, report what you found (diseases, grades/subtypes, access level) and stop. Then propose options, such as applying for controlled access, a summary-statistics-only analysis, or a practice run on the tutorial's 1000 Genomes data with a simulated phenotype. Wait for my decision.

## Phase 2: Download plan (I run the downloads)
- Do not download any file larger than about 100 MB. Instead, give me the `wget`/`curl` (or `prefetch`/`fasterq-dump`, `aspera`, EGA client) commands, plus file sizes, MD5/SHA checksums, and the total disk space needed.
- Use this target layout: `<project_root>/data/{raw,reference,interim,results}`. I will download into it and grant you access to `<project_root>`.
- For controlled-access data I will handle credentials myself. Never ask me to paste keys or tokens into the chat.
- Also list the GRCh38 reference files needed (FASTA, .fai, .dict, BWA-MEM2 index, dbSNP, known-indel and HapMap/1000G resource VCFs for BQSR/VQSR).

## Phase 3: Code to write (after my approval of Phases 1–2)
Write modular, re-runnable Python scripts (per-sample, idempotent, logged), as follows:

1. Read QC on FASTQ files: FastQC per file, aggregated with MultiQC. Optionally add adapter/quality trimming with fastp, then run QC again. Outputs: the MultiQC HTML report plus a per-sample summary table (read count, %Q30, GC content, adapter content, duplication).
2. FASTQ to VCF on GRCh38:
   1. Align with BWA-MEM2 (with read groups).
   2. Sort and index with samtools, and mark duplicates with GATK MarkDuplicates.
   3. Run BQSR.
   4. Call per-sample gVCFs with GATK HaplotypeCaller.
   5. Joint-genotype with GenomicsDBImport and GenotypeGVCFs.
   6. Filter variants with VQSR, or hard filters for small cohorts.
   7. Normalise with `bcftools norm` (split multiallelics, left-align).
   8. Outputs: final bgzipped and tabix-indexed VCF, plus alignment QC (samtools flagstat/stats, mosdepth coverage) and variant QC (`bcftools stats`, Ti/Tv).
3. VCF to PLINK, with QC following the tutorial: convert (`plink2 --vcf`); sample and variant missingness; MAF and HWE filters; heterozygosity; sex check; relatedness (KING); LD pruning; PCA for population structure.
   Note: if the chosen dataset already contains genotypes (array or VCF), skip steps 1–2 and start here.

Out of scope until I confirm: imputation, the association test itself, and downstream steps (Manhattan/QQ plots, clumping, annotation).

## Working style
Wait for my approval at the end of each phase. State every assumption explicitly, and cite the accession for every dataset you propose.
# Project brief: two human GWAS (pancreatic cancer and Alzheimer's disease)

## Goal
Set up two human genome-wide association studies (GWAS), one for pancreatic cancer and one for Alzheimer's disease (AD). Run every association and genotype-QC step in PLINK (1.9 and/or 2.0), following the GWASTutorial (https://cloufield.github.io/GWASTutorial).

## Fixed conventions
- Reference genome: GRCh38 (use the no-ALT analysis set, e.g. `GCA_000001405.15_GRCh38_no_alt_analysis_set`). If a dataset is only available on GRCh37/hg19, flag this and propose a liftover (CrossMap or `bcftools +liftover`) rather than mixing builds.
- Environment: uv with Python 3.12. Wrap every shell command in Python (`subprocess.run(..., check=True)`, with logging). Command-line tools that pip/uv cannot install (FastQC, fastp, BWA-MEM2, samtools, bcftools, GATK, PLINK) go into one separate conda/mamba env (bioconda) or are installed as pinned static binaries. Record every tool version.
- Hardware: local Linux machine, 20 cores, 63 GB RAM, no GPU in the sandbox. Size threads and memory to fit.

Installation using conda

conda create -n genome_env \
    --override-channels \
    -c conda-forge -c bioconda \
    --strict-channel-priority \
    fastqc fastp samtools bcftools bwa-mem2 gatk4 plink


conda activate genome_env

(genome_env) flavio@flavio:~$ conda list -n genome_env '^(fastqc|fastp|samtools|bcftools|bwa-mem2|gatk4.*|plink2?)$'
# packages in environment at /home/flavio/miniforge3/envs/genome_env:
#
# Name                     Version          Build            Channel
bcftools                   1.24             h118bc1c_2       bioconda
bwa-mem2                   2.3              he70b90d_0       bioconda
fastp                      1.3.7            h43da1c4_0       bioconda
fastqc                     0.12.1           hdfd78af_0       bioconda
gatk4                      4.6.2.0          py310hdfd78af_1  bioconda
plink                      1.90b7.7         h18e278d_1       bioconda
samtools                   1.24             h9dcdb79_1       bioconda


## Phase 1: Find data (search only, no downloads)
Note: a GWAS needs a cohort of cases and controls (usually hundreds to thousands of people). "One dataset" here means one cohort per disease.

For each disease, find candidate datasets and report them in a table with these columns: accession, repository, data type (genotyping array / WES / WGS FASTQ / VCF / summary statistics), number of cases and controls, ancestry, genome build, access level (open vs. controlled, e.g. dbGaP/EGA/ADNI/NIAGADS), approximate size, and link.

1. Pancreatic cancer. Prefer low-grade disease and state which definition each dataset uses: tumour grade G1 (well-differentiated PDAC), low-grade precursor lesions (PanIN, low-grade IPMN), or low-grade pancreatic neuroendocrine tumour (PanNET G1). If no grade-stratified data exist, fall back to all-grade PDAC and say so.
2. Alzheimer's disease. Note whether cases are clinical, autopsy-confirmed, or proxy (family history), and whether the dataset is late-onset or early-onset AD.

Rank sources in this order:
1. Open individual-level raw data (FASTQ/BAM/CRAM).
2. Open individual-level genotypes (VCF or PLINK files).
3. Controlled-access individual-level data. List the application steps; do not apply.
4. Open GWAS summary statistics (e.g. GWAS Catalog). These can support clumping, look-ups and comparison, but not a de novo PLINK association test.

Stop condition: if neither disease has openly downloadable individual-level data, report what you found (diseases, grades/subtypes, access level) and stop. Then propose options, such as applying for controlled access, a summary-statistics-only analysis, or a practice run on the tutorial's 1000 Genomes data with a simulated phenotype. Wait for my decision.

## Phase 2: Download plan (I run the downloads)
- Do not download any file larger than about 100 MB. Instead, give me the `wget`/`curl` (or `prefetch`/`fasterq-dump`, `aspera`, EGA client) commands, plus file sizes, MD5/SHA checksums, and the total disk space needed.
- Use this target layout: `<project_root>/data/{raw,reference,interim,results}`. I will download into it and grant you access to `<project_root>`.
- For controlled-access data I will handle credentials myself. Never ask me to paste keys or tokens into the chat.
- Also list the GRCh38 reference files needed (FASTA, .fai, .dict, BWA-MEM2 index, dbSNP, known-indel and HapMap/1000G resource VCFs for BQSR/VQSR).

## Phase 3: Code to write (after my approval of Phases 1–2)
Write modular, re-runnable Python scripts (per-sample, idempotent, logged), as follows:

1. Read QC on FASTQ files: FastQC per file, aggregated with MultiQC. Optionally add adapter/quality trimming with fastp, then run QC again. Outputs: the MultiQC HTML report plus a per-sample summary table (read count, %Q30, GC content, adapter content, duplication).
2. FASTQ to VCF on GRCh38:
   1. Align with BWA-MEM2 (with read groups).
   2. Sort and index with samtools, and mark duplicates with GATK MarkDuplicates.
   3. Run BQSR.
   4. Call per-sample gVCFs with GATK HaplotypeCaller.
   5. Joint-genotype with GenomicsDBImport and GenotypeGVCFs.
   6. Filter variants with VQSR, or hard filters for small cohorts.
   7. Normalise with `bcftools norm` (split multiallelics, left-align).
   8. Outputs: final bgzipped and tabix-indexed VCF, plus alignment QC (samtools flagstat/stats, mosdepth coverage) and variant QC (`bcftools stats`, Ti/Tv).
3. VCF to PLINK, with QC following the tutorial: convert (`plink2 --vcf`); sample and variant missingness; MAF and HWE filters; heterozygosity; sex check; relatedness (KING); LD pruning; PCA for population structure.
   Note: if the chosen dataset already contains genotypes (array or VCF), skip steps 1–2 and start here.

Out of scope until I confirm: imputation, the association test itself, and downstream steps (Manhattan/QQ plots, clumping, annotation).

## Working style
Wait for my approval at the end of each phase. State every assumption explicitly, and cite the accession for every dataset you propose.
