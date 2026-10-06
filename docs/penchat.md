**Yes — PennCATH is openly downloadable, with no application.** It is the only real-patient, individual-level GWAS dataset I've found that is genuinely open, so it meets rank 2 of your priority list. The catch is that the disease is CAD, not pancreatic cancer or Alzheimer's.

I verified everything below against the live archive using HTTP range requests, reading only the zip directory, the 7.6 MB `.bim` and the two tiny text files. I did not download the 302 MB `.bed`.

## What it is
| Field | Value |
|---|---|
| Study | PennCATH (Penn Catheterization study), CAD genetic risk factors — PMID 21239051 |
| Tutorial | Reed et al. 2015 (doi:10.1002/sim.6605); maintained version by Breheny & Peter |
| Access | Open, direct HTTP, no credentials |
| Archive | 142.8 MB (`penncath.tar.gz` or `.zip`, identical), last modified 2023-05-16 |
| Samples | **1,401 individuals** |
| Variants | **861,473** autosomal |
| Cases / controls | **933 CAD cases / 468 controls** |
| Sex | 937 male / 464 female |
| Age | mean 55.7, range 22–87 |
| Covariates | age, sex, triglycerides, HDL, LDL (92–119 missing values in the lipid columns) |
| Build | **GRCh37/hg19** |

The variant count is self-consistent: the `.bed` is 302,377,026 bytes, which is exactly 3 + 351 × 861,473, and the `.bim` has 861,473 lines.

The build is confirmed rather than assumed. rs3094315 sits at 1:752566 in the `.bim`, which is its GRCh37 position; in GRCh38 it is at 1:817186. Two further SNPs gave the same result.

## Disk space
| Item | Size |
|---|---|
| Downloaded archive | 143 MB |
| Unzipped: bed 302.38 MB + bim 24.05 MB + csv 38 KB + fam 29 KB | 327 MB |
| QC, PCA and association outputs | ~400 MB |
| **Subtotal** | **~900 MB** |
| Optional GRCh38 liftover (second copy + intermediate VCF) | ~700 MB |

Budget **1.6 GB** if you want the liftover, 900 MB without. Either is trivial against your 49 GB free.

```bash
cd <project_root>/data/raw
wget -c https://d1ypx1ckp5bo16.cloudfront.net/penncath/penncath.tar.gz
tar -xzf penncath.tar.gz     # creates data/penncath.{bed,bim,csv,fam}
```
No published checksum exists for this file, so I'd compute an `sha256sum` after download and record it for reproducibility.

## Four things to handle in QC
1. **The `.fam` phenotype column is miscoded and must not be used.** Controls are coded 0 (468 individuals) and cases are split across 1 (470) and 2 (463). Under PLINK's convention — 1 control, 2 case, 0 missing — that would invert the cases and discard the controls. The phenotype has to come from `penncath.csv` via `--pheno`. The 1-versus-2 split nests perfectly inside CAD=1, so it is probably an unlabelled case subtype.
2. **26,154 variants carry a `-9` missing allele code.** This is exactly what makes the tutorial's first `--maf` call fail; the files need rewriting with `--make-bed` first.
3. **No sex chromosomes** (chr1–22 only), so `--check-sex` is impossible. The only sex check available is `.fam` against the csv, and they agree for all 1,401 samples.
4. **GRCh37.** Your brief fixes GRCh38, so either lift over (convert to VCF, `CrossMap`/`Picard LiftoverVcf` with the hg19→hg38 chain, then back to PLINK, losing a small share of variants) or record GRCh37 as a declared exception for this dataset. For QC and association the build is irrelevant; it matters only when annotating hits or merging with other data.

## How this compares with the options
PennCATH is a better practice vehicle than the 1000 Genomes track in one important way: the phenotype is real. The tutorial's simulated trait has 80% heritability and 50% prevalence, so the association step is unrealistically easy. PennCATH gives real effect sizes, real confounding and a real case/control imbalance. A 933/468 study is small by modern standards and will likely recover only the 9p21 locus, which is the honest outcome for a cohort this size.

It does not, however, contain FASTQ files, so Phase 3 steps 1 and 2 still need separate data.

Three ways forward, and I'd suggest the first:
- **Run the full PLINK workflow on PennCATH** (real CAD data), and separately fetch one small open FASTQ sample for the quality-metrics and VCF-calling code.
- Keep PennCATH as practice and still apply for dbGaP/NIAGADS access for the two diseases you actually want.
- Drop CAD and go with the summary-statistics track for pancreatic cancer and AD.

Which do you want? If it's the first, tell me whether to include the liftover to GRCh38 or to stay on GRCh37.

[phase1_dataset_survey.md]({{artifact:bbc02b56-af16-4b64-a2d0-7a13d5891e83}}) — now version 2, with the full PennCATH fact sheet added.


# Phase 1: Dataset survey (search only, nothing downloaded)

Date: 2026-10-02. Build is GRCh37/hg19 unless stated otherwise. Items marked "(verify)" are from prior knowledge and were not confirmed by a search in this session.

## Pancreatic cancer

| Accession | Repository | Data type | Cases / controls | Ancestry | Build | Access | Grade / subtype |
|---|---|---|---|---|---|---|---|
| phs000206 (PanScan I–III) | dbGaP | Illumina SNP arrays (HumanHap550, 610-Quad, OmniExpress) | PanScan I+II cases and controls + PanScan III cases, 9,437 subjects | European | hg19 | Controlled (dbGaP DAR) | All-grade PDAC, not stratified by grade |
| phs000648 (PanC4) | dbGaP | Illumina OmniExpressExome array | >8,000 participants, 9 case-control studies | European (US, Europe, Australia) | hg19 | Controlled (dbGaP DAR) | All-grade PDAC |
| PanScan I–III summary statistics (Wolpin 2014; Klein 2018) | dbGaP (inside phs000206) | Summary statistics | 5,117 / 8,845 (PanScan I–III) | European | hg19 | Controlled within dbGaP (verify whether open) | All-grade PDAC |
| IPMN progression GWAS (PMID 39639588) | Not deposited openly (verify) | Array | 338 IPMN patients (Cox model) | Not checked | Not checked | Not public | Includes low-grade IPMN (the closest match to "low grade") |
| IPMN GWAS, Mass General Brigham (MGB) Biobank (PMID 41060111) | MGB Biobank | Array | 2,525 IPMN / ~66,400 non-IPMN | Mostly European | Not checked | Institutional only | IPMN of any grade |
| NET GWAS (Endocr Relat Cancer 2016) | Not deposited openly (verify) | Affymetrix 6.0 | Mostly small-intestine NET; few pancreatic NET (PanNET) cases | European | hg19 | Not public | No grade-stratified PanNET data |
| FinnGen pancreatic cancer endpoint (verify) | FinnGen | Summary statistics | Several hundred to ~1,000 cases (verify) | Finnish | GRCh38 | Open download after registration (verify) | All-grade |

## Alzheimer's disease

| Accession | Repository | Data type | Cases / controls | Ancestry | Build | Access | Phenotype |
|---|---|---|---|---|---|---|---|
| ADNI | LONI IDA / NIAGADS | Illumina arrays + WGS | 819 subjects (229 cognitively normal, 398 MCI, 192 AD) in the ADNI-1 genotyped set | Mostly European | hg19 (arrays), GRCh38 (WGS) | Controlled (ADNI application) | Clinical diagnosis, late-onset AD (LOAD) |
| NG00022 / NG00023 (ADGC ADC1 / ADC2) | NIAGADS DSS | SNP arrays | 2,768 / 925 samples | Mostly European | hg19 | Controlled (NIAGADS DAR) | Clinical LOAD, part autopsy-confirmed |
| ADSP Release 4 WGS | NIAGADS DSS | WGS (CRAM / pVCF) | 36,361 participants, 45% non-European | Multi-ancestry | GRCh38 | Controlled (NIAGADS DAR) | Clinical + autopsy, LOAD |
| NG00075 (Kunkle 2019, IGAP) (verify) | NIAGADS Open Access | Summary statistics | 21,982 / 41,944 (stage 1) | European | hg19 | Open | Clinical + autopsy LOAD |
| GCST90027158 (Bellenguez 2022) (verify) | EBI GWAS Catalog | Summary statistics | ~111k cases (incl. proxy) / ~677k controls | European | GRCh38 | Open | Clinical + proxy (family history) |

## Addendum: PennCATH (coronary artery disease) — open individual-level genotypes

Checked 2026-10-02 by HTTP range requests against the live archive (no full download). This is a different disease from the two in the brief, but it is the only open, individual-level, real-patient GWAS dataset found so far, and it meets rank 2 of the priority list.

| Field | Value |
|---|---|
| Study | PennCATH (Penn Catheterization study), CAD genetic risk factors; PMID 21239051 |
| Tutorial | Reed et al. 2015, Stat Med, doi:10.1002/sim.6605; maintained version by Breheny & Peter (`pbreheny/adv-gwas-tutorial`) |
| Access | Open, no application, direct HTTP |
| URL | `https://d1ypx1ckp5bo16.cloudfront.net/penncath/penncath.tar.gz` (or `.zip`, identical contents) |
| Archive size | 142,819,504 bytes (tar.gz) / 142,818,657 bytes (zip); last modified 2023-05-16 |
| Data type | Genotyping array (PLINK bed/bim/fam), autosomes only |
| Samples | 1,401 individuals |
| Variants | 861,473 (verified: .bed is 302,377,026 bytes = 3 + 351 × 861,473; .bim has 861,473 lines) |
| Chromosomes | 1–22 only; no X, Y or MT |
| Phenotype | `penncath.csv`: FamID, CAD, sex, age, tg, hdl, ldl |
| Cases / controls | 933 CAD cases / 468 controls |
| Sex | 937 male / 464 female (.fam SEX agrees with the csv for all 1,401) |
| Age | mean 55.7, range 22–87 |
| Build | **GRCh37/hg19** — verified against Ensembl: rs3094315 .bim 1:752566 = GRCh37 1:752566, GRCh38 1:817186; same for rs12565286 and rs2980300 |

Unzipped contents: `penncath.bed` 302.38 MB, `penncath.bim` 24.05 MB, `penncath.csv` 38 KB, `penncath.fam` 29 KB — 326.5 MB total.

Known issues to handle in QC:
- The `.fam` phenotype column is miscoded for PLINK. Controls are 0 (468) and cases are split across 1 (470) and 2 (463); nesting within CAD is perfect. Read under PLINK's 1=control / 2=case / 0=missing convention this inverts and drops the phenotype, so the phenotype must be supplied from `penncath.csv` with `--pheno`.
- 26,154 variants carry a `-9` missing allele code, which is what makes the tutorial's first `--maf` call fail until the files are rewritten with `--make-bed`.
- No sex chromosomes, so `--check-sex` is not possible; sex can only be compared between the `.fam` and the csv (they agree).
- GRCh37, so it needs a liftover to meet the project's GRCh38 convention, or the build has to be declared as an exception.

## Result for the stop condition

The stop condition is met. No openly downloadable individual-level data (FASTQ/BAM/CRAM, VCF or PLINK) with case/control labels were found for either disease. Every cohort suitable for a PLINK association test is controlled-access (dbGaP or NIAGADS DAR, or an ADNI application). No grade-stratified pancreatic cancer GWAS cohort exists in public repositories. The closest match to "low grade" is IPMN germline GWAS, but those cohorts are small and not deposited openly. Open data exist only as summary statistics, which cannot be used for a de novo PLINK association test.

## Options for the decision
- A. Apply for controlled access: PanScan and PanC4 via a dbGaP DAR (eRA Commons account, institutional signing official, Data Use Certification, IRB statement), and ADNI and/or ADSP via a NIAGADS DSS DAR. Typical turnaround is weeks to months.
- B. Practice pipeline now: run the GWASTutorial 1000 Genomes data (open, about 500 samples) with a simulated case/control phenotype, through the full PLINK QC and association workflow. For FASTQ → VCF practice, use a few open 1000 Genomes high-coverage samples (ENA/SRA), lifted over or already on GRCh38.
- C. Summary-statistics track: download the open AD statistics (Bellenguez 2022 GCST90027158 / Kunkle NG00075) and an open pancreatic set (FinnGen). Then do PLINK `--clump` with the 1000 Genomes LD reference, look-ups, cross-disease comparison, and a polygenic risk score (PRS) applied to the 1000 Genomes samples.
- B and C can run while A is pending.



