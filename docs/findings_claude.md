The stop condition is met: I found no openly downloadable individual-level case/control data for either disease, so I've stopped at Phase 1. Every cohort that could support a PLINK association test is controlled-access, and none of them is stratified by tumour grade.

## Pancreatic cancer
- **PanScan (dbGaP phs000206)** is a GWAS of PDAC susceptibility. The project's genotype and phenotype data cover 9,437 individuals (PanScan I and II cases and controls, plus PanScan III cases only). Access is through a dbGaP request, with a Data Use Certification (DUC) Agreement.
- **PanC4 (phs000648)** brings together over 8,000 participants from 9 case-control studies in the United States, Europe and Australia. It is also on dbGaP and also controlled.
- **Grade:** these are the largest European-ancestry pancreatic cancer GWAS, and both cover all-grade PDAC.
- **Closest match to "low grade": IPMN germline GWAS.**
  - One is the first GWAS of IPMN progression, with a polygenic hazard score, in 338 patients.
  - The other used the Mass General Brigham Biobank: of 68,931 individuals, 2,525 had IPMNs.
  - Neither cohort is openly deposited.
- **PanNET:** the only germline studies I found are a 2016 neuroendocrine tumour GWAS that was mostly small-intestine tumours, and a candidate-SNP study of 320 PanNET cases and 4,436 controls. Neither provides grade-stratified PanNET data.

## Alzheimer's disease (late-onset)
- **NIAGADS** holds the main AD cohorts. Individual-level GWAS genotypes and ADSP whole genomes require a Data Access Request under the NIH Genomic Data Sharing Policy. Examples:
  - ADGC sets such as NG00022 (ADC1, 2,768 samples) and NG00023 (ADC2, 925 samples), both genotyping arrays.
  - ADSP Release 4 WGS: 36,361 participants, 45% of non-European ancestry.
- **ADNI** is also controlled. It contains 819 subjects: 229 cognitively normal, 398 with MCI and 192 with AD.
- **Open data are summary statistics only.** NIAGADS Open Access provides knowledge bases and downloadable annotation data, and the GWAS Catalog has large AD meta-analyses. These cannot be used for a new PLINK association test.

The survey table also lists build and phenotype for each source. Some accessions in it came from prior knowledge rather than this search, such as the Kunkle 2019 AD summary statistics (NG00075) and Bellenguez 2022 (GCST90027158). They are marked "(verify)" and need confirming before any download.

## Options
- **A. Apply for controlled access:** dbGaP for PanScan and PanC4; NIAGADS for ADSP and ADGC; the ADNI application for ADNI. Expect weeks to months.
- **B. Practice run now:** the GWASTutorial 1000 Genomes data with a simulated phenotype, through the full PLINK QC and association workflow. For FASTQ → VCF practice, add a few open 1000 Genomes high-coverage samples (GRCh38).
- **C. Summary-statistics track:** open AD and pancreatic summary statistics, with PLINK clumping against a 1000 Genomes LD reference, look-ups and polygenic risk scores.

B and C can run while A is pending.

## Gaps in your conda environment (genome_env)
It has PLINK 1.9 only. The brief also needs `plink2` (for `--vcf` conversion and KING relatedness), `multiqc` and `mosdepth`, so these should be added before Phase 3.

Which option or options do you want to pursue?

[phase1_dataset_survey.md]({{artifact:002195b5-e5ef-43f4-9881-ea9a3324ebd6}}) has the full survey table and the options.




### if we choose GWASTutorial 1000 Genome

how much space in disc must I have which diseases it covers?

## Disk space

Sizes below are from the GitHub API and the tutorial's own file listings.

| Item | Size |
|---|---|
| Repository clone (`Cloufield/GWASTutorial`) | ~458 MB (468,559 KB reported by GitHub) |
| → of which the sample dataset zip `1KG.EAS.auto...missing.zip` | 63.8 MB (63,756,897 bytes), shipped inside the repo |
| Unzipped PLINK binary set | ~190 MB: 149M .bed, 40M .bim, 13K .fam |
| Tutorial intermediate outputs (QC-cleaned set, PCA, GRM, LD matrices, association results) | ~300 MB |
| **Core PLINK track total** | **~1 GB** |

So about **1 GB**, and nothing needs downloading separately: the sample dataset zip is included in `01_Dataset` when you clone the repository, so there is no need to download it again. You have ~46 GB free, which is ample.

Two caveats on the optional modules:
- The visualization module uses gwaslab, whose reference panel is larger: the processed 1000 Genomes EAS dataset (hg19) is about 2.8 GB. Its own sample summary statistics are small (~40 MB).
- Later modules (LD score regression, functional annotation with VEP/ANNOVAR, PRS) need extra reference databases I have not sized. Those can add several GB each; I'd check them only if you want to run those modules.

## Which diseases it covers

**None — the phenotype is simulated, not a real disease.** The genotypes are 504 East Asian individuals from the 1000 Genomes Project, a healthy population reference with no clinical labels, carrying 1,235,116 variants. Phenotypes were simulated with GCTA from the 1KG EAS genotypes, using `--simu-cc 250 254 --simu-hsq 0.8 --simu-k 0.5`: 250 cases and 254 controls, 80% heritability, 50% prevalence, with causal SNPs picked from a supplied `causal.snplist`. In the phenotype file, `B1` is the phenotype name, where 1 means control and 2 means case.

A 50% prevalence and 0.8 heritability are deliberately unrealistic — they make the signal easy to detect for teaching. The trait is labelled "B1", not any named condition.

The only real disease data in the tutorial is summary statistics: the gwaslab module uses processed type 2 diabetes summary statistics from Biobank Japan, restricted to HapMap3 variants on chromosome 7.

## Two consequences for our project

1. **Build mismatch.** The tutorial data is GRCh37, listed among its own key terms as "1000 Genomes Project (1KG), EAS, GRCh37". Our brief fixes GRCh38. For a practice run this does not matter, as long as we keep the tracks separate and never mix the two builds.
2. **No FASTQ files.** The tutorial starts from PLINK binaries, so Phase 3 steps 1 and 2 (FASTQ quality metrics and FASTQ → VCF) cannot be practiced on it. For those we need separate open data, for example a 1000 Genomes high-coverage sample from ENA. Whole-genome FASTQ runs 40–120 GB per sample, so I would restrict it to one chromosome or a downsampled subset to keep it within your free space.

If you want option B, tell me whether to include a small FASTQ exercise alongside the PLINK track, and I'll size the download before you fetch anything.

