# PennCATH GWAS — PLINK track

Genotype QC pipeline for the PennCATH coronary artery disease cohort
(1,401 samples × 861,473 autosomal variants, GRCh37/hg19).

## Layout

    scripts/gwaslib.py              config, logging, subprocess wrapper, PLINK helpers
    scripts/01_setup.py             create directories, record tool versions
    scripts/02_download_penncath.py print / run / verify the 143 MB download
    scripts/03_qc_penncath.py       the QC pipeline
    data/PennCATH/raw/              penncath.{bed,bim,fam,csv}
    data/PennCATH/interim/          one PLINK fileset per QC stage
    data/PennCATH/results/          final fileset, PCA, qc_summary.csv
    data/PennCATH/figures/          penncath_qc.png
    logs/                           one log per script run

Set `GWAS_DATASET=Alzheimer` to point the same scripts at
`data/Alzheimer/`; `GWAS_ROOT` relocates the whole tree.

## Run

    conda activate genome_env          # needs plink and plink2
    python 01_setup.py
    python 02_download_penncath.py     # prints the wget command
    python 02_download_penncath.py --run   # or download yourself, then this verifies
    python 03_qc_penncath.py

`plink2` is not in `genome_env` yet:

    conda install -n genome_env -c conda-forge -c bioconda plink2

Paths and thresholds are overridable by environment variable, so nothing is
hard-coded to one machine:

    GWAS_ROOT=/path/to/project PLINK=/path/to/plink PLINK2=/path/to/plink2 python 03_qc_penncath.py

## QC stages

| Stage | Operation | Why |
|---|---|---|
| 00 | read raw fileset | baseline counts |
| 00 | build `.pheno` / `.covar` from `penncath.csv` | the `.fam` phenotype column is miscoded |
| 01a | `--update-ids` | strips the literal quotes from sample IDs; also rewrites the 26,154 `-9` allele codes that block `--maf` |
| 01b | `--pheno`, then verify | a non-matching ID file makes `--pheno` a silent no-op |
| 02a | `--geno 0.1` | drop poorly genotyped variants |
| 02b | `--mind 0.1` | drop poorly genotyped samples (after 02a, so fewer are lost) |
| 03a | `--maf 0.01` | rare variants have no power at n=1,401 |
| 03b | `--hwe 1e-10 midp` | genotyping error; controls only, by PLINK default for case/control |
| — | `--indep-pairwise 50 5 0.2` | LD-pruned subset for the steps that assume independence |
| 04 | `--het`, drop \|F − mean\| > 3 SD | contamination or poor DNA |
| 05 | `plink2 --king-cutoff 0.0884` | remove one of each pair closer than 2nd-degree |
| — | `plink2 --pca 10` | population-structure outliers and association covariates |

Final output `data/results/penncath_qc.{bed,bim,fam}` keeps all QC-passing
variants, not only the pruned subset — pruning exists for the QC steps above,
not for the association test.

## Known properties of this dataset

* **Sample IDs contain literal double quotes.** The `.fam` first column is
  `"10002"`, quotes included, for all 1,401 samples, and PLINK treats them as
  part of the ID. Any ID file written with the quotes stripped matches nothing,
  and `--pheno` then fails *silently* — PLINK reports `0 phenotype values
  present` in its log and carries on, so `--hwe` is applied to all samples
  instead of controls and every later step runs phenotype-blind. Stage 01a
  strips the quotes once with `--update-ids`; stage 01b attaches the phenotype
  and `assert_phenotype()` then verifies the counts in the `.fam`.
* `.fam` phenotype column: controls 0, cases split across 1 and 2. Unusable
  directly; the phenotype comes from `penncath.csv` (933 cases / 468 controls).
* No X or Y chromosome, so `--check-sex` is impossible. `01` cross-checks the
  `.fam` sex column against the clinical file instead (they agree for all 1,401).
* Build GRCh37/hg19. Liftover to GRCh38 is a separate step, not needed for QC
  or association.
* Lipid covariates have missing values (tg 92, hdl 92, ldl 119), written as `NA`.

## Not included yet

Association testing, Manhattan/QQ plots, clumping, annotation, imputation, and
the FASTQ → VCF track.
