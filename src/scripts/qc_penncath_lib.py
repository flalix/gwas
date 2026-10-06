#!/usr/bin/env python3
"""Step 3 - genotype quality control for PennCATH, in PLINK.

    python 03_qc_penncath.py

Pipeline, each stage written to data/interim/ so it can be inspected or resumed:

    00  phenotype + covariate files built from penncath.csv
    01  normalise the fileset and attach the phenotype
    02  call-rate filters   (--geno, --mind)
    03  MAF and HWE filters (--maf, --hwe)
    --  LD pruning          (--indep-pairwise) -> SNP set for the steps below
    04  heterozygosity outliers removed  (|F - mean| > 3 SD)
    05  relatedness         (plink2 --king-cutoff)
    --  PCA                 (plink2 --pca)

Two PennCATH-specific points, both established by inspecting the files:

* The .fam phenotype column is miscoded.  Controls are 0 and cases are split
  across 1 and 2.  Read with PLINK's 1=control/2=case/0=missing convention that
  would invert the cases and drop every control, so the phenotype is taken from
  penncath.csv instead.
* 26,154 variants carry a '-9' missing allele code, which makes PLINK refuse
  --maf until the fileset has been rewritten.  Stage 01 does that rewrite.

There is no X chromosome in this dataset, so --check-sex is impossible; stage 01
instead cross-checks the .fam sex against the clinical file.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import gwaslib

RAW_STEM = gwaslib.ROOT_RAW / "penncath"
CLINICAL = gwaslib.ROOT_RAW / "penncath.csv"

STAGES: list[tuple[str, int, int]] = []     # (label, n_variants, n_samples)


def record(label: str, stem) -> None:
    n_var, n_sam = gwaslib.count_bfile(stem)
    STAGES.append((label, n_var, n_sam))
    gwaslib.LOG.info("[%s] %d variants x %d samples", label, n_var, n_sam)


# --------------------------------------------------------------------------
# Stage 00 - phenotype and covariates
# --------------------------------------------------------------------------
def build_phenotype() -> tuple[gwaslib.Path, gwaslib.Path, gwaslib.Path, pd.DataFrame]:
    # quoting=QUOTE_NONE is essential: PennCATH ships FIDs wrapped in literal
    # double quotes ("10002"), and PLINK treats the quotes as part of the ID.
    # Let pandas strip them and every ID file we write silently matches nothing
    # -- which is exactly how the phenotype gets lost without an error.
    fam = pd.read_csv(f"{RAW_STEM}.fam", sep=r"\s+", header=None, dtype=str,
                      quoting=csv.QUOTE_NONE,
                      names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"])
    fam["FID_CLEAN"] = fam["FID"].str.strip('"')
    fam["IID_CLEAN"] = fam["IID"].str.strip('"')
    n_quoted = int((fam["FID"] != fam["FID_CLEAN"]).sum())
    gwaslib.LOG.info("%d of %d sample IDs carry literal quote characters", n_quoted, len(fam))

    # Mapping consumed by --update-ids: old FID, old IID, new FID, new IID.
    idmap_path = gwaslib.ROOT_INTERIM / "penncath.update_ids"
    fam[["FID", "IID", "FID_CLEAN", "IID_CLEAN"]].to_csv(
        idmap_path, sep="\t", index=False, header=False, quoting=csv.QUOTE_NONE)

    clinical = pd.read_csv(CLINICAL, dtype={"FamID": str})
    merged = fam.merge(clinical, left_on="FID_CLEAN", right_on="FamID",
                       how="inner", validate="one_to_one")
    if len(merged) != len(fam):
        raise ValueError(f"{len(fam)} samples in .fam but {len(merged)} matched penncath.csv")

    # PLINK case/control coding: 1 = control, 2 = case, 0/-9 = missing.
    merged["CAD_PLINK"] = merged["CAD"] + 1
    merged["SEX"] = merged["SEX"].astype(int)

    pheno_path = gwaslib.ROOT_INTERIM / "penncath.pheno"
    merged[["FID_CLEAN", "IID_CLEAN", "CAD_PLINK"]].to_csv(
        pheno_path, sep="\t", index=False, header=["FID", "IID", "CAD"],
        quoting=csv.QUOTE_NONE)

    covar_path = gwaslib.ROOT_INTERIM / "penncath.covar"
    covar = merged[["FID_CLEAN", "IID_CLEAN", "sex", "age", "tg", "hdl", "ldl"]].copy()
    covar.columns = ["FID", "IID", "sex", "age", "tg", "hdl", "ldl"]
    covar.to_csv(covar_path, sep="\t", index=False, na_rep="NA", quoting=csv.QUOTE_NONE)

    n_case = int((merged["CAD"] == 1).sum())
    n_ctrl = int((merged["CAD"] == 0).sum())
    gwaslib.LOG.info("phenotype: %d cases, %d controls", n_case, n_ctrl)

    # Sex check by proxy: the .fam column against the clinical file.
    disagree = int((merged["SEX"] != merged["sex"]).sum())
    gwaslib.LOG.info("sex .fam vs clinical: %d disagreements (no X chromosome, so "
               "--check-sex is not available)", disagree)

    missing_covar = covar[["age", "tg", "hdl", "ldl"]].isna().sum().to_dict()
    gwaslib.LOG.info("covariate missingness: %s", missing_covar)
    return pheno_path, covar_path, idmap_path, merged


def assert_phenotype(stem, n_case: int, n_ctrl: int) -> None:
    """A mismatched ID file makes --pheno a silent no-op, so verify it took."""
    fam = pd.read_csv(f"{stem}.fam", sep=r"\s+", header=None, dtype=str,
                      quoting=csv.QUOTE_NONE,
                      names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"])
    counts = fam["PHENO"].value_counts().to_dict()
    got_case, got_ctrl = int(counts.get("2", 0)), int(counts.get("1", 0))
    gwaslib.LOG.info("phenotype in .fam: %d cases, %d controls, %d missing",
               got_case, got_ctrl, len(fam) - got_case - got_ctrl)
    if (got_case, got_ctrl) != (n_case, n_ctrl):
        raise ValueError(
            f"--pheno did not take: .fam has {got_case}/{got_ctrl} "
            f"but penncath.csv has {n_case}/{n_ctrl}. Check sample-ID formatting.")



def write_summary() -> None:
    summary = pd.DataFrame(STAGES, columns=["stage", "variants", "samples"])
    summary["variants_lost"] = -summary["variants"].diff().fillna(0).astype(int)
    summary["samples_lost"] = -summary["samples"].diff().fillna(0).astype(int)
    path = gwaslib.ROOT_RESULTS / "qc_summary.csv"
    summary.to_csv(path, index=False)
    gwaslib.LOG.info("wrote %s\n%s", path, summary.to_string(index=False))


def make_figures(clinical: pd.DataFrame, pca_stem) -> None:
    kw = dict(sep=r"\s+", quoting=csv.QUOTE_NONE)
    imiss = pd.read_csv(gwaslib.ROOT_INTERIM / "rep_missing.imiss", **kw)
    lmiss = pd.read_csv(gwaslib.ROOT_INTERIM / "rep_missing.lmiss", **kw)
    freq = pd.read_csv(gwaslib.ROOT_INTERIM / "rep_freq.frq", **kw)
    hardy = pd.read_csv(gwaslib.ROOT_INTERIM / "rep_hardy.hwe", **kw)
    het = pd.read_csv(gwaslib.ROOT_INTERIM / "rep_het.het", **kw)
    eigenvec = pd.read_csv(f"{pca_stem}.eigenvec", **kw)
    eigenval = np.loadtxt(f"{pca_stem}.eigenval")

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))

    ax = axes[0, 0]
    ax.hist(imiss["F_MISS"], bins=50, color="#4C72B0")
    ax.axvline(gwaslib.QC["mind"], color="crimson", ls="--", lw=1)
    ax.set(xlabel="per-sample missing rate", ylabel="samples",
           title=f"Sample call rate (n={len(imiss)})")

    ax = axes[0, 1]
    ax.hist(lmiss["F_MISS"], bins=50, color="#4C72B0")
    ax.axvline(gwaslib.QC["geno"], color="crimson", ls="--", lw=1)
    ax.set(xlabel="per-variant missing rate", ylabel="variants", yscale="log",
           title=f"Variant call rate (n={len(lmiss):,})")

    ax = axes[0, 2]
    ax.hist(freq["MAF"], bins=50, color="#4C72B0")
    ax.axvline(gwaslib.QC["maf"], color="crimson", ls="--", lw=1)
    ax.set(xlabel="minor allele frequency", ylabel="variants",
           title="MAF spectrum")

    ax = axes[1, 0]
    merged = het.merge(imiss, on=["FID", "IID"])
    mean_f, sd_f = het["F"].mean(), het["F"].std()
    ax.scatter(merged["F_MISS"], merged["F"], s=10, alpha=0.6, color="#4C72B0")
    for sign in (-1, 1):
        ax.axhline(mean_f + sign * gwaslib.QC["het_sd"] * sd_f, color="crimson", ls="--", lw=1)
    ax.set(xlabel="per-sample missing rate", ylabel="heterozygosity F",
           title=f"Heterozygosity vs missingness (±{gwaslib.QC['het_sd']:.0f} SD)")

    ax = axes[1, 1]
    # PLINK writes one row per TEST: ALL / AFF / UNAFF when a case-control
    # phenotype is loaded, and a single ALL(NP) row when it is not.
    tests = set(hardy["TEST"].unique())
    which = "UNAFF" if "UNAFF" in tests else sorted(tests)[0]
    subset = hardy[hardy["TEST"] == which].dropna(subset=["P"])
    ax.hist(-np.log10(subset["P"].clip(lower=1e-300)), bins=50, color="#4C72B0")
    ax.axvline(-np.log10(gwaslib.QC["hwe"]), color="crimson", ls="--", lw=1)
    ax.set(xlabel=f"-log10 HWE p ({which})", ylabel="variants", yscale="log",
           title="Hardy-Weinberg equilibrium")

    ax = axes[1, 2]
    status = clinical.set_index("FID_CLEAN")["CAD"]
    pcs = eigenvec.rename(columns={"#FID": "FID", "#IID": "IID"})
    pcs["CAD"] = pcs["FID"].astype(str).map(status)
    for value, label, col in [(0, "control", "#4C72B0"), (1, "CAD case", "#C44E52")]:
        subset = pcs[pcs["CAD"] == value]
        ax.scatter(subset["PC1"], subset["PC2"], s=10, alpha=0.6, label=label, color=col)
    # Only 10 PCs were computed, so eigenval.sum() is not the total variance;
    # quoting a "% variance explained" from it would overstate each PC.  Report
    # the raw eigenvalues instead.
    ax.set(xlabel=f"PC1 (eigenvalue {eigenval[0]:.1f})",
           ylabel=f"PC2 (eigenvalue {eigenval[1]:.1f})",
           title=f"Population structure ({gwaslib.QC['n_pcs']} PCs computed)")
    ax.legend(frameon=False)

    fig.tight_layout()
    path = gwaslib.ROOT_FIGURES / "penncath_qc.png"
    fig.savefig(path, dpi=200)
    gwaslib.LOG.info("wrote %s", path)




