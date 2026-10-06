#!/usr/bin/env python3
"""Step 4 - case/control association testing for PennCATH.

Companion to qc_penncath_lib.py; consumes its output fileset
data/PennCATH/results/penncath_qc.{bed,bim,fam} (696,660 variants x 1,382
samples, 923 cases / 459 controls, build GRCh37/hg19).

    from scripts.assoc_penncath_lib import *

    covar  = build_assoc_covar(n_pcs=4)
    glm    = run_glm(covar)
    df     = load_glm(glm)
    lam    = lambda_gc(df["P"])
    clumps = run_clump(glm)
    genes  = annotate_loci(top_hits(df))

Model: logistic regression on allele dosage, adjusted for sex, age and the
first N principal components, with Firth fallback for sparse tables.

Covariate choice, deliberately:
* sex and age are confounders of CAD and belong in the model.
* PCs 1-4 only.  The QC eigenvalues were 6.24, 2.19, 1.85, 1.69 and then flat
  at ~1.35, so PC5 onward is noise and would only cost degrees of freedom.
* the lipid measurements (tg, hdl, ldl) are deliberately EXCLUDED.  They sit on
  the causal pathway between genotype and coronary disease, so conditioning on
  them would attenuate exactly the effects we are trying to detect, and they
  are missing for 92-119 samples, which would shrink the analysed set.

Build note: PennCATH coordinates are GRCh37/hg19, so gene annotation queries
the GRCh37 Ensembl endpoint.  Using the current GRCh38 endpoint would return
the wrong genes -- positions differ by tens of kilobases.
"""
from __future__ import annotations

import csv
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy import stats

import gwaslib

QC_STEM = gwaslib.ROOT_RESULTS / "penncath_qc"
PCA_STEM = gwaslib.ROOT_RESULTS / "penncath_pca"
COVAR_IN = gwaslib.ROOT_INTERIM / "penncath.covar"

GENOME_WIDE = 5e-8        # conventional threshold
SUGGESTIVE = 1e-5         # conventional "worth a look" threshold

# Ensembl GRCh37 endpoint -- PennCATH is hg19, NOT hg38.
ENSEMBL_GRCH37 = "https://grch37.rest.ensembl.org"


# --------------------------------------------------------------------------
# Covariates
# --------------------------------------------------------------------------
def build_assoc_covar(n_pcs: int = 4, out: Path | None = None) -> Path:
    """Merge the clinical covariates with the PCA eigenvectors."""
    out = out or gwaslib.ROOT_INTERIM / "penncath_assoc.covar"
    kw = dict(sep=r"\s+", quoting=csv.QUOTE_NONE)

    covar = pd.read_csv(COVAR_IN, dtype={"FID": str, "IID": str}, **kw)
    pcs = pd.read_csv(f"{PCA_STEM}.eigenvec", dtype={"#FID": str, "IID": str}, **kw)
    pcs = pcs.rename(columns={"#FID": "FID", "#IID": "IID"})

    pc_cols = [f"PC{i}" for i in range(1, n_pcs + 1)]
    missing = set(pc_cols) - set(pcs.columns)
    if missing:
        raise ValueError(f"{PCA_STEM}.eigenvec lacks {sorted(missing)}")

    merged = pcs[["FID", "IID", *pc_cols]].merge(
        covar[["FID", "IID", "sex", "age"]], on=["FID", "IID"],
        how="inner", validate="one_to_one")
    if len(merged) != len(pcs):
        raise ValueError(f"{len(pcs)} samples in PCA but {len(merged)} matched covariates")

    merged.to_csv(out, sep="\t", index=False, na_rep="NA", quoting=csv.QUOTE_NONE)
    gwaslib.LOG.info("covariates for %d samples: %s", len(merged),
                     ", ".join(["sex", "age", *pc_cols]))
    return out


# --------------------------------------------------------------------------
# Association
# --------------------------------------------------------------------------
def run_glm(covar: Path, bfile: Path = QC_STEM, out: Path | None = None,
            threads: int = 8) -> Path:
    """plink2 --glm logistic regression, Firth fallback, covariates hidden.

    Returns the path to the .PHENO1.glm.logistic.hybrid table.
    """
    out = out or gwaslib.ROOT_RESULTS / "penncath_assoc"
    gwaslib.plink([
        "--bfile", bfile,
        "--covar", covar,
        # put covariates on a common scale so the optimiser behaves; this
        # changes covariate coefficients, not the SNP effect we report
        "--covar-variance-standardize",
        "--glm", "firth-fallback", "hide-covar",
        "--ci", "0.95",
        "--threads", str(threads),
    ], out=out, v2=True)

    produced = sorted(Path(out).parent.glob(f"{Path(out).name}*.glm.logistic*"))
    if not produced:
        raise FileNotFoundError(f"no --glm output next to {out}")
    gwaslib.LOG.info("association table: %s", produced[0].name)
    return produced[0]


def load_glm(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", quoting=csv.QUOTE_NONE,
                     dtype={"#CHROM": str, "ID": str})
    df = df.rename(columns={"#CHROM": "CHROM"})
    if "TEST" in df.columns:                     # keep only the additive term
        df = df[df["TEST"] == "ADD"]
    df["CHROM"] = pd.to_numeric(df["CHROM"], errors="coerce")
    df = df.dropna(subset=["P", "CHROM"]).copy()
    df["CHROM"] = df["CHROM"].astype(int)
    df = df.sort_values(["CHROM", "POS"]).reset_index(drop=True)

    n_err = int((df.get("ERRCODE", pd.Series(["."] * len(df))) != ".").sum())
    gwaslib.LOG.info("%d variants tested, %d with a non-'.' ERRCODE", len(df), n_err)
    if "FIRTH?" in df.columns:
        gwaslib.LOG.info("Firth fallback used for %d variants",
                         int((df["FIRTH?"] == "Y").sum()))
    return df


def lambda_gc(p: pd.Series) -> float:
    """Genomic inflation factor: observed median chi2 over its null expectation.

    1.0 means no inflation.  Much above ~1.05 at this sample size points to
    residual structure or relatedness rather than polygenicity.
    """
    chi2 = stats.chi2.isf(p.clip(lower=1e-300), df=1)
    return float(np.median(chi2) / stats.chi2.ppf(0.5, df=1))


def top_hits(df: pd.DataFrame, threshold: float = SUGGESTIVE) -> pd.DataFrame:
    cols = [c for c in ["CHROM", "POS", "ID", "A1", "A1_FREQ", "OBS_CT", "OR",
                        "L95", "U95", "Z_STAT", "P", "FIRTH?"] if c in df.columns]
    return df.loc[df["P"] < threshold, cols].sort_values("P").reset_index(drop=True)


# --------------------------------------------------------------------------
# Independent loci
# --------------------------------------------------------------------------
def run_clump(glm_path: Path, bfile: Path = QC_STEM, out: Path | None = None,
              p1: float = SUGGESTIVE, p2: float = 1e-3,
              r2: float = 0.1, kb: int = 500) -> Path | None:
    """LD-clump the association results into independent signals.

    A single causal variant drags its whole LD block over the threshold, so the
    raw hit count badly overstates the number of findings.  Clumping keeps the
    lead SNP per block.  PLINK 1.9 is used because --clump was dropped from
    early plink2 builds.
    """
    out = out or gwaslib.ROOT_RESULTS / "penncath_clump"
    gwaslib.plink([
        "--bfile", bfile,
        "--clump", glm_path,
        "--clump-snp-field", "ID",
        "--clump-field", "P",
        "--clump-p1", p1, "--clump-p2", p2,
        "--clump-r2", r2, "--clump-kb", kb,
    ], out=out)
    clumped = Path(f"{out}.clumped")
    if not clumped.exists():
        gwaslib.LOG.info("no clumps formed at p1=%g (no variant passed)", p1)
        return None
    return clumped


def load_clumps(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", quoting=csv.QUOTE_NONE)
    return df.dropna(subset=["P"]).reset_index(drop=True)


# --------------------------------------------------------------------------
# Gene annotation (GRCh37)
# --------------------------------------------------------------------------
def _ensembl(path: str, retries: int = 3) -> object:
    url = f"{ENSEMBL_GRCH37}{path}"
    for attempt in range(retries):
        try:
            request = urllib.request.Request(
                url, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(request, timeout=45) as response:
                return json.loads(response.read().decode())
        except Exception as exc:                               # noqa: BLE001
            if attempt == retries - 1:
                gwaslib.LOG.warning("Ensembl query failed (%s): %s", path, exc)
                return None
            time.sleep(2 ** attempt)
    return None


def power_params(df: pd.DataFrame) -> tuple[int, float]:
    """(n, case fraction) read off the tested fileset, not hard-coded."""
    fam = pd.read_csv(f"{QC_STEM}.fam", sep=r"\s+", header=None, dtype=str,
                      quoting=csv.QUOTE_NONE,
                      names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"])
    n = len(fam)
    n_case = int((fam["PHENO"] == "2").sum())
    return n, n_case / n


def detectable_or(maf: float, alpha: float = GENOME_WIDE, power: float = 0.80,
                  n: int = 1382, case_frac: float = 923 / 1382) -> float:
    """Per-allele OR detectable at `power`, from the Wald variance of log-OR.

    Var(beta) ~ 1 / (2 N f(1-f) phi(1-phi)) for a per-allele log-odds, with f
    the minor allele frequency and phi the case fraction.  Asymptotic, so it is
    optimistic for very rare variants.
    """
    var_beta = 1.0 / (2 * n * maf * (1 - maf) * case_frac * (1 - case_frac))
    beta = (stats.norm.isf(alpha / 2) + stats.norm.ppf(power)) * np.sqrt(var_beta)
    return float(np.exp(beta))


def power_for(or_value: float, maf: float, alpha: float = GENOME_WIDE,
              n: int = 1382, case_frac: float = 923 / 1382) -> float:
    """Power to detect `or_value` at `alpha` -- the inverse of detectable_or."""
    var_beta = 1.0 / (2 * n * maf * (1 - maf) * case_frac * (1 - case_frac))
    z = abs(np.log(or_value)) / np.sqrt(var_beta)
    return float(stats.norm.sf(stats.norm.isf(alpha / 2) - z))


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------
# Typography is set inside each figure function rather than inherited from the
# caller's rcParams.  A figure that only lays out correctly when the notebook
# happens to have applied a style first is a figure that breaks silently: the
# tick labels collide at larger default sizes.  Three sizes, mapped to role.
FIG_RC = {
    "figure.dpi": 150, "savefig.dpi": 300,
    "font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9,
    "legend.fontsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.direction": "out", "ytick.direction": "out",
    "axes.titlelocation": "left", "axes.titlepad": 6,
}


def make_manhattan_qq(df: pd.DataFrame, clumps: pd.DataFrame,
                      out: Path | None = None,
                      labels: dict[str, str] | None = None) -> Path:
    """Manhattan over a QQ panel, lead SNPs circled."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = out or gwaslib.ROOT_FIGURES / "penncath_manhattan_qq.png"
    ctx = matplotlib.rc_context(FIG_RC)
    ctx.__enter__()
    df = df.copy()
    offset = df.groupby("CHROM")["POS"].max().cumsum().shift(fill_value=0)
    df["gpos"] = df["POS"] + df["CHROM"].map(offset)
    centres = df.groupby("CHROM")["gpos"].median()
    lam = lambda_gc(df["P"])

    fig = plt.figure(figsize=(7.2, 5.8))
    # hspace must clear the staggered chromosome row plus the x-axis label,
    # otherwise "chromosome" lands on the QQ panel title
    gs = fig.add_gridspec(2, 1, height_ratios=[1.35, 1], hspace=0.58)
    ax = fig.add_subplot(gs[0])
    grey, dark, focal = "#BFC4CC", "#8A9099", "#C44E52"

    for chrom, sub in df.groupby("CHROM"):
        ax.scatter(sub["gpos"], -np.log10(sub["P"]), s=1.6, linewidths=0,
                   color=grey if chrom % 2 else dark, rasterized=True)
    lead = df[df["ID"].isin(clumps["SNP"])]
    ax.scatter(lead["gpos"], -np.log10(lead["P"]), s=26, facecolor="none",
               edgecolor=focal, linewidths=1.2, zorder=5)

    ax.axhline(-np.log10(GENOME_WIDE), color="#333333", lw=0.9)
    ax.axhline(-np.log10(SUGGESTIVE), color="#333333", lw=0.8, ls=":")
    ax.text(df["gpos"].max(), -np.log10(GENOME_WIDE) + 0.12,
            "genome-wide 5\u00d710$^{-8}$", ha="right", va="bottom", fontsize=7)
    ax.text(df["gpos"].max(), -np.log10(SUGGESTIVE) + 0.12,
            "suggestive 10$^{-5}$", ha="right", va="bottom", fontsize=7)

    for _, row in lead.iterrows():
        if labels and row["ID"] in labels:
            ax.annotate(f"{labels[row['ID']]}\n{row['ID']}",
                        xy=(row["gpos"], -np.log10(row["P"])),
                        xytext=(row["gpos"], -np.log10(row["P"]) + 1.35),
                        ha="center", va="bottom", fontsize=7,
                        arrowprops=dict(arrowstyle="-", lw=0.7, color=focal,
                                        shrinkA=0, shrinkB=3))

    ax.set_xticks(centres.values)
    ticks = ax.set_xticklabels([str(c) for c in centres.index])
    # chromosomes 18+ are narrow enough that labels collide; drop alternate
    # ones to a second row rather than shrinking the whole size ladder
    for lab, chrom in zip(ticks, centres.index):
        if chrom >= 18 and chrom % 2 == 1:
            lab.set_y(lab.get_position()[1] - 0.055)
    ax.set_xlim(df["gpos"].min() - 2e7, df["gpos"].max() + 2e7)
    ax.set_ylim(0, max(9.6, -np.log10(df["P"].min()) + 2.5))
    ax.set_xlabel("chromosome", labelpad=14)
    ax.set_ylabel("$-\\log_{10}$ $P$")
    n_gw = int((df["P"] < GENOME_WIDE).sum())
    n, case_frac = power_params(df)
    ax.set_title(
        f"{'No locus reaches' if n_gw == 0 else f'{n_gw} variants reach'} "
        f"genome-wide significance in {round(n * case_frac)} cases / "
        f"{round(n * (1 - case_frac))} controls", loc="left")

    axq = fig.add_subplot(gs[1])
    p = np.sort(df["P"].values)
    expected = -np.log10((np.arange(1, len(p) + 1) - 0.5) / len(p))
    keep = np.r_[np.arange(0, min(20000, len(p))), np.arange(20000, len(p), 37)]
    axq.plot([0, expected.max()], [0, expected.max()], color="#333333", lw=0.9)
    axq.scatter(expected[keep], -np.log10(p)[keep], s=3, linewidths=0,
                color=dark, rasterized=True)
    axq.set_xlabel("expected $-\\log_{10}$ $P$")
    axq.set_ylabel("observed $-\\log_{10}$ $P$")
    axq.set_title(f"Test statistics are not inflated "
                  f"($\\lambda_{{GC}}$ = {lam:.3f})", loc="left")
    axq.margins(0.04)

    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    ctx.__exit__(None, None, None)
    gwaslib.LOG.info("wrote %s", out)
    return out


def make_locus_plots(df: pd.DataFrame, specs: list[tuple[str, str, str]],
                     out: Path | None = None) -> Path:
    """Regional plots coloured by LD to each lead.

    `specs` is [(lead_rsid, ld_file_stem, panel_title), ...]; the .ld files come
    from `plink --r2 --ld-snp <lead>` (see the notebook).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = out or gwaslib.ROOT_FIGURES / "penncath_loci.png"
    ctx = matplotlib.rc_context(FIG_RC)
    ctx.__enter__()
    bins = [(0.0, 0.2, "#4C72B0", "<0.2"), (0.2, 0.4, "#64B5CD", "0.2-0.4"),
            (0.4, 0.6, "#55A868", "0.4-0.6"), (0.6, 0.8, "#DD8452", "0.6-0.8"),
            (0.8, 1.01, "#C44E52", ">0.8")]

    fig, axes = plt.subplots(1, len(specs), figsize=(3.7 * len(specs), 3.2),
                             squeeze=False)
    for ax, (lead_id, stem, title) in zip(axes[0], specs):
        ld = pd.read_csv(gwaslib.ROOT_INTERIM / f"{stem}.ld", sep=r"\s+")
        sub = df.merge(ld[["SNP_B", "R2"]], left_on="ID", right_on="SNP_B")
        for lo, hi, colour, label in bins:
            m = (sub["R2"] >= lo) & (sub["R2"] < hi) & (sub["ID"] != lead_id)
            ax.scatter(sub.loc[m, "POS"] / 1e6, -np.log10(sub.loc[m, "P"]),
                       s=22, color=colour, linewidths=0, label=label, alpha=0.9)
        led = sub[sub["ID"] == lead_id]
        ax.scatter(led["POS"] / 1e6, -np.log10(led["P"]), s=90, marker="D",
                   facecolor="#8E44AD", edgecolor="white", linewidths=0.8, zorder=6)
        ax.annotate(lead_id, xy=(led["POS"].iloc[0] / 1e6,
                                 -np.log10(led["P"].iloc[0])),
                    xytext=(6, 4), textcoords="offset points", fontsize=7)
        ax.axhline(-np.log10(GENOME_WIDE), color="#333333", lw=0.8)
        ax.axhline(-np.log10(SUGGESTIVE), color="#333333", lw=0.8, ls=":")
        chrom = int(led["CHROM"].iloc[0])
        ax.set_xlabel(f"chr{chrom} position (Mb, GRCh37)")
        ax.set_ylabel("$-\\log_{10}$ $P$")
        ax.set_title(title, loc="left")
        ax.margins(0.04)

    handles, labels_ = axes[0][0].get_legend_handles_labels()
    axes[0][-1].legend(handles, labels_, title="$r^2$ to lead", frameon=False,
                       fontsize=6.5, title_fontsize=6.5, loc="upper right",
                       handletextpad=0.3, borderaxespad=0.2, labelspacing=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    ctx.__exit__(None, None, None)
    gwaslib.LOG.info("wrote %s", out)
    return out


def annotate_loci(hits: pd.DataFrame, window_kb: int = 200) -> pd.DataFrame:
    """For each lead variant, list overlapping and nearby protein-coding genes.

    Reports the gene containing the variant when there is one, plus the nearest
    genes within +/- window_kb.  A GWAS hit is not evidence that the nearest
    gene is the causal one -- most lead variants are non-coding and act through
    regulatory elements that can skip over intervening genes -- so this is a
    locus label, not a causal assignment.
    """
    rows = []
    for _, hit in hits.iterrows():
        chrom, pos = int(hit["CHROM"]), int(hit["POS"])
        start, end = max(1, pos - window_kb * 1000), pos + window_kb * 1000
        genes = _ensembl(
            f"/overlap/region/human/{chrom}:{start}-{end}"
            f"?feature=gene;content-type=application/json") or []
        coding = [g for g in genes if g.get("biotype") == "protein_coding"]
        other = [g for g in genes if g.get("biotype") != "protein_coding"]

        def describe(gene: dict) -> str:
            name = gene.get("external_name") or gene.get("gene_id")
            if gene["start"] <= pos <= gene["end"]:
                return f"{name} (within)"
            distance = min(abs(pos - gene["start"]), abs(pos - gene["end"]))
            return f"{name} ({distance // 1000} kb)"

        within = [g for g in genes if g["start"] <= pos <= g["end"]]
        rows.append({
            "ID": hit["ID"], "CHROM": chrom, "POS": pos, "P": hit["P"],
            "OR": hit.get("OR"),
            "gene_at_variant": ", ".join(
                (g.get("external_name") or g["gene_id"]) for g in within) or "-",
            "protein_coding_within_window": ", ".join(
                describe(g) for g in sorted(coding, key=lambda g: g["start"])) or "-",
            "other_biotypes": ", ".join(
                describe(g) for g in sorted(other, key=lambda g: g["start"])) or "-",
        })
        time.sleep(0.2)                       # stay polite to the REST endpoint
    return pd.DataFrame(rows)
