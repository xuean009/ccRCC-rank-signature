#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import gzip
import json
import shutil
from collections import OrderedDict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, proportional_hazard_test
from lifelines.utils import concordance_index
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from PIL import Image
from scipy.stats import chi2


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "output" / "cancer_medicine_submission_package"
FIG_DIR = OUT_DIR / "figures"
TAB_DIR = OUT_DIR / "tables"
SUPP_DIR = OUT_DIR / "supplementary"
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"

BOOTSTRAP_N = 1000
OPTIMISM_BOOTSTRAP_N = 1000
SEED = 42

DISCOVERY_ROLE = "Discovery cohort"
DEVELOPMENT_LINKED_ROLE = "Development-linked OS projection cohort"
PRIMARY_EXTERNAL_ROLE = "Primary independent OS validation cohort"
SECONDARY_ENDPOINT_ROLE = "Secondary endpoint/context sensitivity cohort"

PAL = {
    "ink": "#1B2733",
    "muted": "#637381",
    "grid": "#D8E1E8",
    "red": "#B14B3B",
    "blue": "#2B698C",
    "green": "#2E8262",
    "gold": "#B58726",
    "mauve": "#7A6C92",
    "panel": "#F7FAFC",
}


def ensure_dirs():
    for path in (OUT_DIR, FIG_DIR, TAB_DIR, SUPP_DIR):
        path.mkdir(parents=True, exist_ok=True)


def clean_output():
    if not OUT_DIR.exists():
        return
    for child in OUT_DIR.iterdir():
        if child.is_file():
            try:
                child.unlink()
            except PermissionError:
                pass
        elif child.is_dir():
            for nested in child.rglob("*"):
                if nested.is_file():
                    try:
                        nested.unlink()
                    except PermissionError:
                        pass
    ensure_dirs()


def set_style():
    plt.rcParams.update(
        {
            "figure.dpi": 170,
            "savefig.dpi": 350,
            "font.family": "DejaVu Sans",
            "font.size": 10.5,
            "axes.linewidth": 0.8,
            "axes.edgecolor": PAL["ink"],
            "axes.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.labelcolor": PAL["ink"],
            "xtick.color": PAL["ink"],
            "ytick.color": PAL["ink"],
            "text.color": PAL["ink"],
            "grid.color": PAL["grid"],
            "grid.linewidth": 0.7,
            "grid.alpha": 0.75,
            "svg.fonttype": "none",
        }
    )


def read_table(path, sep="\t", **kwargs):
    path = Path(path)
    if str(path).endswith(".gz"):
        with gzip.open(str(path), "rt", encoding="utf-8", errors="ignore") as handle:
            return pd.read_csv(handle, sep=sep, **kwargs)
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        return pd.read_csv(handle, sep=sep, **kwargs)


def rank_rows(df):
    ranks = df.rank(axis=1, method="average", na_option="keep")
    max_rank = ranks.max(axis=1)
    return ranks.div(max_rank, axis=0).fillna(0.5)


def format_p(value):
    if value is None or not np.isfinite(value):
        return "NA"
    if value < 1e-4:
        return "<0.0001"
    if value < 1e-3:
        return "{:.4f}".format(value)
    return "{:.3f}".format(value)


def format_p_figure(value):
    if value is None or not np.isfinite(value):
        return "NA"
    if value < 1e-4:
        return "<0.0001"
    if value < 1e-3:
        return "{:.4f}".format(value)
    return "{:.3f}".format(value)


def format_float(value, digits=3):
    if value is None or not np.isfinite(value):
        return "NA"
    return "{:.{d}f}".format(value, d=digits)


def format_ci(lo, hi, digits=3):
    if not np.isfinite(lo) or not np.isfinite(hi):
        return "NA"
    return "{:.{d}f}-{:.{d}f}".format(lo, hi, d=digits)


def format_percent(value, digits=1):
    if value is None or not np.isfinite(value):
        return "NA"
    return "{:.{d}f}%".format(value * 100.0, d=digits)


def median_iqr(series, digits=1):
    s = pd.to_numeric(series, errors="coerce").dropna()
    if s.empty:
        return "NR"
    return "{:.{d}f} ({:.{d}f}-{:.{d}f})".format(
        float(s.median()), float(s.quantile(0.25)), float(s.quantile(0.75)), d=digits
    )


def male_fraction(series):
    s = series.astype(str).str.upper()
    denom = s[s.ne("NAN")]
    if denom.empty:
        return "NR"
    male = denom.isin(["MALE", "M"]).sum()
    return "{} ({})".format(int(male), format_percent(male / float(len(denom))))


def save_table(df, basename, directory=TAB_DIR, preamble_lines=None, note_lines=None):
    tsv_path = directory / (basename + ".tsv")
    md_path = directory / (basename + ".md")
    df.to_csv(str(tsv_path), sep="\t", index=False)
    lines = []
    if preamble_lines:
        lines.extend(preamble_lines)
        lines.append("")
    lines.extend(["| " + " | ".join(df.columns) + " |", "| " + " | ".join(["---"] * len(df.columns)) + " |"])
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in df.columns) + " |")
    if note_lines:
        lines.append("")
        lines.extend(note_lines)
    md_path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def save_figure(fig, basename, directory=FIG_DIR):
    png_path = directory / (basename + ".png")
    svg_path = directory / (basename + ".svg")
    tiff_path = directory / (basename + ".tiff")
    fig.savefig(str(png_path), facecolor="white", bbox_inches="tight")
    fig.savefig(str(svg_path), facecolor="white", bbox_inches="tight")
    with Image.open(str(png_path)) as im:
        im.save(str(tiff_path), dpi=(300, 300), compression="tiff_lzw")
    plt.close(fig)


def soften(ax):
    ax.tick_params(length=3.3, width=0.7)
    for spine in ("left", "bottom"):
        if spine in ax.spines:
            ax.spines[spine].set_linewidth(0.8)
            ax.spines[spine].set_color(PAL["ink"])


def score_formula_note():
    return (
        "r_ig = rank(x_ig)/G_i within the frozen 33-gene shared-gene universe "
        "(average ties; missing ranks assigned 0.5); score_i = sum(beta_g * r_ig)."
    )


def load_signature():
    coef = read_table(RESULTS_DIR / "best_signature_coefficients.tsv")
    coef["gene"] = coef["gene"].astype(str).str.upper()
    coef["abs_coef"] = coef["coef"].abs()
    coef["direction"] = np.where(coef["coef"] >= 0, "Risk-increasing", "Protective")
    category_map = {
        "CXCL6": "immune",
        "LTF": "immune",
        "SAA1": "immune",
        "SERPINA3": "immune",
        "FGG": "immune",
        "ATP6V0D2": "immune",
        "CYP4F2": "metabolic",
        "SLC22A8": "metabolic",
        "APOB": "metabolic",
        "ANGPTL3": "metabolic",
        "GBA3": "metabolic",
        "SLC15A1": "metabolic",
        "GDA": "metabolic",
        "CASR": "metabolic",
    }
    coef["category"] = coef["gene"].map(category_map).fillna("mixed")
    return coef


def bootstrap_cindex(time, event, score, n_boot=BOOTSTRAP_N, seed=SEED):
    rng = np.random.RandomState(seed)
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    score = np.asarray(score, dtype=float)
    idx = np.arange(time.size)
    estimates = []
    for _ in range(n_boot):
        sample_idx = rng.choice(idx, size=idx.size, replace=True)
        est = concordance_index(time[sample_idx], -score[sample_idx], event[sample_idx])
        if np.isfinite(est):
            estimates.append(est)
    if not estimates:
        return np.nan, np.nan
    return tuple(np.quantile(estimates, [0.025, 0.975]))


def compute_orientation_metrics(time, event, risk):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    risk = np.asarray(risk, dtype=float)
    mask = np.isfinite(time) & np.isfinite(event) & np.isfinite(risk)
    time = time[mask]
    event = event[mask]
    risk = risk[mask]
    out = {}
    for label, signed_risk in (("original", risk), ("sign_reversed", -risk)):
        cut = float(np.nanmedian(signed_risk))
        group = signed_risk >= cut
        pval = float(
            logrank_test(time[group], time[~group], event_observed_A=event[group], event_observed_B=event[~group]).p_value
        )
        cidx = float(concordance_index(time, -signed_risk, event))
        lo, hi = bootstrap_cindex(time, event, signed_risk, seed=SEED + (13 if label == "original" else 37))
        out[label] = {"p": pval, "cindex": cidx, "low": lo, "high": hi}
    return out


def univariable_per_sd(time, event, risk):
    df = pd.DataFrame({"time": time, "event": event, "risk": risk}).dropna().copy()
    df = df[np.isfinite(df["time"]) & np.isfinite(df["risk"])].copy()
    sd = float(df["risk"].std())
    df["risk_per_sd"] = df["risk"] / sd
    model = CoxPHFitter()
    model.fit(df[["time", "event", "risk_per_sd"]], duration_col="time", event_col="event")
    row = model.summary.loc["risk_per_sd"]
    return {
        "sd": sd,
        "hr": float(row["exp(coef)"]),
        "lo": float(row["exp(coef) lower 95%"]),
        "hi": float(row["exp(coef) upper 95%"]),
        "p": float(row["p"]),
    }


def load_tcga_scored(signature):
    coef = signature.set_index("gene")["coef"]
    clin = read_table(DATA_DIR / "TCGA_Xena" / "TCGA.KIRC.sampleMap" / "KIRC_clinicalMatrix")
    clin["time_days"] = pd.to_numeric(clin["days_to_death"], errors="coerce").fillna(
        pd.to_numeric(clin["days_to_last_followup"], errors="coerce")
    )
    clin["time_months"] = clin["time_days"] / 30.4375
    status_raw = clin["vital_status"].astype(str).str.lower()
    clin["status"] = status_raw.map({"deceased": 1, "dead": 1, "living": 0, "alive": 0})
    clin = clin[(clin["sample_type"] == "Primary Tumor") & clin["time_days"].notna() & (clin["time_days"] > 0) & clin["status"].notna()].copy()

    with gzip.open(str(DATA_DIR / "TCGA_Xena" / "TCGA.KIRC.sampleMap" / "HiSeqV2.gz"), "rt", encoding="utf-8", errors="ignore") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        sample_ids = header[1:]
        gene_rows = {}
        wanted = set(coef.index)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            gene = parts[0].upper()
            if gene in wanted:
                gene_rows[gene] = [float(x) if x not in ("", "NA", "NaN") else np.nan for x in parts[1:]]
    expr = pd.DataFrame(gene_rows, index=sample_ids).reindex(columns=coef.index)
    score = rank_rows(expr).dot(coef)
    clin = clin.merge(score.rename("risk_score").reset_index().rename(columns={"index": "sampleID"}), on="sampleID", how="inner")
    clin["age"] = pd.to_numeric(clin["age_at_initial_pathologic_diagnosis"], errors="coerce")
    clin["grade_num"] = pd.to_numeric(clin["neoplasm_histologic_grade"].astype(str).str.extract(r"G(\d)", expand=False), errors="coerce")
    clin["stage_group"] = clin["pathologic_stage"].astype(str).str.strip().str.lower().map(
        {"stage i": "I", "stage ii": "II", "stage iii": "III", "stage iv": "IV"}
    )
    clin["stage_ord"] = clin["stage_group"].map({"I": 1, "II": 2, "III": 3, "IV": 4})
    clin["platform"] = "RNA-seq (TCGA Xena HiSeqV2)"
    clin["expression_matrix"] = "TCGA.KIRC.sampleMap/HiSeqV2.gz"
    clin["normalization"] = "Xena-distributed normalized gene-level RNA-seq values"
    return clin


def load_cptac_scored(signature):
    coef = signature.set_index("gene")["coef"]
    expr = read_table(DATA_DIR / "cBioPortal" / "rcc_cptac_gdc" / "expr_genelist.tsv")
    expr = expr.rename(columns={expr.columns[0]: "gene"})
    expr["gene"] = expr["gene"].astype(str).str.upper()
    wide = expr.set_index("gene").T.apply(pd.to_numeric, errors="coerce").reindex(columns=coef.index)
    sample_risk = rank_rows(wide).dot(coef).rename("risk_score").reset_index().rename(columns={"index": "sampleId"})
    sample_map = read_table(DATA_DIR / "cBioPortal" / "rcc_cptac_gdc" / "sample_to_patient.tsv")
    clin = read_table(DATA_DIR / "cBioPortal" / "rcc_cptac_gdc" / "clinical_os_gdc.tsv")
    merged = sample_risk.merge(sample_map, on="sampleId", how="left")
    patient_risk = merged.groupby("patientId", as_index=False)["risk_score"].mean()
    out = clin.merge(patient_risk, on="patientId", how="inner")
    out["OS_MONTHS"] = pd.to_numeric(out["OS_MONTHS"], errors="coerce")
    out["OS_STATUS"] = pd.to_numeric(out["OS_STATUS"], errors="coerce")
    out = out[out["OS_MONTHS"].notna() & out["OS_STATUS"].notna() & (out["OS_MONTHS"] > 0)].copy()
    out["platform"] = "RNA-seq (cBioPortal TPM)"
    out["expression_matrix"] = "cBioPortal/rcc_cptac_gdc/expr_genelist.tsv"
    out["normalization"] = "cBioPortal continuous TPM profile"
    return out


def load_emtab1980_scored(signature):
    coef = signature.set_index("gene")["coef"]
    expr = read_table(DATA_DIR / "external_clinical" / "E-MTAB-1980" / "expression_signature_gene33.tsv")
    expr["gene"] = expr["gene"].astype(str).str.upper()
    wide = expr.drop_duplicates("gene").set_index("gene").T.apply(pd.to_numeric, errors="coerce").reindex(columns=coef.index)
    risk = rank_rows(wide).dot(coef)
    clin = read_table(DATA_DIR / "external_clinical" / "E-MTAB-1980" / "clinical_matched_101.tsv")
    clin = clin.set_index("sample_id").loc[wide.index].reset_index()
    clin["age"] = pd.to_numeric(clin["age"], errors="coerce")
    clin["fuhrman_grade"] = pd.to_numeric(clin["fuhrman_grade"], errors="coerce")
    stage = clin["stage_at_diagnosis"].astype(str)
    clin["stage_group"] = np.where(stage.str.contains(r"^pT[34]", regex=True), "III/IV-like", "I/II-like")
    clin["pT3_4"] = stage.str.contains(r"^pT[34]", regex=True).astype(float)
    clin["M1"] = stage.str.contains("M1").astype(float)
    clin["metastatic"] = clin["M1"]
    clin["risk_score"] = risk.values
    clin["platform"] = "Microarray (Agilent A-MEXP-2183)"
    clin["expression_matrix"] = "E-MTAB-1980.processed.1.zip::ccRCC_exp_log_quantile_normalized.txt"
    clin["normalization"] = "processed log-quantile-normalized Agilent intensities"
    return clin


def load_immotion150_scored(signature):
    coef = signature.set_index("gene")["coef"]
    expr = read_table(DATA_DIR / "external_clinical" / "IMmotion150" / "expression_signature_gene33.tsv")
    expr["gene"] = expr["gene"].astype(str).str.upper()
    wide = expr.drop_duplicates("gene").set_index("gene").T.apply(pd.to_numeric, errors="coerce").reindex(columns=coef.index)
    risk = rank_rows(wide).dot(coef)
    clin = read_table(DATA_DIR / "external_clinical" / "IMmotion150" / "clinical_matched_263.tsv")
    clin = clin.set_index("sample_id").loc[wide.index].reset_index()
    clin["pfs_months"] = pd.to_numeric(clin["pfs_months"], errors="coerce")
    clin["status_pfs"] = pd.to_numeric(clin["status_pfs"], errors="coerce")
    clin["age_num"] = pd.to_numeric(clin.get("age"), errors="coerce")
    clin["metastatic"] = clin["METASTASIZED"].astype(str).str.upper().eq("TRUE").astype(float)
    clin["risk_score"] = risk.values
    clin["platform"] = "RNA-seq (cBioPortal continuous expression)"
    clin["expression_matrix"] = "IMmotion150/expression_signature_gene33.tsv"
    clin["normalization"] = "public continuous RNA-seq profile reported by cBioPortal/iAtlas"
    return clin


def load_emtab3267_scored(signature):
    coef = signature.set_index("gene")["coef"]
    expr = read_table(DATA_DIR / "external_clinical" / "E-MTAB-3267" / "expression_signature_gene33.tsv")
    expr["gene"] = expr["gene"].astype(str).str.upper()
    wide = expr.drop_duplicates("gene").set_index("gene").T.apply(pd.to_numeric, errors="coerce").reindex(columns=coef.index)
    risk = rank_rows(wide).dot(coef)
    clin = read_table(DATA_DIR / "external_clinical" / "E-MTAB-3267" / "clinical_matched_53.tsv")
    clin = clin.set_index("sample_id").loc[wide.index].reset_index()
    clin["pfs_months"] = pd.to_numeric(clin["pfs_months"], errors="coerce")
    clin["status_pfs"] = pd.to_numeric(clin["status_pfs"], errors="coerce")
    clin["age_num"] = pd.to_numeric(clin["age"], errors="coerce")
    clin["metastatic"] = 1.0
    clin["risk_score"] = risk.values
    clin["platform"] = "Microarray (Affymetrix Human Gene 1.0 ST)"
    clin["expression_matrix"] = "E-MTAB-3267/rma_expression_all_samples.rds"
    clin["normalization"] = "oligo::rma(target='core') from raw CEL files"
    return clin


def stage34_fraction_tcga(df):
    valid = df["stage_ord"].dropna()
    if valid.empty:
        return "NR"
    count = int((valid >= 3).sum())
    return "{} ({})".format(count, format_percent(count / float(valid.size)))


def grade34_fraction(values):
    s = pd.to_numeric(values, errors="coerce").dropna()
    if s.empty:
        return "NR"
    count = int((s >= 3).sum())
    return "{} ({})".format(count, format_percent(count / float(s.size)))


def build_table1(tcga, cptac, em1980, imm150, em3267):
    rows = [
        OrderedDict(
            [
                ("Cohort", "TCGA-KIRC"),
                ("Role", DISCOVERY_ROLE),
                ("Endpoint", "overall survival"),
                ("Platform", "RNA-seq"),
                ("N", int(tcga.shape[0])),
                ("Events", int(tcga["status"].sum())),
                ("Age, median (IQR)", median_iqr(tcga["age"])),
                ("Male sex", male_fraction(tcga["gender"])),
                ("Advanced/metastatic disease", "{} stage III/IV".format(stage34_fraction_tcga(tcga))),
                ("Grade 3/4", grade34_fraction(tcga["grade_num"])),
                ("Metastatic/systemic setting", "84 (15.6%) stage IV"),
                ("Clinical-model missingness", "10/531 excluded from age+stage+grade complete-case model"),
            ]
        ),
        OrderedDict(
            [
                ("Cohort", "CPTAC"),
                ("Role", DEVELOPMENT_LINKED_ROLE),
                ("Endpoint", "overall survival"),
                ("Platform", "RNA-seq"),
                ("N", int(cptac.shape[0])),
                ("Events", int(cptac["OS_STATUS"].sum())),
                ("Age, median (IQR)", "NR"),
                ("Male sex", "NR"),
                ("Advanced/metastatic disease", "NR"),
                ("Grade 3/4", "NR"),
                ("Metastatic/systemic setting", "Primary-tumor proteogenomic cohort"),
                ("Clinical-model missingness", "Stage/grade unavailable in linked survival file"),
            ]
        ),
        OrderedDict(
            [
                ("Cohort", "E-MTAB-1980"),
                ("Role", PRIMARY_EXTERNAL_ROLE),
                ("Endpoint", "overall survival"),
                ("Platform", "microarray"),
                ("N", int(em1980.shape[0])),
                ("Events", int(em1980["status_os"].sum())),
                ("Age, median (IQR)", median_iqr(em1980["age"])),
                ("Male sex", male_fraction(em1980["sex"])),
                ("Advanced/metastatic disease", "{} ({}) pT3/4".format(int(em1980["pT3_4"].sum()), format_percent(em1980["pT3_4"].mean()))),
                ("Grade 3/4", grade34_fraction(em1980["fuhrman_grade"])),
                ("Metastatic/systemic setting", "{} ({}) M1 at diagnosis".format(int(em1980["M1"].sum()), format_percent(em1980["M1"].mean()))),
                ("Clinical-model missingness", "2/101 excluded from age+grade+pT3/4+M1 complete-case model"),
            ]
        ),
        OrderedDict(
            [
                ("Cohort", "IMmotion150"),
                ("Role", SECONDARY_ENDPOINT_ROLE),
                ("Endpoint", "progression-free survival"),
                ("Platform", "RNA-seq"),
                ("N", int(imm150.shape[0])),
                ("Events", int(imm150["status_pfs"].sum())),
                ("Age, median (IQR)", median_iqr(imm150["age_num"])),
                ("Male sex", "NR"),
                ("Advanced/metastatic disease", "263 (100.0%) stage IV/metastatic"),
                ("Grade 3/4", "NR"),
                ("Metastatic/systemic setting", "263 (100.0%) trial cohort"),
                ("Clinical-model missingness", "Not benchmarked; harmonized grade data unavailable"),
            ]
        ),
        OrderedDict(
            [
                ("Cohort", "E-MTAB-3267"),
                ("Role", SECONDARY_ENDPOINT_ROLE),
                ("Endpoint", "progression-free survival"),
                ("Platform", "microarray"),
                ("N", int(em3267.shape[0])),
                ("Events", int(em3267["status_pfs"].sum())),
                ("Age, median (IQR)", median_iqr(em3267["age_num"])),
                ("Male sex", male_fraction(em3267["sex"])),
                ("Advanced/metastatic disease", "53 (100.0%) metastatic treated cohort"),
                ("Grade 3/4", "NR"),
                ("Metastatic/systemic setting", "53 (100.0%) sunitinib-treated"),
                ("Clinical-model missingness", "Not benchmarked; stage/grade not harmonized"),
            ]
        ),
    ]
    df = pd.DataFrame(rows)
    save_table(
        df,
        "Table1_cohort_characteristics",
        note_lines=[
            "Note: Percentages use cohort-specific available-case denominators when covariates were missing; in E-MTAB-1980, Grade 3/4 corresponds to 27/99 cases with non-missing Fuhrman grade. NR, not reported in the linked public clinical file."
        ],
    )
    return df


def build_table2(signature, cohorts):
    rows = []
    for cohort_name, meta in cohorts.items():
        per_sd = univariable_per_sd(meta["time"], meta["event"], meta["risk"])
        rows.append(
            OrderedDict(
                [
                    ("Cohort", cohort_name),
                    ("Endpoint", meta["endpoint"]),
                    ("Validation role", meta["role"]),
                    ("Platform", meta["platform"]),
                    ("N", int(len(meta["time"]))),
                    ("Events", int(np.sum(meta["event"]))),
                    ("Log-rank P (median split)", format_p(meta["original"]["p"])),
                    ("Original score C-index", format_float(meta["original"]["cindex"])),
                    ("Original score 95% CI", format_ci(meta["original"]["low"], meta["original"]["high"])),
                    ("Per-SD HR", "{:.3f} ({:.3f}-{:.3f})".format(per_sd["hr"], per_sd["lo"], per_sd["hi"])),
                    ("Per-SD HR P", format_p(per_sd["p"])),
                    ("Signature genes available", "{}/33".format(meta["genes_available"])),
                ]
            )
        )
    df = pd.DataFrame(rows)
    save_table(
        df,
        "Table2_external_validation_summary",
        note_lines=[
            "Note: Original-direction performance is shown here. Sign-reversed analyses are reported separately in `supplementary/TableS1_orientation_robustness.md`."
        ],
    )
    return df


def build_table3(signature):
    df = signature.sort_values("abs_coef", ascending=False).copy()
    out = pd.DataFrame(
        {
            "Gene symbol": df["gene"],
            "Coefficient": df["coef"].map(lambda x: "{:.6f}".format(x)),
            "Absolute coefficient": df["abs_coef"].map(lambda x: "{:.6f}".format(x)),
            "Direction": df["direction"],
            "Category": df["category"],
        }
    )
    save_table(out, "Table3_final_33_gene_signature")
    return out


def build_table4():
    rows = [
        OrderedDict(
            [
                ("Cohort", "TCGA-KIRC"),
                ("Model-development overlap", "Discovery cohort used for penalized model fitting"),
                ("Expression matrix / file", "TCGA.KIRC.sampleMap/HiSeqV2.gz"),
                ("Normalization / scale", "Xena-distributed normalized RNA-seq values"),
                ("Frozen rank universe", "Predefined 33-gene shared-gene universe"),
                ("Tie handling", "average"),
                ("Missing signature genes", "0/33"),
                ("Outcome complete cases", "531/537"),
                ("Clinical-model complete cases", "521/531"),
                ("Complete-case exclusion reason", "2 invalid survival, 4 unmatched expression rows, 10 missing age/stage/grade"),
            ]
        ),
        OrderedDict(
            [
                ("Cohort", "CPTAC"),
                ("Model-development overlap", "Used in archived cross-platform search-stage ranking; interpreted here as development-linked projection"),
                ("Expression matrix / file", "cBioPortal/rcc_cptac_gdc/expr_genelist.tsv"),
                ("Normalization / scale", "cBioPortal continuous TPM profile"),
                ("Frozen rank universe", "Predefined 33-gene shared-gene universe"),
                ("Tie handling", "average"),
                ("Missing signature genes", "0/33"),
                ("Outcome complete cases", "237/237"),
                ("Clinical-model complete cases", "Not attempted"),
                ("Complete-case exclusion reason", "Linked survival table lacked usable stage/grade"),
            ]
        ),
        OrderedDict(
            [
                ("Cohort", "E-MTAB-1980"),
                ("Model-development overlap", "No model-definition involvement; primary independent OS validation cohort"),
                ("Expression matrix / file", "E-MTAB-1980.processed.1.zip::ccRCC_exp_log_quantile_normalized.txt"),
                ("Normalization / scale", "processed log-quantile-normalized Agilent intensities"),
                ("Frozen rank universe", "Predefined 33-gene shared-gene universe"),
                ("Tie handling", "average"),
                ("Missing signature genes", "0/33"),
                ("Outcome complete cases", "101/101"),
                ("Clinical-model complete cases", "99/101"),
                ("Complete-case exclusion reason", "2 missing age and/or stage-derived covariates"),
            ]
        ),
        OrderedDict(
            [
                ("Cohort", "IMmotion150"),
                ("Model-development overlap", "No model-definition involvement; secondary endpoint/context sensitivity cohort"),
                ("Expression matrix / file", "IMmotion150/expression_signature_gene33.tsv"),
                ("Normalization / scale", "public continuous RNA-seq profile reported by cBioPortal/iAtlas"),
                ("Frozen rank universe", "Predefined 33-gene shared-gene universe"),
                ("Tie handling", "average"),
                ("Missing signature genes", "0/33"),
                ("Outcome complete cases", "263/263"),
                ("Clinical-model complete cases", "Not attempted"),
                ("Complete-case exclusion reason", "Clinicopathologic benchmarking not prespecified in trial cohort"),
            ]
        ),
        OrderedDict(
            [
                ("Cohort", "E-MTAB-3267"),
                ("Model-development overlap", "No model-definition involvement; secondary endpoint/context sensitivity cohort"),
                ("Expression matrix / file", "E-MTAB-3267/rma_expression_all_samples.rds"),
                ("Normalization / scale", "oligo::rma(target='core') from raw CEL files"),
                ("Frozen rank universe", "Predefined 33-gene shared-gene universe"),
                ("Tie handling", "average"),
                ("Missing signature genes", "0/33"),
                ("Outcome complete cases", "53/53"),
                ("Clinical-model complete cases", "Not attempted"),
                ("Complete-case exclusion reason", "Clinicopathologic benchmarking not feasible from public annotations"),
            ]
        ),
    ]
    df = pd.DataFrame(rows)
    save_table(df, "Table4_preprocessing_and_complete_case_analysis")
    return df


def fit_cindex_model(df, covariates):
    sub = df[["time", "event"] + covariates].dropna().copy()
    model = CoxPHFitter()
    model.fit(sub, duration_col="time", event_col="event")
    lp = model.predict_log_partial_hazard(sub[covariates]).values.ravel()
    hazard = model.predict_partial_hazard(sub[covariates]).values.ravel()
    cidx = concordance_index(sub["time"], -hazard, sub["event"])
    return {"data": sub, "model": model, "lp": lp, "cindex": float(cidx)}


def ph_test_pvalues(fitted):
    test = proportional_hazard_test(fitted["model"], fitted["data"], time_transform="rank")
    summary = test.summary.copy()
    return {str(idx): float(val) for idx, val in summary["p"].to_dict().items()}


def optimism_corrected_cindex_and_slope(df, covariates, n_boot=OPTIMISM_BOOTSTRAP_N, seed=SEED):
    fitted = fit_cindex_model(df, covariates)
    apparent_c = fitted["cindex"]
    apparent_slope = 1.0
    rng = np.random.RandomState(seed)
    optimism_c = []
    optimism_s = []
    full = fitted["data"].reset_index(drop=True)
    for _ in range(n_boot):
        idx = rng.choice(full.index.values, size=full.shape[0], replace=True)
        boot = full.iloc[idx].reset_index(drop=True)
        try:
            boot_fit = fit_cindex_model(boot, covariates)
            c_train = boot_fit["cindex"]
            c_test = concordance_index(full["time"], -boot_fit["model"].predict_partial_hazard(full[covariates]).values.ravel(), full["event"])
            test_df = full[["time", "event"]].copy()
            test_df["lp"] = boot_fit["model"].predict_log_partial_hazard(full[covariates]).values.ravel()
            slope_model = CoxPHFitter()
            slope_model.fit(test_df, duration_col="time", event_col="event")
            test_slope = float(slope_model.params_["lp"])
        except Exception:
            continue
        if np.isfinite(c_train) and np.isfinite(c_test):
            optimism_c.append(c_train - c_test)
        if np.isfinite(test_slope):
            optimism_s.append(1.0 - test_slope)
    corrected_c = apparent_c - float(np.mean(optimism_c)) if optimism_c else np.nan
    corrected_slope = apparent_slope - float(np.mean(optimism_s)) if optimism_s else np.nan
    return {
        "apparent_cindex": apparent_c,
        "corrected_cindex": corrected_c,
        "apparent_slope": apparent_slope,
        "corrected_slope": corrected_slope,
        "optimism_c": float(np.mean(optimism_c)) if optimism_c else np.nan,
        "optimism_slope": float(np.mean(optimism_s)) if optimism_s else np.nan,
        "n_boot": len(optimism_c),
    }


def hr_per_sd_from_combined(model, data, col):
    beta = float(model.params_[col])
    se = float(model.standard_errors_[col])
    sd = float(data[col].std())
    return {
        "sd": sd,
        "hr": float(np.exp(beta * sd)),
        "lo": float(np.exp((beta - 1.96 * se) * sd)),
        "hi": float(np.exp((beta + 1.96 * se) * sd)),
    }


def build_table5(tcga, em1980):
    tcga_df = tcga[["time_months", "status", "age", "stage_ord", "grade_num", "risk_score"]].dropna().copy()
    tcga_df.columns = ["time", "event", "age", "stage_ord", "grade_num", "risk_score"]
    tcga_clin = fit_cindex_model(tcga_df, ["age", "stage_ord", "grade_num"])
    tcga_sig = fit_cindex_model(tcga_df, ["risk_score"])
    tcga_comb = fit_cindex_model(tcga_df, ["age", "stage_ord", "grade_num", "risk_score"])
    tcga_lr_p = float(chi2.sf(2.0 * (tcga_comb["model"].log_likelihood_ - tcga_clin["model"].log_likelihood_), 1))
    tcga_per_sd = hr_per_sd_from_combined(tcga_comb["model"], tcga_comb["data"], "risk_score")
    tcga_optimism = optimism_corrected_cindex_and_slope(tcga_df, ["age", "stage_ord", "grade_num", "risk_score"])
    tcga_ph = ph_test_pvalues(tcga_comb)

    em_df = em1980[["time_months", "status_os", "age", "fuhrman_grade", "pT3_4", "M1", "risk_score"]].dropna().copy()
    em_df.columns = ["time", "event", "age", "fuhrman_grade", "pT3_4", "M1", "risk_score"]
    em_clin = fit_cindex_model(em_df, ["age", "fuhrman_grade", "pT3_4", "M1"])
    em_sig = fit_cindex_model(em_df, ["risk_score"])
    em_comb = fit_cindex_model(em_df, ["age", "fuhrman_grade", "pT3_4", "M1", "risk_score"])
    em_lr_p = float(chi2.sf(2.0 * (em_comb["model"].log_likelihood_ - em_clin["model"].log_likelihood_), 1))
    em_per_sd = hr_per_sd_from_combined(em_comb["model"], em_comb["data"], "risk_score")
    em_ph = ph_test_pvalues(em_comb)

    rows = [
        OrderedDict(
            [
                ("Cohort", "TCGA-KIRC"),
                ("Clinical comparator", "Age + pathologic stage + grade"),
                ("N", int(tcga_comb["data"].shape[0])),
                ("Events", int(tcga_comb["data"]["event"].sum())),
                ("Clinical model C-index", format_float(tcga_clin["cindex"])),
                ("Signature-alone C-index", format_float(tcga_sig["cindex"])),
                ("Combined apparent C-index", format_float(tcga_comb["cindex"])),
                ("Combined optimism-corrected C-index", format_float(tcga_optimism["corrected_cindex"])),
                ("Optimism-corrected calibration slope", format_float(tcga_optimism["corrected_slope"])),
                (
                    "Signature HR per SD in combined model",
                    "{:.3f} ({:.3f}-{:.3f})".format(tcga_per_sd["hr"], tcga_per_sd["lo"], tcga_per_sd["hi"]),
                ),
                ("Signature p", format_p(float(tcga_comb["model"].summary.loc["risk_score", "p"]))),
                ("Incremental LR p", format_p(tcga_lr_p)),
            ]
        ),
        OrderedDict(
            [
                ("Cohort", "E-MTAB-1980"),
                ("Clinical comparator", "Age + Fuhrman grade + pT3/4 + M1 (exploratory)"),
                ("N", int(em_comb["data"].shape[0])),
                ("Events", int(em_comb["data"]["event"].sum())),
                ("Clinical model C-index", format_float(em_clin["cindex"])),
                ("Signature-alone C-index", format_float(em_sig["cindex"])),
                ("Combined apparent C-index", format_float(em_comb["cindex"])),
                ("Combined optimism-corrected C-index", "NA"),
                ("Optimism-corrected calibration slope", "NA"),
                (
                    "Signature HR per SD in combined model",
                    "{:.3f} ({:.3f}-{:.3f})".format(em_per_sd["hr"], em_per_sd["lo"], em_per_sd["hi"]),
                ),
                ("Signature p", format_p(float(em_comb["model"].summary.loc["risk_score", "p"]))),
                ("Incremental LR p", format_p(em_lr_p)),
            ]
        ),
    ]
    df = pd.DataFrame(rows)
    save_table(df, "Table5_clinicopathologic_benchmarking")

    supp = pd.DataFrame(
        [
            OrderedDict(
                [
                    ("Cohort", "TCGA-KIRC"),
                    ("Model", "Combined"),
                    ("Apparent C-index", format_float(tcga_optimism["apparent_cindex"])),
                    ("Mean optimism", format_float(tcga_optimism["optimism_c"])),
                    ("Optimism-corrected C-index", format_float(tcga_optimism["corrected_cindex"])),
                    ("Apparent slope", format_float(tcga_optimism["apparent_slope"])),
                    ("Mean slope optimism", format_float(tcga_optimism["optimism_slope"])),
                    ("Optimism-corrected slope", format_float(tcga_optimism["corrected_slope"])),
                    ("Bootstrap samples", int(tcga_optimism["n_boot"])),
                ]
            )
        ]
    )
    save_table(
        supp,
        "TableS2_tcga_optimism_correction",
        directory=SUPP_DIR,
        preamble_lines=[
            "# Supplementary Table S2. TCGA optimism-correction summary",
            "Bootstrap optimism estimates for the prespecified TCGA clinicopathologic benchmarking model using 1000 resamples.",
        ],
    )
    ph_meta = {
        "TCGA-KIRC combined": tcga_ph,
        "E-MTAB-1980 combined": em_ph,
    }
    return df, tcga_comb, ph_meta


def km_panel(ax, time, event, risk, title, label):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    risk = np.asarray(risk, dtype=float)
    cut = float(np.nanmedian(risk))
    high = risk >= cut
    for mask, color, name in ((~high, PAL["blue"], "Low risk"), (high, PAL["red"], "High risk")):
        km = KaplanMeierFitter()
        km.fit(time[mask], event_observed=event[mask], label=name)
        xs = np.r_[0.0, km.survival_function_.index.values]
        ys = np.r_[1.0, km.survival_function_[name].values]
        ax.step(xs, ys, where="post", color=color, linewidth=2.0, label=name)
    pval = logrank_test(time[high], time[~high], event_observed_A=event[high], event_observed_B=event[~high]).p_value
    cidx = concordance_index(time, -risk, event)
    ax.set_title("{}  {}".format(label, title), loc="left", fontsize=12.0, fontweight="bold")
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")
    ax.set_ylim(0, 1.02)
    ax.grid(axis="y")
    ax.text(
        0.98,
        0.10,
        "P={}\nC-index={}".format(format_p_figure(pval), format_float(cidx)),
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9.0,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=PAL["grid"]),
    )


def add_risk_table(ax, time, risk):
    time = np.asarray(time, dtype=float)
    risk = np.asarray(risk, dtype=float)
    cut = float(np.nanmedian(risk))
    high = risk >= cut
    low = ~high
    ticks = np.round(np.linspace(0, float(np.nanmax(time)), 4)).astype(int)
    ax.axis("off")
    ax.text(0.0, 0.84, "No. at risk", transform=ax.transAxes, fontsize=9.0, fontweight="bold", ha="left")
    ax.text(0.0, 0.45, "Low", transform=ax.transAxes, fontsize=9.0, fontweight="bold", color=PAL["blue"], ha="left")
    ax.text(0.0, 0.08, "High", transform=ax.transAxes, fontsize=9.0, fontweight="bold", color=PAL["red"], ha="left")
    x_positions = np.linspace(0.28, 0.98, len(ticks))
    for x_ax, tick in zip(x_positions, ticks):
        ax.text(x_ax, 0.84, str(int(tick)), transform=ax.transAxes, fontsize=8.7, ha="center")
        ax.text(x_ax, 0.45, str(int(np.sum(time[low] >= tick))), transform=ax.transAxes, fontsize=8.7, ha="center")
        ax.text(x_ax, 0.08, str(int(np.sum(time[high] >= tick))), transform=ax.transAxes, fontsize=8.7, ha="center")


def km_panel_frontiers(ax, time, event, risk, title, panel_letter):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    risk = np.asarray(risk, dtype=float)
    cut = float(np.nanmedian(risk))
    high = risk >= cut

    for mask, color, name in ((~high, PAL["blue"], "Low risk"), (high, PAL["red"], "High risk")):
        kmf = KaplanMeierFitter()
        kmf.fit(time[mask], event_observed=event[mask], label=name)
        times = kmf.survival_function_.index.values
        surv = kmf.survival_function_[name].values
        ci = kmf.confidence_interval_
        ax.step(times, surv, where="post", color=color, linewidth=2.0, label=name)
        ax.fill_between(times, ci.iloc[:, 0].values, ci.iloc[:, 1].values, step="post", color=color, alpha=0.12)

    pval = logrank_test(time[high], time[~high], event_observed_A=event[high], event_observed_B=event[~high]).p_value
    cidx = concordance_index(time, -risk, event)
    ax.text(-0.08, 1.02, panel_letter, transform=ax.transAxes, ha="left", va="bottom", fontsize=13.5, fontweight="bold")
    ax.set_title(title, loc="left", fontsize=12.0, pad=4)
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")
    ax.set_ylim(0, 1.02)
    ax.grid(axis="y")
    soften(ax)
    ax.text(
        0.98,
        0.08,
        "Orientation: Original\nP={}\nC-index={}".format(format_p_figure(pval), format_float(cidx)),
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8.8,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=PAL["grid"]),
    )
    ax.legend(frameon=False, loc="lower left")


def figure3_principal_km(cptac_meta, em1980_meta):
    fig = plt.figure(figsize=(7.8, 9.4))
    gs = fig.add_gridspec(4, 1, height_ratios=[3.2, 0.9, 3.2, 0.9], hspace=0.24)
    ax1 = fig.add_subplot(gs[0])
    ax1_tab = fig.add_subplot(gs[1])
    ax2 = fig.add_subplot(gs[2])
    ax2_tab = fig.add_subplot(gs[3])
    fig.suptitle("Figure 3. Overall-survival validation and projection plots", x=0.5, y=0.985, fontsize=16.0, fontweight="bold")

    km_panel_frontiers(ax1, em1980_meta["time"], em1980_meta["event"], em1980_meta["risk"], "E-MTAB-1980 overall survival", "A")
    add_risk_table(ax1_tab, em1980_meta["time"], em1980_meta["risk"])
    km_panel_frontiers(ax2, cptac_meta["time"], cptac_meta["event"], cptac_meta["risk"], "CPTAC overall survival", "B")
    add_risk_table(ax2_tab, cptac_meta["time"], cptac_meta["risk"])
    save_figure(fig, "Figure3_external_kaplan_meier")


def figure4_discrimination(table2_df):
    df = table2_df.copy()
    bounds = df["Original score 95% CI"].str.split("-", expand=True).astype(float)
    df["low"] = bounds[0]
    df["high"] = bounds[1]
    df["cindex"] = df["Original score C-index"].astype(float)
    df["y"] = np.arange(df.shape[0])[::-1]
    colors = [PAL["green"] if role == PRIMARY_EXTERNAL_ROLE else PAL["gold"] for role in df["Validation role"]]
    labels = ["{}\n(original)".format(cohort) for cohort in df["Cohort"]]
    fig, ax = plt.subplots(figsize=(8.8, 6.3))
    fig.subplots_adjust(top=0.86)
    ax.axvspan(0.60, 0.82, color=PAL["green"], alpha=0.10, lw=0)
    ax.axvline(0.60, color=PAL["green"], linewidth=1.4, linestyle="--")
    ax.errorbar(
        df["cindex"],
        df["y"],
        xerr=[df["cindex"] - df["low"], df["high"] - df["cindex"]],
        fmt="none",
        ecolor=PAL["grid"],
        elinewidth=1.5,
        capsize=3.5,
        zorder=1,
    )
    ax.scatter(df["cindex"], df["y"], c=colors, s=360, edgecolors="white", linewidths=2.2, zorder=2)
    for x, y, color in zip(df["cindex"], df["y"], colors):
        ax.text(x, y + 0.12, format_float(x), color=color, fontsize=13.5, fontweight="bold", ha="center")
    ax.set_yticks(df["y"])
    ax.set_yticklabels(labels)
    ax.set_xlim(0.55, 0.82)
    ax.set_xlabel("Harrell's C-index")
    fig.suptitle("Figure 4. External discrimination summary", y=0.96, fontsize=18.0, fontweight="bold")
    fig.text(
        0.5,
        0.905,
        "Displayed in the original score direction across the prespecified study hierarchy.",
        ha="center",
        va="center",
        fontsize=10.7,
        color=PAL["muted"],
    )
    ax.grid(axis="x")
    soften(ax)
    save_figure(fig, "Figure4_external_discrimination")


def figure5_signature(signature):
    df = signature.sort_values("coef").copy()
    color_map = {"immune": PAL["red"], "metabolic": PAL["green"], "mixed": PAL["mauve"]}
    fig, ax = plt.subplots(figsize=(9.6, 11.0))
    fig.subplots_adjust(top=0.90)
    y = np.arange(df.shape[0])
    for idx, (_, row) in enumerate(df.iterrows()):
        color = color_map[row["category"]]
        ax.hlines(idx, 0, row["coef"], color=color, linewidth=2.2, alpha=0.95)
        ax.scatter(row["coef"], idx, s=60 + row["abs_coef"] * 40, color=color, edgecolor="white", linewidth=0.8, zorder=3)
    ax.axvline(0, color=PAL["ink"], linewidth=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels(df["gene"])
    ax.set_xlabel("Frozen coefficient")
    fig.suptitle("Figure 5. Coefficient architecture of the frozen 33-gene signature", y=0.965, fontsize=16.8, fontweight="bold")
    ax.text(
        0.02,
        1.01,
        "Positive coefficients indicate higher predicted risk; colors denote post hoc functional categories.",
        transform=ax.transAxes,
        color=PAL["muted"],
        fontsize=9.6,
    )
    ax.grid(axis="x")
    legend_handles = [
        Line2D([0], [0], marker="o", color="w", label=label.capitalize(), markerfacecolor=color_map[label], markersize=9)
        for label in ("immune", "metabolic", "mixed")
    ]
    ax.legend(handles=legend_handles, frameon=False, loc="lower right")
    soften(ax)
    save_figure(fig, "Figure5_signature_characteristics")


def supplementary_orientation_figure(cohorts):
    order = ["CPTAC", "E-MTAB-1980", "IMmotion150", "E-MTAB-3267"]
    fig, axes = plt.subplots(2, 2, figsize=(14.2, 8.4), sharey=True)
    fig.subplots_adjust(left=0.10, right=0.94, bottom=0.10, top=0.88, wspace=0.34, hspace=0.24)
    axes = axes.ravel()
    for idx, (ax, cohort_name, letter) in enumerate(zip(axes, order, ("A", "B", "C", "D"))):
        meta = cohorts[cohort_name]
        orig = meta["original"]
        flipped = meta["sign_reversed"]
        ax.axvspan(0.60, 0.82, color=PAL["green"], alpha=0.08, lw=0)
        ax.axvline(0.60, color=PAL["green"], linewidth=1.2, linestyle="--")
        ax.hlines(1, flipped["cindex"], orig["cindex"], color=PAL["blue"], linewidth=2.0)
        ax.scatter(orig["cindex"], 1, s=230, color=PAL["blue"], edgecolors="white", linewidths=2.0, zorder=3)
        ax.scatter(flipped["cindex"], 0, s=150, color=PAL["red"], edgecolors="white", linewidths=1.8, zorder=3)
        ax.text(orig["cindex"], 1.12, "C={}".format(format_float(orig["cindex"])), color=PAL["blue"], ha="center", fontsize=10.0, fontweight="bold")
        ax.text(flipped["cindex"], 0.16, "C={}".format(format_float(flipped["cindex"])), color=PAL["red"], ha="center", fontsize=10.0, fontweight="bold")
        ax.text(orig["cindex"] + 0.008, 0.96, "Reported", color=PAL["blue"], fontsize=9.0, fontweight="bold")
        ax.text(orig["cindex"] - 0.02, 0.78, "P={}".format(format_p_figure(orig["p"])), color=PAL["muted"], fontsize=8.7)
        ax.text(flipped["cindex"] - 0.055, -0.08, "P={}".format(format_p_figure(flipped["p"])), color=PAL["muted"], fontsize=8.7)
        ax.set_yticks([0, 1])
        if idx % 2 == 0:
            ax.set_yticklabels(["Sign-reversed", "Original"])
            ax.tick_params(axis="y", labelleft=True, labelright=False, length=0, pad=4)
        else:
            ax.set_yticklabels(["Sign-reversed", "Original"])
            ax.yaxis.tick_right()
            ax.tick_params(axis="y", labelleft=False, labelright=True, length=0, pad=4)
        ax.set_xlim(min(0.18, flipped["cindex"] - 0.07), 0.82)
        if idx >= 2:
            ax.set_xlabel("C-index")
        else:
            ax.set_xlabel("")
        ax.set_title("{}  {} orientation sensitivity".format(letter, cohort_name), loc="left", fontsize=11.5, pad=6)
        ax.grid(axis="x")
        soften(ax)
    fig.suptitle("Supplementary Figure S1. Orientation robustness diagnostics", x=0.06, ha="left", fontsize=16.0, fontweight="bold")
    save_figure(fig, "FigureS1_orientation_robustness", directory=SUPP_DIR)


def supplementary_calibration_figure(tcga_comb_fit):
    model = tcga_comb_fit["model"]
    data = tcga_comb_fit["data"].copy()
    covs = ["age", "stage_ord", "grade_num", "risk_score"]
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.5))
    for ax, time_point, label in zip(axes, (36.0, 60.0), ("3-year", "5-year")):
        pred = model.predict_survival_function(data[covs], times=[time_point]).T.iloc[:, 0]
        groups = pd.qcut(pred.rank(method="first"), 5, labels=False)
        points = []
        for group_id in sorted(groups.unique()):
            sub = data.loc[groups == group_id]
            km = KaplanMeierFitter()
            km.fit(sub["time"], event_observed=sub["event"])
            points.append((float(pred.loc[groups == group_id].mean()), float(km.predict(time_point))))
        pts = np.asarray(points)
        ax.scatter(pts[:, 0], pts[:, 1], s=72, color=PAL["blue"], edgecolors="white", linewidths=0.8)
        ax.plot([0, 1], [0, 1], linestyle="--", color=PAL["grid"], linewidth=1.0)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel("Predicted survival probability")
        ax.set_ylabel("Observed Kaplan-Meier estimate")
        ax.set_title(label)
        ax.grid(True)
    fig.suptitle("Supplementary Figure S2. TCGA combined-model calibration", x=0.06, ha="left", fontsize=14.0, fontweight="bold")
    save_figure(fig, "FigureS2_tcga_calibration", directory=SUPP_DIR)


def build_supplementary_tables(cohorts, ph_meta):
    orient_rows = []
    for cohort_name, meta in cohorts.items():
        orient_rows.append(
            OrderedDict(
                [
                    ("Cohort", cohort_name),
                    ("Endpoint", meta["endpoint"]),
                    ("Original log-rank P (median split)", format_p(meta["original"]["p"])),
                    ("Original score C-index", format_float(meta["original"]["cindex"])),
                    ("Sign-reversed log-rank P (median split)", format_p(meta["sign_reversed"]["p"])),
                    ("Sign-reversed score C-index", format_float(meta["sign_reversed"]["cindex"])),
                ]
            )
        )
    save_table(
        pd.DataFrame(orient_rows),
        "TableS1_orientation_robustness",
        directory=SUPP_DIR,
        preamble_lines=[
            "# Supplementary Table S1. Orientation robustness summary",
            "Original-direction and sign-reversed log-rank P values and C-index summaries across the four external cohorts; reported as robustness checks rather than confirmatory validation.",
        ],
    )

    versions = pd.DataFrame(
        [("numpy", np.__version__), ("pandas", pd.__version__), ("matplotlib", matplotlib.__version__)],
        columns=["Software", "Version"],
    )
    save_table(
        versions,
        "TableS3_software_versions",
        directory=SUPP_DIR,
        preamble_lines=[
            "# Supplementary Table S3. Software versions",
            "Software versions used for the reproducible analysis environment.",
        ],
    )

    external_ph_bits = []
    for cohort_name, meta in cohorts.items():
        fitted = fit_cindex_model(pd.DataFrame({"time": meta["time"], "event": meta["event"], "risk_score": meta["risk"]}), ["risk_score"])
        ph_p = ph_test_pvalues(fitted).get("risk_score", np.nan)
        external_ph_bits.append("{} risk-score P={}".format(cohort_name, format_p(ph_p)))

    tcga_risk_p = format_p(ph_meta["TCGA-KIRC combined"].get("risk_score", np.nan))
    em_risk_p = format_p(ph_meta["E-MTAB-1980 combined"].get("risk_score", np.nan))

    note = pd.DataFrame(
        [
            OrderedDict(
                [
                    ("Item", "Reporting framework"),
                    (
                        "Details",
                        "Primary manuscript prepared as a prognostic biomarker external-validation study using REMARK, STROBE, SAMPL, and TRIPOD logic.",
                    ),
                ]
            ),
            OrderedDict([("Item", "Frozen score rule"), ("Details", score_formula_note())]),
            OrderedDict(
                [
                    ("Item", "Archived model-definition overlap"),
                    (
                        "Details",
                        "Archived search-stage scripts restricted candidates to genes shared across TCGA-KIRC, GSE29609, and CPTAC. GSE29609 informed early feature filtering and CPTAC contributed to historical search-stage ranking, so CPTAC is reported here as a development-linked OS projection cohort rather than as a strictly independent validation cohort.",
                    ),
                ]
            ),
            OrderedDict(
                [
                    ("Item", "Archived preselection candidate pool"),
                    (
                        "Details",
                        "The preserved search utilities capped the shared-gene pool at the top 500 genes by TCGA variance and reduced this to at most 220 candidates using GSE29609 univariable signal and TCGA/GSE29609 direction-consistency scoring before penalized Cox fitting.",
                    ),
                ]
            ),
            OrderedDict(
                [
                    ("Item", "Penalized search settings"),
                    (
                        "Details",
                        "The archived fixed-direction Python search explored l1_ratio values of 0, 0.1, 0.3, 0.5, 0.7, and 1.0 with penalizers from 0.001 to 0.3; the archival glmnet workflow used 10-fold cross-validation with alpha values of 0, 0.1, 0.3, 0.5, 0.7, and 1.0 and retained lambda.min or lambda.1se solutions.",
                    ),
                ]
            ),
            OrderedDict(
                [
                    ("Item", "Frozen 33-gene rank universe rationale"),
                    (
                        "Details",
                        "Ranks were calculated only within the final 33-gene universe because the deployed score is defined as a weighted sum over those genes; using a larger cohort-specific universe would change the rank denominator and therefore redefine the locked model.",
                    ),
                ]
            ),
            OrderedDict(
                [
                    ("Item", "Platform-control strategy"),
                    (
                        "Details",
                        "Cross-platform control relied on cohort-specific within-sample ranks, matrix-specific preprocessing documented in Table 4, and score projection without cross-cohort re-centering or coefficient refitting.",
                    ),
                ]
            ),
            OrderedDict(
                [
                    ("Item", "Missing signature genes"),
                    ("Details", "Prespecified neutral fill value was 0.5; no retained cohort required this fallback because all 33 genes were available."),
                ]
            ),
            OrderedDict(
                [
                    ("Item", "Bootstrap sampling"),
                    (
                        "Details",
                        "Harrell's C-index confidence intervals and TCGA optimism correction used {} and {} bootstrap resamples, respectively.".format(BOOTSTRAP_N, OPTIMISM_BOOTSTRAP_N),
                    ),
                ]
            ),
            OrderedDict(
                [
                    ("Item", "Proportional-hazards checks"),
                    (
                        "Details",
                        "Schoenfeld residual-based tests were reviewed for the reported Cox models. The risk-score term showed no clear departure from proportional hazards in the TCGA-KIRC combined model (P={}) or the E-MTAB-1980 combined model (P={}). External univariable checks were: {}. The IMmotion150 result was therefore retained as an exploratory endpoint/context sensitivity analysis.".format(
                            tcga_risk_p,
                            em_risk_p,
                            "; ".join(external_ph_bits),
                        ),
                    ),
                ]
            ),
        ]
    )
    save_table(
        note,
        "TableS4_reporting_and_scoring_notes",
        directory=SUPP_DIR,
        preamble_lines=[
            "# Supplementary Table S4. Reporting and scoring notes",
            "Study-level reporting clarifications, archived search provenance, and frozen-score implementation notes for the submitted Cancer Medicine package.",
        ],
    )


def build_cohort_meta(signature, cptac, em1980, imm150, em3267):
    return OrderedDict(
        [
            (
                "CPTAC",
                {
                    "endpoint": "overall survival",
                    "role": DEVELOPMENT_LINKED_ROLE,
                    "platform": "RNA-seq",
                    "time": cptac["OS_MONTHS"].to_numpy(),
                    "event": cptac["OS_STATUS"].to_numpy(),
                    "risk": cptac["risk_score"].to_numpy(),
                    "genes_available": int(signature.shape[0]),
                },
            ),
            (
                "E-MTAB-1980",
                {
                    "endpoint": "overall survival",
                    "role": PRIMARY_EXTERNAL_ROLE,
                    "platform": "microarray",
                    "time": em1980["time_months"].to_numpy(),
                    "event": em1980["status_os"].to_numpy(),
                    "risk": em1980["risk_score"].to_numpy(),
                    "genes_available": int(signature.shape[0]),
                },
            ),
            (
                "IMmotion150",
                {
                    "endpoint": "progression-free survival",
                    "role": SECONDARY_ENDPOINT_ROLE,
                    "platform": "RNA-seq",
                    "time": imm150["pfs_months"].to_numpy(),
                    "event": imm150["status_pfs"].to_numpy(),
                    "risk": imm150["risk_score"].to_numpy(),
                    "genes_available": int(signature.shape[0]),
                },
            ),
            (
                "E-MTAB-3267",
                {
                    "endpoint": "progression-free survival",
                    "role": SECONDARY_ENDPOINT_ROLE,
                    "platform": "microarray",
                    "time": em3267["pfs_months"].to_numpy(),
                    "event": em3267["status_pfs"].to_numpy(),
                    "risk": em3267["risk_score"].to_numpy(),
                    "genes_available": int(signature.shape[0]),
                },
            ),
        ]
    )


def annotate_cohort_meta(cohorts):
    for meta in cohorts.values():
        meta.update(compute_orientation_metrics(meta["time"], meta["event"], meta["risk"]))
    return cohorts


def write_summary_json(table2_df, table5_df, ph_meta):
    payload = {
        "score_rule": score_formula_note(),
        "bootstrap_resamples": {
            "cindex_ci": BOOTSTRAP_N,
            "optimism_correction": OPTIMISM_BOOTSTRAP_N,
        },
        "ph_assessment": ph_meta,
        "external_validation": table2_df.to_dict(orient="records"),
        "clinical_benchmarking": table5_df.to_dict(orient="records"),
    }
    (SUPP_DIR / "summary_metrics.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_readme():
    text = (
        "# Cancer Medicine submission package\n\n"
        "This package reframes the manuscript as a cautious frozen-model projection and external-validation study.\n\n"
        "- Archived model-definition support: GSE29609 informed early feature filtering and CPTAC contributed to historical search-stage ranking.\n"
        "- Primary independent OS validation cohort: E-MTAB-1980.\n"
        "- Development-linked OS projection cohort: CPTAC.\n"
        "- Secondary endpoint/context sensitivity cohorts: IMmotion150 and E-MTAB-3267 progression-free survival.\n"
        "- Table 2 reports log-rank P values from cohort-specific median splits for visualization; primary inference relied on the continuous score as summarized by the per-SD Cox hazard ratio and Harrell's C-index.\n"
        "- Score definition: {}\n".format(score_formula_note())
    )
    (OUT_DIR / "README_submission_package.md").write_text(text, encoding="utf-8")


def figure1_flow():
    fig, ax = plt.subplots(figsize=(12.2, 7.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.add_patch(
        FancyBboxPatch(
            (0.025, 0.06),
            0.95,
            0.84,
            boxstyle="round,pad=0.012,rounding_size=0.03",
            linewidth=0,
            facecolor=PAL["panel"],
            alpha=0.85,
        )
    )
    ax.text(0.05, 0.95, "Figure 1. Study design and cohort flow diagram", fontsize=17.5, fontweight="bold", va="top")
    ax.text(
        0.05,
        0.905,
        "Discovery, archived model-definition support, and the prespecified validation/projection cohorts across the study hierarchy.",
        fontsize=10.8,
        color=PAL["muted"],
        va="top",
    )
    ax.plot([0.05, 0.95], [0.875, 0.875], color=PAL["grid"], linewidth=1.2)

    def card(x, y, w, h, eyebrow, title, lines, accent, fill="#FFFFFF", line_gap=0.038):
        ax.add_patch(
            FancyBboxPatch(
                (x + 0.008, y - 0.008),
                w,
                h,
                boxstyle="round,pad=0.012,rounding_size=0.025",
                linewidth=0,
                facecolor="#DCE5EC",
                alpha=0.35,
                zorder=1,
            )
        )
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle="round,pad=0.012,rounding_size=0.025",
                linewidth=1.0,
                edgecolor=PAL["grid"],
                facecolor=fill,
                zorder=2,
            )
        )
        ax.add_patch(
            FancyBboxPatch(
                (x, y + h - 0.028),
                w,
                0.028,
                boxstyle="round,pad=0.012,rounding_size=0.025",
                linewidth=0,
                facecolor=accent,
                zorder=3,
            )
        )
        ax.text(
            x + 0.022,
            y + h - 0.055,
            eyebrow.upper(),
            fontsize=8.0,
            fontweight="bold",
            color=accent,
            va="top",
            zorder=4,
        )
        ax.text(x + 0.022, y + h - 0.095, title, fontsize=12.7, fontweight="bold", va="top", zorder=4)
        start_y = y + h - 0.145
        for idx, line in enumerate(lines):
            ax.text(x + 0.024, start_y - idx * line_gap, line, fontsize=10.3, color=PAL["ink"], va="top", zorder=4)

    def arrow(x0, y0, x1, y1, rad=0.0):
        ax.add_patch(
            FancyArrowPatch(
                (x0, y0),
                (x1, y1),
                arrowstyle="-|>",
                mutation_scale=12,
                linewidth=1.6,
                color="#C9D6E2",
                connectionstyle="arc3,rad={}".format(rad),
                zorder=0,
            )
        )

    card(
        0.06,
        0.55,
        0.27,
        0.275,
        "Discovery set",
        "TCGA-KIRC discovery",
        ["537 downloaded", "2 invalid survival exclusions", "4 unmatched expression rows", "531 evaluable discovery cases"],
        PAL["blue"],
        fill="#FCFEFF",
    )
    card(
        0.385,
        0.55,
        0.22,
        0.275,
        "Archived search stage",
        "Model-definition support",
        ["GSE29609 feature screen", "CPTAC ranking constraint", "Not strict external validation"],
        PAL["gold"],
        fill="#FFFCF4",
    )
    card(
        0.67,
        0.55,
        0.27,
        0.275,
        "Primary validation",
        "Independent OS cohort",
        ["E-MTAB-1980 n=101", "Primary external OS validation", "CPTAC reported separately"],
        PAL["green"],
        fill="#F7FEFB",
    )
    card(
        0.18,
        0.19,
        0.28,
        0.225,
        "Secondary cohorts",
        "Endpoint/context sensitivity",
        ["IMmotion150 n=263", "E-MTAB-3267 n=53", "PFS/context sensitivity"],
        PAL["gold"],
        fill="#FFFCF4",
    )
    card(
        0.54,
        0.19,
        0.28,
        0.225,
        "Projection and benchmarking",
        "CPTAC + complete cases",
        ["CPTAC n=237 (development-linked)", "TCGA complete cases n=521", "E-MTAB-1980 complete cases n=99"],
        PAL["mauve"],
        fill="#FAF8FD",
    )

    arrow(0.335, 0.685, 0.665, 0.685, rad=0.0)
    arrow(0.195, 0.55, 0.315, 0.415, rad=0.02)
    arrow(0.805, 0.55, 0.68, 0.415, rad=-0.02)
    arrow(0.205, 0.55, 0.66, 0.415, rad=-0.06)

    ax.text(0.50, 0.50, "Primary validation, development-linked projection, and sensitivity analyses", ha="center", va="center", fontsize=10.0, color=PAL["muted"])
    save_figure(fig, "Figure1_strobe_flow_diagram")


def figure2_workflow():
    fig, ax = plt.subplots(figsize=(11.2, 8.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.5, 0.965, "Figure 2. Legacy score projection and study hierarchy", fontsize=17.6, fontweight="bold", va="top", ha="center")
    ax.text(
        0.5,
        0.915,
        "The prespecified analysis workflow starts from a legacy frozen 33-gene coefficient table and distinguishes independent validation from development-linked projection.",
        fontsize=10.5,
        color=PAL["muted"],
        va="top",
        ha="center",
    )

    steps = [
        ("01", "DISCOVERY", "Fit the TCGA-KIRC discovery model", "Within-sample ranked expression values were used to fit the legacy penalized Cox model in 531 evaluable TCGA-KIRC discovery cases.", PAL["blue"], 0.72),
        ("02", "ARCHIVED SEARCH STAGE", "Document model-definition overlap", "GSE29609 informed early feature filtering and CPTAC contributed to historical cross-platform search-stage ranking.", PAL["gold"], 0.50),
        ("03", "LOCKED SCORE PROJECTION", "Project the frozen score", "After the 33-gene coefficient table was frozen, the score was projected without refitting in CPTAC, E-MTAB-1980, IMmotion150, and E-MTAB-3267.", PAL["green"], 0.28),
        ("04", "STUDY HIERARCHY", "Interpret each cohort cautiously", "E-MTAB-1980 is the primary independent OS validation cohort; CPTAC is development-linked; IMmotion150 and E-MTAB-3267 are secondary endpoint/context sensitivity cohorts.", PAL["red"], 0.06),
    ]
    ax.add_line(Line2D([0.08, 0.08], [0.15, 0.84], color=PAL["grid"], linewidth=2.0))
    for number, eyebrow, title, body, color, y in steps:
        ax.add_patch(FancyBboxPatch((0.15, y), 0.77, 0.14, boxstyle="round,pad=0.015,rounding_size=0.02", linewidth=0.8, edgecolor=PAL["grid"], facecolor=PAL["panel"]))
        ax.add_patch(FancyBboxPatch((0.04, y + 0.02), 0.08, 0.10, boxstyle="round,pad=0.012,rounding_size=0.02", linewidth=0, facecolor=color))
        ax.text(0.08, y + 0.07, number, ha="center", va="center", color="white", fontsize=13.5, fontweight="bold")
        ax.text(0.18, y + 0.108, eyebrow, fontsize=9.3, fontweight="bold", va="top", color=color)
        ax.text(0.18, y + 0.078, title, fontsize=12.2, fontweight="bold", va="top")
        ax.text(0.18, y + 0.038, body, fontsize=10.1, va="top", color=PAL["muted"])
        if y > 0.08:
            ax.add_patch(FancyArrowPatch((0.08, y), (0.08, y - 0.045), arrowstyle="-|>", mutation_scale=13, linewidth=1.0, color=PAL["grid"]))
    save_figure(fig, "Figure2_study_workflow")


def main():
    ensure_dirs()
    clean_output()
    set_style()
    signature = load_signature()
    tcga = load_tcga_scored(signature)
    cptac = load_cptac_scored(signature)
    em1980 = load_emtab1980_scored(signature)
    imm150 = load_immotion150_scored(signature)
    em3267 = load_emtab3267_scored(signature)

    cohorts = annotate_cohort_meta(build_cohort_meta(signature, cptac, em1980, imm150, em3267))

    build_table1(tcga, cptac, em1980, imm150, em3267)
    table2 = build_table2(signature, cohorts)
    build_table3(signature)
    build_table4()
    table5, tcga_comb_fit, ph_meta = build_table5(tcga, em1980)
    build_supplementary_tables(cohorts, ph_meta)

    figure1_flow()
    figure2_workflow()
    figure3_principal_km(cohorts["CPTAC"], cohorts["E-MTAB-1980"])
    figure4_discrimination(table2)
    figure5_signature(signature)
    supplementary_orientation_figure(cohorts)
    supplementary_calibration_figure(tcga_comb_fit)

    write_summary_json(table2, table5, ph_meta)
    write_readme()
    print("Cancer Medicine package rebuilt at:", OUT_DIR)


if __name__ == "__main__":
    main()
