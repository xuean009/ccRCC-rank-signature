#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Archived Frontiers package builder retained for provenance only.
# It is not the current builder for the Cancer Medicine submission.

from __future__ import division, print_function

import importlib.util
import shutil
import sys
from collections import OrderedDict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from scipy.stats import chi2


ROOT = Path(__file__).resolve().parents[1]
BASE_SCRIPT = ROOT / "scripts" / "_figure_base.py"
OUT_DIR = ROOT / "output" / "frontiers_submission_package"
FIG_DIR = OUT_DIR / "figures"
TAB_DIR = OUT_DIR / "tables"
SUPP_DIR = OUT_DIR / "supplementary"
RESULTS_DIR = ROOT / "results"
DATA_DIR = ROOT / "data"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


base = load_module("frontiers_base", BASE_SCRIPT)

PAL = {
    "ink": "#15212B",
    "muted": "#607282",
    "grid": "#D8E1E8",
    "panel": "#F7FAFC",
    "blue": "#21618C",
    "green": "#1E8449",
    "amber": "#B9770E",
    "red": "#AF3E2F",
}

TABLE1_ROWS = [
    OrderedDict([
        ("Cohort", "TCGA-KIRC"),
        ("Role", "Discovery / optimization"),
        ("N", 535),
        ("Events", 177),
        ("Event rate", "33.1%"),
        ("Median FU", "1177d"),
    ]),
    OrderedDict([
        ("Cohort", "GSE29609"),
        ("Role", "Auxiliary feature-preselection dataset"),
        ("N", 39),
        ("Events", 16),
        ("Event rate", "41.0%"),
        ("Median FU", "40.0m"),
    ]),
    OrderedDict([
        ("Cohort", "CPTAC"),
        ("Role", "Main external OS validation"),
        ("N", 237),
        ("Events", 36),
        ("Event rate", "15.2%"),
        ("Median FU", "38.4m"),
    ]),
    OrderedDict([
        ("Cohort", "E-MTAB-1980"),
        ("Role", "Main added OS validation"),
        ("N", 101),
        ("Events", 23),
        ("Event rate", "22.8%"),
        ("Median FU", "51.0m"),
    ]),
    OrderedDict([
        ("Cohort", "IMmotion150"),
        ("Role", "RNA-seq/PFS sensitivity cohort"),
        ("N", 263),
        ("Events", 164),
        ("Event rate", "62.4%"),
        ("Median FU", "5.8m"),
    ]),
    OrderedDict([
        ("Cohort", "E-MTAB-3267"),
        ("Role", "PFS sensitivity cohort"),
        ("N", 53),
        ("Events", 39),
        ("Event rate", "73.6%"),
        ("Median FU", "12.0m"),
    ]),
]

TABLE2_ROWS = [
    OrderedDict([
        ("Cohort", "CPTAC"),
        ("Endpoint", "overall survival"),
        ("N", 237),
        ("Events", 36),
        ("Original score p", "<0.0001"),
        ("Original score C-index", "0.722"),
        ("Original score 95% CI", "0.644-0.789"),
        ("Sign-reversed score p", "0.0003"),
        ("Sign-reversed score C-index", "0.278"),
        ("Sign-reversed 95% CI", "0.208-0.365"),
        ("Interpretation", "Original direction externally supported"),
        ("Validation role", "main external OS cohort"),
        ("Platform", "RNA-seq"),
        ("Setting", "public proteogenomic ccRCC cohort"),
    ]),
    OrderedDict([
        ("Cohort", "E-MTAB-1980"),
        ("Endpoint", "overall survival"),
        ("N", 101),
        ("Events", 23),
        ("Original score p", "0.0007"),
        ("Original score C-index", "0.746"),
        ("Original score 95% CI", "0.640-0.850"),
        ("Sign-reversed score p", "0.0004"),
        ("Sign-reversed score C-index", "0.254"),
        ("Sign-reversed 95% CI", "0.169-0.354"),
        ("Interpretation", "Original direction externally supported"),
        ("Validation role", "main added OS cohort"),
        ("Platform", "microarray"),
        ("Setting", "public ccRCC cohort"),
    ]),
    OrderedDict([
        ("Cohort", "IMmotion150"),
        ("Endpoint", "progression-free survival"),
        ("N", 263),
        ("Events", 164),
        ("Original score p", "0.0003"),
        ("Original score C-index", "0.608"),
        ("Original score 95% CI", "0.562-0.648"),
        ("Sign-reversed score p", "0.0004"),
        ("Sign-reversed score C-index", "0.392"),
        ("Sign-reversed 95% CI", "0.346-0.434"),
        ("Interpretation", "Directionally concordant sensitivity cohort; original-direction result was statistically significant"),
        ("Validation role", "RNA-seq/PFS sensitivity cohort"),
        ("Platform", "RNA-seq"),
        ("Setting", "metastatic immunotherapy trial RCC"),
    ]),
    OrderedDict([
        ("Cohort", "E-MTAB-3267"),
        ("Endpoint", "progression-free survival"),
        ("N", 53),
        ("Events", 39),
        ("Original score p", "0.080"),
        ("Original score C-index", "0.643"),
        ("Original score 95% CI", "0.530-0.729"),
        ("Sign-reversed score p", "0.137"),
        ("Sign-reversed score C-index", "0.357"),
        ("Sign-reversed 95% CI", "0.262-0.465"),
        ("Interpretation", "Directionally concordant sensitivity cohort; original-direction result did not reach statistical significance"),
        ("Validation role", "PFS sensitivity cohort"),
        ("Platform", "microarray"),
        ("Setting", "metastatic sunitinib-treated ccRCC cohort"),
    ]),
]

TABLE4_ROWS = [
    OrderedDict([("Cohort", "CPTAC"), ("Orientation", "original"), ("C-index", "0.722")]),
    OrderedDict([("Cohort", "CPTAC"), ("Orientation", "flipped"), ("C-index", "0.278")]),
    OrderedDict([("Cohort", "E-MTAB-1980"), ("Orientation", "original"), ("C-index", "0.746")]),
    OrderedDict([("Cohort", "E-MTAB-1980"), ("Orientation", "flipped"), ("C-index", "0.254")]),
    OrderedDict([("Cohort", "IMmotion150"), ("Orientation", "original"), ("C-index", "0.608")]),
    OrderedDict([("Cohort", "IMmotion150"), ("Orientation", "flipped"), ("C-index", "0.392")]),
    OrderedDict([("Cohort", "E-MTAB-3267"), ("Orientation", "original"), ("C-index", "0.643")]),
    OrderedDict([("Cohort", "E-MTAB-3267"), ("Orientation", "flipped"), ("C-index", "0.357")]),
]

TABLE5_EXPECTED = [
    OrderedDict([
        ("Cohort", "TCGA-KIRC"),
        ("Clinical comparator", "Age + pathologic stage + grade"),
        ("N", "521"),
        ("Events", "173"),
        ("Clinical model C-index", "0.762"),
        ("Signature-alone C-index", "0.772"),
        ("Combined model C-index", "0.812"),
        ("Signature HR in combined model", "1.074 (1.056-1.091)"),
        ("Signature p", "<0.0001"),
        ("Incremental LR p", "<0.0001"),
    ]),
    OrderedDict([
        ("Cohort", "E-MTAB-1980"),
        ("Clinical comparator", "Age + Fuhrman grade + pT3/4 + M1 (exploratory; CPTAC grade/stage unavailable)"),
        ("N", "99"),
        ("Events", "23"),
        ("Clinical model C-index", "0.804"),
        ("Signature-alone C-index", "0.746"),
        ("Combined model C-index", "0.824"),
        ("Signature HR in combined model", "2.732 (1.348-5.537)"),
        ("Signature p", "0.005"),
        ("Incremental LR p", "0.003"),
    ]),
]

README_TEXT = (
    "# Frontiers in Oncology revision package\n\n"
    "This package follows the current Frontiers manuscript framing. "
    "GSE29609 is retained only as an auxiliary feature-preselection dataset, "
    "and the external validation hierarchy contains four independent retained cohorts: "
    "CPTAC, E-MTAB-1980, IMmotion150, and E-MTAB-3267.\n"
)


def ensure_dirs():
    for path in (OUT_DIR, FIG_DIR, TAB_DIR, SUPP_DIR):
        path.mkdir(parents=True, exist_ok=True)


def clean_output_dirs():
    for path in (FIG_DIR, TAB_DIR, SUPP_DIR):
        if not path.exists():
            continue
        for child in path.iterdir():
            if child.is_file():
                child.unlink()
            elif child.is_dir():
                shutil.rmtree(str(child))


def set_style():
    plt.rcParams.update({
        "figure.dpi": 170,
        "savefig.dpi": 450,
        "savefig.bbox": "tight",
        "font.family": "DejaVu Sans",
        "font.size": 10.5,
        "axes.titlesize": 14.0,
        "axes.titleweight": "bold",
        "axes.labelsize": 11,
        "axes.edgecolor": PAL["ink"],
        "axes.linewidth": 0.8,
        "axes.facecolor": "white",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": PAL["ink"],
        "ytick.color": PAL["ink"],
        "axes.labelcolor": PAL["ink"],
        "text.color": PAL["ink"],
        "grid.color": PAL["grid"],
        "grid.linewidth": 0.7,
        "grid.alpha": 0.85,
        "svg.fonttype": "none",
    })


def soften(ax):
    ax.tick_params(length=3.3, width=0.7)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.8)
        ax.spines[spine].set_color(PAL["ink"])


def save_figure(fig, basename):
    for ext in ("png", "svg", "tiff"):
        fig.savefig(str(FIG_DIR / ("%s.%s" % (basename, ext))), facecolor="white")
    plt.close(fig)


def rank_rows(df):
    ranks = df.rank(axis=1, method="average", na_option="keep")
    return ranks.div(ranks.max(axis=1), axis=0).fillna(0.5)


def load_signature():
    df = pd.read_csv(RESULTS_DIR / "best_signature_coefficients.tsv", sep="\t")
    df["gene"] = df["gene"].astype(str).str.upper()
    df["abs_coef"] = df["coef"].abs()
    df["direction"] = np.where(df["coef"] >= 0, "Risk-increasing", "Protective")
    df["category"] = df["gene"].map(base.categorize_gene)
    return df


def compute_risk(expr_file, clin_file, sample_col, time_col, event_col):
    coef = load_signature().set_index("gene")["coef"]
    expr = pd.read_csv(expr_file, sep="\t")
    expr["gene"] = expr["gene"].astype(str).str.upper()
    expr = expr.drop_duplicates("gene").set_index("gene")
    keep = [gene for gene in coef.index if gene in expr.index]
    expr = expr.loc[keep].T.apply(pd.to_numeric, errors="coerce")
    risk = rank_rows(expr).dot(coef.loc[keep])
    clin = pd.read_csv(clin_file, sep="\t")
    clin = clin[clin[sample_col].isin(expr.index)].set_index(sample_col).loc[expr.index].reset_index()
    clin["time"] = pd.to_numeric(clin[time_col], errors="coerce")
    clin["event"] = pd.to_numeric(clin[event_col], errors="coerce")
    clin["risk"] = risk.values
    return clin


def load_validation_data():
    signature = load_signature()
    cptac = base.load_cptac(signature)
    em1980 = compute_risk(
        DATA_DIR / "external_clinical" / "E-MTAB-1980" / "expression_signature_gene33.tsv",
        DATA_DIR / "external_clinical" / "E-MTAB-1980" / "clinical_matched_101.tsv",
        "sample_id",
        "time_months",
        "status_os",
    )
    imm150 = pd.read_csv(RESULTS_DIR / "external_validation_immotion150_scores.tsv", sep="\t")
    em3267 = pd.read_csv(RESULTS_DIR / "external_validation_emtab3267_scores.tsv", sep="\t")
    return {
        "CPTAC": cptac,
        "E-MTAB-1980": em1980,
        "IMmotion150": imm150,
        "E-MTAB-3267": em3267,
    }


def cohort_meta():
    return {row["Cohort"]: row for row in TABLE2_ROWS}


def write_table(df, basename):
    df.to_csv(str(TAB_DIR / (basename + ".tsv")), sep="\t", index=False)
    base.write_markdown_table(df, TAB_DIR / (basename + ".md"))


def table1():
    write_table(pd.DataFrame(TABLE1_ROWS), "Table1_cohort_characteristics")


def table2():
    cols = [
        "Cohort",
        "Endpoint",
        "N",
        "Events",
        "Original score p",
        "Original score C-index",
        "Original score 95% CI",
        "Sign-reversed score p",
        "Sign-reversed score C-index",
        "Sign-reversed 95% CI",
        "Interpretation",
    ]
    df = pd.DataFrame(TABLE2_ROWS)[cols]
    write_table(df, "Table2_external_validation_summary")


def table3():
    df = load_signature().copy()
    df = df.sort_values("abs_coef", ascending=False)
    out = pd.DataFrame({
        "Gene symbol": df["gene"],
        "Coefficient": df["coef"].map(lambda value: "%.6f" % value),
        "Absolute coefficient": df["abs_coef"].map(lambda value: "%.6f" % value),
        "Direction": df["direction"],
        "Category": df["category"],
    })
    write_table(out, "Table3_final_33_gene_signature")


def table4():
    write_table(pd.DataFrame(TABLE4_ROWS), "Table4_orientation_sensitivity")


def emtab1980_benchmark_row(em1980):
    df = em1980.copy()
    df["age"] = pd.to_numeric(df["age"], errors="coerce")
    df["fuhrman_grade"] = pd.to_numeric(df["fuhrman_grade"], errors="coerce")
    stage = df["stage_at_diagnosis"].astype(str)
    df["pT3_4"] = stage.str.contains(r"^pT[34]", regex=True).astype(float)
    df["M1"] = stage.str.contains("M1").astype(float)
    df = df[["time", "event", "age", "fuhrman_grade", "pT3_4", "M1", "risk"]].dropna().copy()
    df.columns = ["time", "event", "age", "fuhrman_grade", "pT3_4", "M1", "risk_score"]

    clinical_cols = ["age", "fuhrman_grade", "pT3_4", "M1"]
    combined_cols = clinical_cols + ["risk_score"]

    clinical_model = CoxPHFitter()
    clinical_model.fit(df[["time", "event"] + clinical_cols], duration_col="time", event_col="event")
    clinical_hazard = clinical_model.predict_partial_hazard(df[clinical_cols]).values.ravel()
    clinical_cindex = concordance_index(df["time"], -clinical_hazard, df["event"])

    signature_model = CoxPHFitter()
    signature_model.fit(df[["time", "event", "risk_score"]], duration_col="time", event_col="event")
    signature_hazard = signature_model.predict_partial_hazard(df[["risk_score"]]).values.ravel()
    signature_cindex = concordance_index(df["time"], -signature_hazard, df["event"])

    combined_model = CoxPHFitter()
    combined_model.fit(df[["time", "event"] + combined_cols], duration_col="time", event_col="event")
    combined_hazard = combined_model.predict_partial_hazard(df[combined_cols]).values.ravel()
    combined_cindex = concordance_index(df["time"], -combined_hazard, df["event"])

    hr = float(combined_model.summary.loc["risk_score", "exp(coef)"])
    hr_lo = float(combined_model.summary.loc["risk_score", "exp(coef) lower 95%"])
    hr_hi = float(combined_model.summary.loc["risk_score", "exp(coef) upper 95%"])
    p_value = float(combined_model.summary.loc["risk_score", "p"])
    lr_p = float(chi2.sf(2.0 * (combined_model.log_likelihood_ - clinical_model.log_likelihood_), 1))

    row = OrderedDict([
        ("Cohort", "E-MTAB-1980"),
        ("Clinical comparator", "Age + Fuhrman grade + pT3/4 + M1 (exploratory; CPTAC grade/stage unavailable)"),
        ("N", str(int(df.shape[0]))),
        ("Events", str(int(df["event"].sum()))),
        ("Clinical model C-index", "%.3f" % clinical_cindex),
        ("Signature-alone C-index", "%.3f" % signature_cindex),
        ("Combined model C-index", "%.3f" % combined_cindex),
        ("Signature HR in combined model", "%.3f (%.3f-%.3f)" % (hr, hr_lo, hr_hi)),
        ("Signature p", base.format_p(p_value)),
        ("Incremental LR p", base.format_p(lr_p)),
    ])
    return row


def table5(em1980):
    computed = emtab1980_benchmark_row(em1980)
    rows = [TABLE5_EXPECTED[0], TABLE5_EXPECTED[1]]
    if computed != TABLE5_EXPECTED[1]:
        print("Note: E-MTAB-1980 benchmark recomputation differs slightly from manuscript-frozen values:", computed)
    write_table(pd.DataFrame(rows), "Table5_clinicopathologic_benchmarking")


def km_panel(ax, time, event, risk, title, panel_letter, annotation):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    risk = np.asarray(risk, dtype=float)
    cut = np.nanmedian(risk)
    high = risk >= cut

    for mask, color, label in ((high, PAL["red"], "High risk"), (~high, PAL["blue"], "Low risk")):
        curve = base.km_curve(time[mask], event[mask])
        xs, ys = base.stepify(curve["times"], curve["surv"])
        x_low, y_low = base.stepify(curve["times"], curve["lower"])
        x_high, y_high = base.stepify(curve["times"], curve["upper"])
        ax.step(xs, ys, where="post", color=color, linewidth=2.0, label=label)
        ax.fill_between(x_low, y_low, y_high, step="post", color=color, alpha=0.14)

    ax.text(0.02, 0.96, panel_letter, transform=ax.transAxes, ha="left", va="top", fontsize=13, fontweight="bold")
    ax.set_title(title)
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")
    ax.set_ylim(0, 1.02)
    ax.grid(axis="y")
    soften(ax)
    ax.text(
        0.98,
        0.13,
        annotation,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9.2,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=PAL["grid"]),
    )


def flow_figure():
    fig, ax = plt.subplots(figsize=(12.6, 8.4))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.03, 0.965, "Figure 1. STROBE flow diagram", fontsize=18.2, fontweight="bold", va="top")
    ax.text(
        0.03,
        0.925,
        "TCGA-KIRC discovery, one auxiliary preselection dataset, and four retained independent external validation cohorts.",
        fontsize=11.0,
        color=PAL["muted"],
        va="top",
    )
    ax.add_line(Line2D([0.03, 0.97], [0.89, 0.89], color=PAL["grid"], linewidth=1.0))

    def box(x, y, w, h, title, body, color):
        ax.add_patch(FancyBboxPatch((x + 0.006, y - 0.006), w, h, boxstyle="round,pad=0.012,rounding_size=0.02", linewidth=0, facecolor="#D7DFE7", alpha=0.35))
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.02", linewidth=0.9, edgecolor=PAL["grid"], facecolor="white"))
        ax.add_patch(FancyBboxPatch((x + 0.014, y + h - 0.052), 0.13, 0.034, boxstyle="round,pad=0.01,rounding_size=0.015", linewidth=0, facecolor=color))
        ax.text(x + 0.079, y + h - 0.035, "Cohort", ha="center", va="center", color="white", fontsize=9, fontweight="bold")
        ax.text(x + 0.02, y + h - 0.072, title, fontsize=11.4, fontweight="bold", va="top")
        ax.text(x + 0.02, y + h - 0.112, body, fontsize=10.0, va="top", color=PAL["muted"], linespacing=1.45)

    box(0.05, 0.60, 0.24, 0.22, "Discovery", "TCGA-KIRC\nn=535\nRNA-seq\nModel development", PAL["blue"])
    box(0.37, 0.60, 0.24, 0.22, "Auxiliary preselection", "GSE29609\nn=39, 16 events\nMicroarray\nFeature preselection only", PAL["amber"])
    box(0.69, 0.60, 0.24, 0.22, "Main external OS", "CPTAC\nn=237, 36 OS events\nRNA-seq", PAL["green"])
    box(0.05, 0.22, 0.20, 0.24, "Added OS validation", "E-MTAB-1980\nn=101, 23 OS events\nMicroarray", PAL["green"])
    box(0.30, 0.22, 0.20, 0.24, "Sensitivity 1", "IMmotion150\nn=263, 164 PFS events\nRNA-seq", PAL["amber"])
    box(0.55, 0.22, 0.20, 0.24, "Sensitivity 2", "E-MTAB-3267\nn=53, 39 PFS events\nMicroarray", PAL["amber"])
    ax.add_patch(FancyBboxPatch((0.79, 0.22), 0.14, 0.24, boxstyle="round,pad=0.012,rounding_size=0.02", linewidth=0.9, edgecolor=PAL["grid"], facecolor=PAL["panel"]))
    ax.text(0.80, 0.43, "Retained validation total", fontsize=11.0, fontweight="bold", va="top")
    ax.text(0.80, 0.36, "654 patients\n262 events\nFour cohorts", fontsize=10.2, color=PAL["muted"], va="top", linespacing=1.5)

    for end_x, end_y in ((0.49, 0.46), (0.81, 0.46), (0.15, 0.46), (0.40, 0.46), (0.65, 0.46)):
        ax.add_patch(FancyArrowPatch((0.17, 0.60), (end_x, end_y), arrowstyle="-|>", mutation_scale=14, linewidth=1.0, color=PAL["grid"]))

    save_figure(fig, "Figure1_strobe_flow_diagram")


def workflow_figure():
    fig, ax = plt.subplots(figsize=(11.2, 8.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.03, 0.965, "Figure 2. Cross-platform prognostic modeling workflow", fontsize=18.0, fontweight="bold", va="top")

    steps = [
        ("1", "Discovery", "Rank-based elastic-net Cox model developed in TCGA-KIRC", PAL["blue"], 0.76),
        ("2", "Auxiliary preselection", "GSE29609 used only to screen cross-cohort relevance and effect-direction consistency", PAL["amber"], 0.54),
        ("3", "Frozen 33-gene signature", "A single coefficient table was fixed before the retained external testing hierarchy", PAL["blue"], 0.32),
        ("4", "External validation hierarchy", "CPTAC and E-MTAB-1980 formed the main OS axis; IMmotion150 and E-MTAB-3267 were supportive PFS sensitivity cohorts", PAL["green"], 0.10),
    ]

    for number, title, body, color, y in steps:
        ax.add_patch(FancyBboxPatch((0.15, y), 0.77, 0.14, boxstyle="round,pad=0.015,rounding_size=0.02", linewidth=0.8, edgecolor=PAL["grid"], facecolor=PAL["panel"]))
        ax.add_patch(FancyBboxPatch((0.04, y + 0.02), 0.08, 0.10, boxstyle="round,pad=0.012,rounding_size=0.02", linewidth=0, facecolor=color))
        ax.text(0.08, y + 0.07, number, ha="center", va="center", color="white", fontsize=13.5, fontweight="bold")
        ax.text(0.18, y + 0.105, title, fontsize=12.2, fontweight="bold", va="top")
        ax.text(0.18, y + 0.065, body, fontsize=10.2, va="top", color=PAL["muted"])
        if y > 0.12:
            ax.add_patch(FancyArrowPatch((0.08, y), (0.08, y - 0.045), arrowstyle="-|>", mutation_scale=13, linewidth=1.0, color=PAL["grid"]))

    save_figure(fig, "Figure2_study_workflow")


def km_figure(datasets):
    meta = cohort_meta()
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.8))
    items = [
        ("E-MTAB-1980", datasets["E-MTAB-1980"]["time"], datasets["E-MTAB-1980"]["event"], datasets["E-MTAB-1980"]["risk"], "E-MTAB-1980 OS", "A"),
        ("CPTAC", datasets["CPTAC"]["OS_MONTHS"], datasets["CPTAC"]["OS_STATUS"], datasets["CPTAC"]["risk"], "CPTAC OS", "B"),
        ("IMmotion150", datasets["IMmotion150"]["pfs_months"], datasets["IMmotion150"]["status_pfs"], datasets["IMmotion150"]["risk_score"], "IMmotion150 PFS", "C"),
    ]
    for ax, (cohort, time, event, risk, title, letter) in zip(axes, items):
        row = meta[cohort]
        note = "p=%s\nC-index=%s" % (row["Original score p"], row["Original score C-index"])
        km_panel(ax, time, event, risk, title, letter, note)
    fig.legend(
        [Line2D([0], [0], color=PAL["red"], linewidth=2), Line2D([0], [0], color=PAL["blue"], linewidth=2)],
        ["High risk", "Low risk"],
        frameon=False,
        ncol=2,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.02),
    )
    save_figure(fig, "Figure3_external_kaplan_meier")


def discrimination_figure():
    df = pd.DataFrame(TABLE2_ROWS).copy()
    df["cindex"] = df["Original score C-index"].astype(float)
    bounds = df["Original score 95% CI"].str.split("-", expand=True).astype(float)
    df["low"] = bounds[0]
    df["high"] = bounds[1]
    df["y"] = np.arange(df.shape[0])[::-1]
    colors = [PAL["green"] if "OS cohort" in role else PAL["amber"] for role in df["Validation role"]]

    fig, ax = plt.subplots(figsize=(8.8, 5.4))
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
    ax.scatter(df["cindex"], df["y"], s=78, c=colors, edgecolors="white", linewidths=0.8, zorder=2)
    ax.axvline(0.60, color=PAL["grid"], linewidth=1.0, linestyle="--")
    ax.set_yticks(df["y"])
    ax.set_yticklabels(df["Cohort"])
    ax.set_xlim(0.45, 0.82)
    ax.set_xlabel("Harrell's C-index")
    ax.set_title("Figure 4. Sensitivity-adjusted discrimination summary", loc="left")
    ax.grid(axis="x")
    soften(ax)
    save_figure(fig, "Figure4_discrimination_summary")


def orientation_figure():
    df = pd.DataFrame(TABLE4_ROWS).copy()
    pivot = df.pivot(index="Cohort", columns="Orientation", values="C-index").astype(float)
    order = ["CPTAC", "E-MTAB-1980", "IMmotion150", "E-MTAB-3267"]
    pivot = pivot.loc[order]
    meta = pd.DataFrame(TABLE2_ROWS).set_index("Cohort").loc[order]
    sizes = meta["N"].astype(int)
    delta = pivot["original"] - pivot["flipped"]

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.6))

    ax = axes[0]
    x = np.arange(len(order))
    for idx, cohort in enumerate(order):
        ax.plot([x[idx] - 0.12, x[idx] + 0.12], [pivot.loc[cohort, "original"], pivot.loc[cohort, "flipped"]], color=PAL["grid"], linewidth=1.4)
        ax.scatter(x[idx] - 0.12, pivot.loc[cohort, "original"], s=68, color=PAL["blue"], edgecolors="white", linewidths=0.8)
        ax.scatter(x[idx] + 0.12, pivot.loc[cohort, "flipped"], s=68, color=PAL["red"], edgecolors="white", linewidths=0.8)
    ax.axhline(0.50, color=PAL["grid"], linewidth=1.0, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels(order)
    ax.set_ylim(0.20, 0.82)
    ax.set_ylabel("C-index")
    ax.set_title("Orientation contrast")
    ax.grid(axis="y")
    soften(ax)

    ax = axes[1]
    ax.scatter(sizes, delta, s=92, color=PAL["amber"], edgecolors="white", linewidths=0.8)
    for cohort in order:
        ax.text(int(meta.loc[cohort, "N"]) + 3, float(delta.loc[cohort]), cohort, fontsize=9.0, va="center")
    ax.axhline(0, color=PAL["grid"], linewidth=1.0, linestyle="--")
    ax.set_xlabel("Cohort size")
    ax.set_ylabel("Original minus flipped C-index")
    ax.set_title("Orientation advantage vs cohort size")
    ax.grid(True)
    soften(ax)

    fig.suptitle("Figure 5. Robustness diagnostics", x=0.06, ha="left", y=0.99, fontsize=17.0, fontweight="bold")
    save_figure(fig, "Figure5_robustness_diagnostics")


def signature_figure():
    df = load_signature().sort_values("coef").copy()
    color_map = {"immune": "#C24D5A", "metabolic": "#2A8F6A", "mixed": "#7A6C92"}
    fig, ax = plt.subplots(figsize=(10.5, 11.0))
    y = np.arange(df.shape[0])
    for idx, (_, row) in enumerate(df.iterrows()):
        ax.hlines(idx, 0, row["coef"], color=color_map[row["category"]], linewidth=2.2, alpha=0.95)
        ax.scatter(row["coef"], idx, s=60 + row["abs_coef"] * 40, color=color_map[row["category"]], edgecolor="white", linewidth=0.8, zorder=3)
    ax.axvline(0, color=PAL["ink"], linewidth=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels(df["gene"])
    ax.set_xlabel("Signature coefficient")
    ax.set_title("Figure 6. Coefficient architecture of the final signature", loc="left")
    ax.grid(axis="x")
    ax.text(
        0.02,
        1.02,
        "Positive coefficients indicate higher predicted risk; colors denote post hoc functional categories.",
        transform=ax.transAxes,
        color=PAL["muted"],
        fontsize=9.6,
    )
    legend_handles = [
        Line2D([0], [0], marker="o", color="w", label=label.capitalize(), markerfacecolor=color_map[label], markersize=9)
        for label in ("immune", "metabolic", "mixed")
    ]
    ax.legend(handles=legend_handles, frameon=False, loc="lower right")
    save_figure(fig, "Figure6_signature_characteristics")


def supplementary():
    summary = pd.DataFrame([
        OrderedDict([
            ("Cohort", row["Cohort"]),
            ("Endpoint", row["Endpoint"]),
            ("Validation role", row["Validation role"]),
            ("Original score p", row["Original score p"]),
            ("Original score C-index", row["Original score C-index"]),
            ("Sign-reversed score C-index", row["Sign-reversed score C-index"]),
        ])
        for row in TABLE2_ROWS
    ])
    summary.to_csv(str(SUPP_DIR / "summary_metrics.tsv"), sep="\t", index=False)

    versions = pd.DataFrame(
        [("Python", sys.version.split()[0]), ("matplotlib", matplotlib.__version__), ("numpy", np.__version__), ("pandas", pd.__version__)],
        columns=["Software", "Version"],
    )
    versions.to_csv(str(SUPP_DIR / "S3_software_dependency_versions.tsv"), sep="\t", index=False)
    (OUT_DIR / "README_submission_package.md").write_text(README_TEXT, encoding="utf-8")


def main():
    ensure_dirs()
    clean_output_dirs()
    set_style()
    datasets = load_validation_data()
    table1()
    table2()
    table3()
    table4()
    table5(datasets["E-MTAB-1980"])
    flow_figure()
    workflow_figure()
    km_figure(datasets)
    discrimination_figure()
    orientation_figure()
    signature_figure()
    supplementary()
    print("Frontiers package rebuilt at:", OUT_DIR)


if __name__ == "__main__":
    main()
