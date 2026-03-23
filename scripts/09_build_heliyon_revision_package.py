#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import gzip
import importlib.util
import sys
from collections import OrderedDict
from io import StringIO
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import chi2
from statsmodels.duration.hazard_regression import PHReg


ROOT = Path(__file__).resolve().parents[1]
BASE_SCRIPT = ROOT / "scripts" / "_figure_base.py"
OUT_DIR = ROOT / "output" / "heliyon_revision_package"
FIG_DIR = OUT_DIR / "figures"
TAB_DIR = OUT_DIR / "tables"
SUPP_DIR = OUT_DIR / "supplementary"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


base = load_module("submission_base_v5", BASE_SCRIPT)
base.OUT_DIR = OUT_DIR
base.FIG_DIR = FIG_DIR
base.TAB_DIR = TAB_DIR
base.SUPP_DIR = SUPP_DIR


COLORS = {
    "ink": "#13212D",
    "muted": "#647381",
    "grid": "#DBE2E9",
    "navy": "#1F5C81",
    "coral": "#C05A4B",
    "gold": "#B88B25",
    "teal": "#2E8577",
    "mauve": "#7A6E96",
    "panel": "#F7F9FB",
    "line": "#CCD5DD",
    "soft_blue": "#EAF1F6",
    "soft_red": "#F7ECE9",
}


def ensure_dirs():
    for path in (OUT_DIR, FIG_DIR, TAB_DIR, SUPP_DIR):
        path.mkdir(parents=True, exist_ok=True)


def set_style():
    plt.rcParams.update({
        "figure.dpi": 180,
        "savefig.dpi": 500,
        "savefig.bbox": "tight",
        "font.family": "DejaVu Sans",
        "font.size": 10.4,
        "axes.titlesize": 14.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 11.0,
        "axes.edgecolor": COLORS["ink"],
        "axes.linewidth": 0.85,
        "axes.facecolor": "white",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": COLORS["ink"],
        "ytick.color": COLORS["ink"],
        "text.color": COLORS["ink"],
        "axes.labelcolor": COLORS["ink"],
        "grid.color": COLORS["grid"],
        "grid.linewidth": 0.65,
        "grid.alpha": 0.9,
        "svg.fonttype": "none",
    })


def soften(ax):
    ax.tick_params(length=3.8, width=0.75)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.8)
        ax.spines[side].set_color(COLORS["ink"])


def nice_ticks(max_time):
    if max_time <= 18:
        step = 3
    elif max_time <= 36:
        step = 6
    elif max_time <= 72:
        step = 12
    else:
        step = 24
    ticks = np.arange(0, max_time + step, step)
    if len(ticks) < 5:
        ticks = np.linspace(0, max_time, 5)
    return ticks


def read_gse29609_clinicopath():
    with gzip.open(str(base.DATA_DIR / "GEO" / "GSE29609_series_matrix.txt.gz"), "rt", encoding="utf-8", errors="ignore") as handle:
        series_lines = handle.read().splitlines()
    accession_line = [x for x in series_lines if x.startswith("!Sample_geo_accession")][0]
    sample_ids = [x.replace('"', "") for x in accession_line.split("\t")[1:]]
    char_lines = [x for x in series_lines if x.startswith("!Sample_characteristics_ch1")]
    char_dt = pd.read_csv(StringIO("\n".join(char_lines)), sep="\t", header=None).iloc[:, 1:]
    char_dt.columns = sample_ids

    def get_field(values, prefix):
        prefix = prefix.lower()
        for item in values:
            if isinstance(item, str):
                val = item.strip('"')
                if val.lower().startswith(prefix):
                    return val.split(":", 1)[1].strip()
        return np.nan

    rows = []
    for sid in sample_ids:
        values = char_dt[sid].tolist()
        rows.append({
            "sample": sid,
            "age": pd.to_numeric(get_field(values, "age at diagnosis (y)"), errors="coerce"),
            "t_stage": pd.to_numeric(get_field(values, "t (tnm stage)"), errors="coerce"),
            "n_stage": pd.to_numeric(get_field(values, "n (tnm stage)"), errors="coerce"),
            "m_stage": pd.to_numeric(get_field(values, "m (tnm stage)"), errors="coerce"),
            "fuhrman_grade": pd.to_numeric(get_field(values, "fuhrman grade"), errors="coerce"),
            "tumor_necrosis": pd.to_numeric(get_field(values, "tumor necrosis"), errors="coerce"),
        })
    return pd.DataFrame(rows)


def load_tcga_signature_scores(signature):
    clin = base.load_tcga_baseline().copy()
    clin = clin[clin["sample_type"] == "Primary Tumor"].copy()
    clin = clin.drop_duplicates("bcr_patient_barcode")

    coef = signature.set_index("gene")["coef"]
    expr_path = base.DATA_DIR / "TCGA_Xena" / "TCGA.KIRC.sampleMap" / "HiSeqV2.gz"
    with gzip.open(str(expr_path), "rt", encoding="utf-8", errors="ignore") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        sample_ids = header[1:]
        gene_rows = {}
        wanted = set(coef.index)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            gene = parts[0].upper()
            if gene in wanted:
                gene_rows[gene] = [float(x) if x not in ("", "NA", "NaN") else np.nan for x in parts[1:]]

    expr = pd.DataFrame(gene_rows, index=sample_ids)
    for gene in expr.columns:
        expr[gene] = expr[gene].fillna(expr[gene].median())
    expr["signature_risk"] = expr[list(coef.index)].dot(coef)
    merged = clin.merge(expr[["signature_risk"]].reset_index().rename(columns={"index": "sampleID"}), on="sampleID", how="inner")
    merged["age"] = pd.to_numeric(merged["age_at_initial_pathologic_diagnosis"], errors="coerce")
    merged["grade_num"] = pd.to_numeric(merged["neoplasm_histologic_grade"].astype(str).str.extract(r"G(\d)", expand=False), errors="coerce")
    merged["stage_group"] = merged["pathologic_stage"].astype(str).str.strip().str.lower().map({
        "stage i": "I",
        "stage ii": "II",
        "stage iii": "III",
        "stage iv": "IV",
    })
    merged["stage_ord"] = merged["stage_group"].map({"I": 1, "II": 2, "III": 3, "IV": 4})
    return merged


def adapted_ssign_score(df):
    def t_points(x):
        if not np.isfinite(x):
            return np.nan
        return {1: 0, 2: 3, 3: 4, 4: 5}.get(int(x), np.nan)

    def n_points(x):
        if not np.isfinite(x):
            return np.nan
        return 0 if int(x) == 0 else 2

    def m_points(x):
        if not np.isfinite(x):
            return np.nan
        return 4 if int(x) == 1 else 0

    def grade_points(x):
        if not np.isfinite(x):
            return np.nan
        if int(x) in (1, 2):
            return 0
        if int(x) == 3:
            return 1
        if int(x) == 4:
            return 3
        return np.nan

    return (
        df["t_stage"].map(t_points)
        + df["n_stage"].map(n_points)
        + df["m_stage"].map(m_points)
        + df["fuhrman_grade"].map(grade_points)
        + 2 * df["tumor_necrosis"]
    )


def fit_cox_model(df, time_col, event_col, covars):
    sub = df[[time_col, event_col] + covars].dropna().copy()
    exog = sub[covars].astype(float).values
    model = PHReg(sub[time_col].astype(float).values, exog, status=sub[event_col].astype(int).values)
    result = model.fit(disp=0)
    lp = np.dot(exog, result.params)
    return {
        "n": int(sub.shape[0]),
        "events": int(sub[event_col].sum()),
        "llf": float(result.llf),
        "cindex": base.concordance_index(sub[time_col], sub[event_col], lp),
        "params": pd.Series(result.params, index=covars),
        "bse": pd.Series(result.bse, index=covars),
        "pvalues": pd.Series(result.pvalues, index=covars),
    }


def hr_ci(beta, se):
    return np.exp(beta), np.exp(beta - 1.96 * se), np.exp(beta + 1.96 * se)


def get_orientation_metrics(stat_bundle, key, seed_offset):
    stat = stat_bundle[key]
    lo, hi = base.bootstrap_cindex(stat["time"], stat["event"], stat["risk"], seed=base.SEED + seed_offset)
    return OrderedDict([
        ("p", stat["pval"]),
        ("cindex", stat["cindex"]),
        ("low", lo),
        ("high", hi),
    ])


def build_table_2_revised(gse_stats, cptac_stats):
    gse_orig = get_orientation_metrics(gse_stats["orientation_all"], "original", 101)
    gse_flip = OrderedDict([
        ("p", gse_stats["pval"]),
        ("cindex", gse_stats["cindex"]),
        ("low", gse_stats["cindex_low"]),
        ("high", gse_stats["cindex_high"]),
    ])
    cptac_orig = OrderedDict([
        ("p", cptac_stats["pval"]),
        ("cindex", cptac_stats["cindex"]),
        ("low", cptac_stats["cindex_low"]),
        ("high", cptac_stats["cindex_high"]),
    ])
    cptac_flip = get_orientation_metrics(cptac_stats["orientation_all"], "flipped", 104)

    table = pd.DataFrame([
        OrderedDict([
            ("Cohort", "GSE29609"),
            ("Endpoint", "death from cancer"),
            ("N", int(gse_stats["n"])),
            ("Events", int(gse_stats["n_events"])),
            ("Original score p", base.format_p(gse_orig["p"])),
            ("Original score C-index", base.format_float(gse_orig["cindex"])),
            ("Original score 95% CI", base.format_ci(gse_orig["low"], gse_orig["high"])),
            ("Sign-reversed score p", base.format_p(gse_flip["p"])),
            ("Sign-reversed score C-index", base.format_float(gse_flip["cindex"])),
            ("Sign-reversed 95% CI", base.format_ci(gse_flip["low"], gse_flip["high"])),
            ("Interpretation", "Directionally unstable; sign-reversed result treated as sensitivity analysis"),
        ]),
        OrderedDict([
            ("Cohort", "CPTAC"),
            ("Endpoint", "overall survival"),
            ("N", int(cptac_stats["n"])),
            ("Events", int(cptac_stats["n_events"])),
            ("Original score p", base.format_p(cptac_orig["p"])),
            ("Original score C-index", base.format_float(cptac_orig["cindex"])),
            ("Original score 95% CI", base.format_ci(cptac_orig["low"], cptac_orig["high"])),
            ("Sign-reversed score p", base.format_p(cptac_flip["p"])),
            ("Sign-reversed score C-index", base.format_float(cptac_flip["cindex"])),
            ("Sign-reversed 95% CI", base.format_ci(cptac_flip["low"], cptac_flip["high"])),
            ("Interpretation", "Validated in the original score direction"),
        ]),
    ])
    table.to_csv(str(TAB_DIR / "Table2_external_validation_performance.tsv"), sep="\t", index=False)
    base.write_markdown_table(table, TAB_DIR / "Table2_external_validation_performance.md")
    return table


def build_table_4_revised(gse_data, gse_stats, cptac_data, cptac_stats, cutoff_df):
    rows = []
    for endpoint in ("death", "death_cancer"):
        _, all_stats = base.choose_orientation(gse_data["time"], gse_data[endpoint], gse_data["risk"])
        for label, stat in all_stats.items():
            rows.append(OrderedDict([
                ("Analysis", "GSE29609 endpoint/orientation"),
                ("Cohort", "GSE29609"),
                ("Endpoint", endpoint),
                ("Orientation", label),
                ("N", stat["n"]),
                ("Log-rank p", base.format_p(stat["pval"])),
                ("C-index", base.format_float(stat["cindex"])),
            ]))
    for label, stat in cptac_stats["orientation_all"].items():
        rows.append(OrderedDict([
            ("Analysis", "CPTAC orientation"),
            ("Cohort", "CPTAC"),
            ("Endpoint", "OS"),
            ("Orientation", label),
            ("N", stat["n"]),
            ("Log-rank p", base.format_p(stat["pval"])),
            ("C-index", base.format_float(stat["cindex"])),
        ]))
    for cohort in ("GSE29609", "CPTAC"):
        sub = cutoff_df[cutoff_df["cohort"] == cohort]
        best_row = sub.sort_values("minus_log10_p", ascending=False).iloc[0]
        rows.append(OrderedDict([
            ("Analysis", "Cutoff sensitivity peak"),
            ("Cohort", cohort),
            ("Endpoint", "reported orientation"),
            ("Orientation", "best quantile={}".format(best_row["quantile"])),
            ("N", "NA"),
            ("Log-rank p", base.format_p(best_row["logrank_p"])),
            ("C-index", "NA"),
        ]))
    table = pd.DataFrame(rows)
    table.to_csv(str(TAB_DIR / "Table4_sensitivity_analyses.tsv"), sep="\t", index=False)
    base.write_markdown_table(table, TAB_DIR / "Table4_sensitivity_analyses.md")
    return table


def build_table_5_clinical_comparison(signature, gse_data):
    tcga = load_tcga_signature_scores(signature)
    tcga = tcga[tcga[["age", "stage_ord", "grade_num", "signature_risk"]].notna().all(axis=1)].copy()
    model_tcga_clin = fit_cox_model(tcga, "time", "status", ["age", "stage_ord", "grade_num"])
    model_tcga_sig = fit_cox_model(tcga, "time", "status", ["signature_risk"])
    model_tcga_comb = fit_cox_model(tcga, "time", "status", ["age", "stage_ord", "grade_num", "signature_risk"])
    tcga_lr_p = chi2.sf(2 * (model_tcga_comb["llf"] - model_tcga_clin["llf"]), 1)
    tcga_hr, tcga_lo, tcga_hi = hr_ci(model_tcga_comb["params"]["signature_risk"], model_tcga_comb["bse"]["signature_risk"])

    gse_clin = read_gse29609_clinicopath().merge(gse_data[["sample", "time", "death_cancer", "risk"]], on="sample", how="inner")
    gse_clin["adapted_ssign"] = adapted_ssign_score(gse_clin)
    gse_clin["risk_flipped"] = -gse_clin["risk"]
    gse_clin = gse_clin[gse_clin[["time", "death_cancer", "adapted_ssign", "risk_flipped"]].notna().all(axis=1)].copy()
    model_gse_clin = fit_cox_model(gse_clin, "time", "death_cancer", ["adapted_ssign"])
    model_gse_sig = fit_cox_model(gse_clin, "time", "death_cancer", ["risk_flipped"])
    model_gse_comb = fit_cox_model(gse_clin, "time", "death_cancer", ["adapted_ssign", "risk_flipped"])
    gse_lr_p = chi2.sf(2 * (model_gse_comb["llf"] - model_gse_clin["llf"]), 1)
    gse_hr, gse_lo, gse_hi = hr_ci(model_gse_comb["params"]["risk_flipped"], model_gse_comb["bse"]["risk_flipped"])

    table = pd.DataFrame([
        OrderedDict([
            ("Cohort", "TCGA-KIRC"),
            ("Clinical comparator", "Age + pathologic stage + grade"),
            ("N", model_tcga_comb["n"]),
            ("Events", model_tcga_comb["events"]),
            ("Clinical model C-index", base.format_float(model_tcga_clin["cindex"])),
            ("Signature-alone C-index", base.format_float(model_tcga_sig["cindex"])),
            ("Combined model C-index", base.format_float(model_tcga_comb["cindex"])),
            ("Signature HR in combined model", "{:.3f} ({:.3f}-{:.3f})".format(tcga_hr, tcga_lo, tcga_hi)),
            ("Signature p", base.format_p(model_tcga_comb["pvalues"]["signature_risk"])),
            ("Incremental LR p", base.format_p(tcga_lr_p)),
        ]),
        OrderedDict([
            ("Cohort", "GSE29609"),
            ("Clinical comparator", "Adapted SSIGN-like score (T/N/M, grade, necrosis; size unavailable)"),
            ("N", model_gse_comb["n"]),
            ("Events", model_gse_comb["events"]),
            ("Clinical model C-index", base.format_float(model_gse_clin["cindex"])),
            ("Signature-alone C-index", base.format_float(model_gse_sig["cindex"])),
            ("Combined model C-index", base.format_float(model_gse_comb["cindex"])),
            ("Signature HR in combined model", "{:.3f} ({:.3f}-{:.3f})".format(gse_hr, gse_lo, gse_hi)),
            ("Signature p", base.format_p(model_gse_comb["pvalues"]["risk_flipped"])),
            ("Incremental LR p", base.format_p(gse_lr_p)),
        ]),
    ])
    table.to_csv(str(TAB_DIR / "Table5_clinical_model_comparison.tsv"), sep="\t", index=False)
    base.write_markdown_table(table, TAB_DIR / "Table5_clinical_model_comparison.md")

    detail = pd.DataFrame([
        ("TCGA-KIRC", "Age + stage + grade", model_tcga_clin["n"], model_tcga_clin["events"], base.format_float(model_tcga_clin["cindex"]), "Baseline clinicopathologic model"),
        ("TCGA-KIRC", "Signature only", model_tcga_sig["n"], model_tcga_sig["events"], base.format_float(model_tcga_sig["cindex"]), "Submission-package score"),
        ("TCGA-KIRC", "Combined", model_tcga_comb["n"], model_tcga_comb["events"], base.format_float(model_tcga_comb["cindex"]), "Signature added independently to clinicopathologic model"),
        ("GSE29609", "Adapted SSIGN-like only", model_gse_clin["n"], model_gse_clin["events"], base.format_float(model_gse_clin["cindex"]), "Exploratory external comparison"),
        ("GSE29609", "Sign-reversed signature only", model_gse_sig["n"], model_gse_sig["events"], base.format_float(model_gse_sig["cindex"]), "Sensitivity orientation"),
        ("GSE29609", "Combined", model_gse_comb["n"], model_gse_comb["events"], base.format_float(model_gse_comb["cindex"]), "Exploratory external combination"),
    ], columns=["Cohort", "Model", "N", "Events", "C-index", "Notes"])
    detail.to_csv(str(SUPP_DIR / "clinical_model_details.tsv"), sep="\t", index=False)
    return table


def plot_workflow(signature, tcga_primary, gse_stats, cptac_stats):
    fig, ax = plt.subplots(figsize=(9.2, 14.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    # ── Title block ──────────────────────────────────────────────────
    ax.text(0.50, 0.975, "Figure 1.  Cross-platform prognostic modeling workflow",
            fontsize=17.5, fontweight="bold", va="top", ha="center",
            color=COLORS["ink"])
    ax.text(0.50, 0.946, "Vertical timeline from cohort assembly to transportability assessment",
            fontsize=10.0, color=COLORS["muted"], va="top", ha="center",
            fontstyle="italic")
    # thin decorative rule below subtitle
    ax.add_line(Line2D([0.18, 0.82], [0.935, 0.935],
                       color=COLORS["line"], lw=0.9, zorder=1))

    # ── Main panel background ────────────────────────────────────────
    ax.add_patch(mpatches.FancyBboxPatch(
        (0.075, 0.065), 0.85, 0.845,
        boxstyle="round,pad=0.015,rounding_size=0.028",
        linewidth=0, facecolor=COLORS["panel"]))

    # ── Timeline spine (gradient-like segments) ──────────────────────
    tl_x = 0.185
    seg_colors = [COLORS["navy"], COLORS["gold"], COLORS["teal"], COLORS["coral"]]
    y_tops = [0.855, 0.685, 0.505, 0.325]
    y_bots = [0.695, 0.515, 0.335, 0.155]
    for sc, yt, yb in zip(seg_colors, y_tops, y_bots):
        ax.add_line(Line2D([tl_x, tl_x], [yt, yb],
                           color=sc, lw=2.4, alpha=0.35, zorder=1))
    # thin overlay center-line for crispness
    ax.add_line(Line2D([tl_x, tl_x], [0.155, 0.855],
                       color=COLORS["line"], lw=0.7, zorder=2))

    # ── Steps definition ─────────────────────────────────────────────
    steps = [
        ("01", "COHORTS", "Discovery & validation sets",
         ["TCGA-KIRC primary tumors",
          "GSE29609 and CPTAC external cohorts"],
         COLORS["navy"], "#EDF3F8"),
        ("02", "HARMONIZATION", "Signal processing across platforms",
         ["Probe-level GEO signals collapsed to genes",
          "Within-cohort missing values repaired"],
         COLORS["gold"], "#F8F4EA"),
        ("03", "TRANSFER", "Signature projection & adjudication",
         ["Final 33-gene score projected into each cohort",
          "Original direction reported; sign reversal examined as sensitivity"],
         COLORS["teal"], "#EBF5F2"),
        ("04", "BENCHMARKING", "External performance benchmark",
         ["CPTAC validated in original direction; GSE29609 required sign-reversal sensitivity",
          "GSE29609 original C=0.262, reversed C=0.738 | CPTAC original C=0.722"],
         COLORS["coral"], "#F8EDEB"),
    ]

    y_positions = [0.765, 0.575, 0.385, 0.195]
    card_w = 0.58
    card_h = 0.145
    card_x = 0.265
    circle_size_outer = 1200
    circle_size_mid = 760
    circle_size_inner = 430

    for idx, ((num, kicker, title, lines, accent, fill), y) in enumerate(
            zip(steps, y_positions)):

        # ── Shadow layer behind card ─────────────────────────────────
        ax.add_patch(mpatches.FancyBboxPatch(
            (card_x + 0.004, y - 0.004), card_w, card_h,
            boxstyle="round,pad=0.012,rounding_size=0.022",
            linewidth=0, facecolor="#D7DCE1", alpha=0.38, zorder=2))

        # ── Card body ────────────────────────────────────────────────
        ax.add_patch(mpatches.FancyBboxPatch(
            (card_x, y), card_w, card_h,
            boxstyle="round,pad=0.012,rounding_size=0.022",
            linewidth=0.6, edgecolor=COLORS["line"], facecolor=fill,
            zorder=3))

        # ── Left accent bar ──────────────────────────────────────────
        bar_pad = 0.006
        ax.add_patch(mpatches.FancyBboxPatch(
            (card_x + bar_pad, y + 0.012), 0.008, card_h - 0.024,
            boxstyle="round,pad=0.002,rounding_size=0.004",
            linewidth=0, facecolor=accent, zorder=4))

        # ── Circle node on timeline ──────────────────────────────────
        # outer glow
        cy = y + card_h / 2
        ax.scatter([tl_x], [cy], s=circle_size_outer, color=accent, alpha=0.12, zorder=3)
        ax.scatter([tl_x], [cy], s=circle_size_mid, facecolors="white", edgecolors=accent, linewidths=2.2, zorder=4)
        ax.scatter([tl_x], [cy], s=circle_size_inner, color=accent, zorder=5)
        ax.text(tl_x, cy, num,
                ha="center", va="center", color="white",
                fontsize=9.0, fontweight="bold", zorder=6)

        # ── Connector dash from circle to card ───────────────────────
        ax.add_line(Line2D(
            [tl_x + 0.035, card_x - 0.008],
            [cy, cy],
            color=accent, lw=1.0, linestyle=(0, (4, 3)),
            alpha=0.55, zorder=3))

        # ── Text inside card ─────────────────────────────────────────
        tx = card_x + 0.040
        ax.text(tx, y + card_h - 0.018, kicker,
                fontsize=8.5, fontweight="bold", color=accent,
                va="top", zorder=5, fontstyle="normal",
                fontfamily="DejaVu Sans")
        ax.text(tx, y + card_h - 0.038, title,
                fontsize=12.5, fontweight="bold", va="top", zorder=5)
        detail_block = "{}\n{}".format(
            base.wrap_label(lines[0], 40),
            base.wrap_label(lines[1], 40),
        )
        ax.text(tx, y + 0.060, detail_block,
                fontsize=8.35, color=COLORS["muted"], va="top", zorder=5, linespacing=1.32)

        # ── Down-arrow between steps ─────────────────────────────────
        if idx < len(steps) - 1:
            arr_y = y - 0.012
            ax.annotate("", xy=(tl_x, arr_y - 0.028),
                        xytext=(tl_x, arr_y + 0.005),
                        arrowprops=dict(
                            arrowstyle="-|>",
                            color=seg_colors[idx],
                            lw=1.6, mutation_scale=11),
                        zorder=6)

    # ── Bottom gene info pill ────────────────────────────────────────
    top_genes = signature.sort_values("abs_coef", ascending=False).head(5)["gene"].tolist()
    pill_y = 0.082
    ax.add_patch(mpatches.FancyBboxPatch(
        (0.14, pill_y), 0.72, 0.048,
        boxstyle="round,pad=0.010,rounding_size=0.020",
        linewidth=1.0, edgecolor=COLORS["mauve"], facecolor="#F4F1F8",
        zorder=3))
    ax.text(0.17, pill_y + 0.024,
            "\u25C6  Top |coef| genes:  " + ",  ".join(top_genes),
            fontsize=9.2, color=COLORS["mauve"], va="center",
            fontweight="bold", zorder=4)

    base.save_figure(fig, "Figure1_study_workflow")


def km_panel(fig, spec, stat, title, tag):
    inner = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=spec, height_ratios=[1.15, 4.35, 1.15], hspace=0.08)
    ax_head = fig.add_subplot(inner[0])
    ax = fig.add_subplot(inner[1])
    ax_tab = fig.add_subplot(inner[2])

    time = stat["time"]
    event = stat["event"]
    group = stat["group_used"]
    max_time = float(np.nanmax(time))
    ticks = nice_ticks(max_time)

    ax_head.axis("off")
    ax_head.text(0.00, 0.90, "{}  {}".format(tag, base.wrap_label(title, 34)), fontsize=12.7, fontweight="bold", va="top")
    ax_head.text(0.00, 0.26, "p={}   C-index={} ({})   n={}   events={}".format(base.format_p(stat["pval"]), base.format_float(stat["cindex"]), base.format_ci(stat["cindex_low"], stat["cindex_high"]), stat["n"], stat["n_events"]), fontsize=9.3, color=COLORS["muted"], va="top")

    for grp, line_color, fill_color in ((0, COLORS["navy"], COLORS["soft_blue"]), (1, COLORS["coral"], COLORS["soft_red"])):
        mask = group == grp
        curve = base.km_curve(time[mask], event[mask])
        xs, ys = base.stepify(curve["times"], curve["surv"])
        xl, yl = base.stepify(curve["times"], curve["lower"])
        xu, yu = base.stepify(curve["times"], curve["upper"])
        ax.fill_between(xl, yl, yu, step="post", color=fill_color, linewidth=0)
        ax.step(xs, ys, where="post", color=line_color, linewidth=2.25)

    ax.set_xlim(0, max_time * 1.03)
    ax.set_ylim(0, 1.02)
    ax.set_xticks(ticks)
    ax.set_xticklabels([])
    ax.set_ylabel("Survival probability")
    ax.grid(axis="y")
    soften(ax)

    ax_tab.axis("off")
    ax_tab.set_xlim(0, 1)
    ax_tab.set_ylim(0, 1)
    ax_tab.add_patch(mpatches.FancyBboxPatch((0.15, 0.06), 0.81, 0.76, boxstyle="round,pad=0.01,rounding_size=0.02", linewidth=0, facecolor=COLORS["panel"]))
    ax_tab.add_line(Line2D([0.22, 0.27], [0.88, 0.88], transform=ax_tab.transAxes, color=COLORS["navy"], lw=2.3))
    ax_tab.text(0.285, 0.88, "Low risk", transform=ax_tab.transAxes, fontsize=9.2, color=COLORS["muted"], va="center")
    ax_tab.add_line(Line2D([0.48, 0.53], [0.88, 0.88], transform=ax_tab.transAxes, color=COLORS["coral"], lw=2.3))
    ax_tab.text(0.545, 0.88, "High risk", transform=ax_tab.transAxes, fontsize=9.2, color=COLORS["muted"], va="center")
    ax_tab.text(0.02, 0.62, "No. at risk", fontsize=9.6, fontweight="bold")
    col_x = np.linspace(0.23, 0.94, len(ticks))
    for x, t in zip(col_x, ticks):
        ax_tab.text(x, 0.70, "{:.0f}".format(t), ha="center", fontsize=8.9, color=COLORS["muted"])
    low_counts = base.risk_counts(time[group == 0], ticks)
    high_counts = base.risk_counts(time[group == 1], ticks)
    ax_tab.text(0.02, 0.40, "Low risk", fontsize=9.0, color=COLORS["navy"], fontweight="bold")
    ax_tab.text(0.02, 0.14, "High risk", fontsize=9.0, color=COLORS["coral"], fontweight="bold")
    for x, n in zip(col_x, low_counts):
        ax_tab.text(x, 0.40, str(n), ha="center", fontsize=8.9)
    for x, n in zip(col_x, high_counts):
        ax_tab.text(x, 0.14, str(n), ha="center", fontsize=8.9)
    ax_tab.text(0.58, -0.02, "Time (months)", ha="center", fontsize=9.8)


def plot_external_km(gse_stats, cptac_stats):
    fig = plt.figure(figsize=(16.2, 8.6))
    outer = gridspec.GridSpec(1, 2, figure=fig, wspace=0.20)
    fig.suptitle("Figure 2. External survival stratification plots", x=0.06, y=0.99, ha="left", fontsize=18, fontweight="bold")
    km_panel(fig, outer[0], gse_stats, "GSE29609 disease-specific survival (sign-reversed sensitivity)", "A")
    km_panel(fig, outer[1], cptac_stats, "CPTAC overall survival", "B")
    base.save_figure(fig, "Figure2_external_kaplan_meier")


def plot_discrimination_summary(gse_stats, cptac_stats):
    df = pd.DataFrame([
        OrderedDict([("cohort", "GSE29609"), ("cindex", gse_stats["cindex"]), ("low", gse_stats["cindex_low"]), ("high", gse_stats["cindex_high"]), ("p", base.format_p(gse_stats["pval"]))]),
        OrderedDict([("cohort", "CPTAC"), ("cindex", cptac_stats["cindex"]), ("low", cptac_stats["cindex_low"]), ("high", cptac_stats["cindex_high"]), ("p", base.format_p(cptac_stats["pval"]))]),
    ])
    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    y = np.arange(df.shape[0])[::-1]
    colors = [COLORS["coral"], COLORS["navy"]]
    ax.axvspan(0.60, 0.90, color="#E6F1EC", zorder=0)
    ax.axvline(0.60, color="#2D7A58", linestyle="--", linewidth=1.0)
    for yi, (_, row), color in zip(y, df.iterrows(), colors):
        ax.hlines(yi, row["low"], row["high"], color=color, linewidth=1.8, zorder=2)
        ax.scatter(row["cindex"], yi, s=58, color=color, edgecolor="white", linewidth=0.9, zorder=3)
        ax.text(row["high"] + 0.010, yi, "p={}".format(row["p"]), va="center", fontsize=9.0, color=COLORS["muted"])
    ax.set_yticks(y)
    ax.set_yticklabels(df["cohort"])
    ax.set_xlim(0.52, max(df["high"]) + 0.10)
    ax.set_xlabel("Harrell's C-index")
    ax.set_title("Figure 3. Sensitivity-adjusted discrimination summary", loc="left")
    ax.grid(axis="x")
    soften(ax)
    base.save_figure(fig, "Figure3_discrimination_summary")


def plot_signature_characteristics(signature):
    df = signature.sort_values("coef").copy().reset_index(drop=True)
    cat_colors = {"immune": COLORS["coral"], "metabolic": COLORS["teal"], "mixed": COLORS["mauve"]}

    fig = plt.figure(figsize=(13.8, 13.2))
    outer = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1.2, 9.0], hspace=0.04)
    ax_head = fig.add_subplot(outer[0])
    inner = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1], width_ratios=[5.3, 1.8], wspace=0.05)
    ax = fig.add_subplot(inner[0])
    ax_cat = fig.add_subplot(inner[1], sharey=ax)

    ax_head.axis("off")
    ax_head.text(0.00, 0.88, "Figure 5. Coefficient architecture of the final signature", fontsize=15.8, fontweight="bold", va="top")
    ax_head.text(0.00, 0.28, "Coefficient direction is shown at left; curated immune/metabolic labels are isolated in a separate right-hand annotation column.", fontsize=9.5, color=COLORS["muted"], va="top")
    ax_head.text(0.89, 0.28, "Category", fontsize=10.4, fontweight="bold", va="top", ha="center")

    y = np.arange(df.shape[0])
    for i in range(len(y)):
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color=COLORS["panel"], zorder=0)
            ax_cat.axhspan(i - 0.5, i + 0.5, color=COLORS["panel"], zorder=0)
    ax.axvspan(df["coef"].min() - 0.22, 0, color=COLORS["soft_red"], zorder=0)
    ax.axvspan(0, df["coef"].max() + 0.22, color=COLORS["soft_blue"], zorder=0)

    for idx, row in df.iterrows():
        color = cat_colors[row["category"]]
        ax.hlines(idx, 0, row["coef"], color=color, linewidth=2.0)
        ax.scatter(row["coef"], idx, s=56 + row["abs_coef"] * 30, color=color, edgecolor="white", linewidth=0.85, zorder=3)

    ax.axvline(0, color=COLORS["ink"], lw=0.95)
    ax.set_yticks(y)
    ax.set_yticklabels(df["gene"], fontsize=9.0)
    ax.set_xlabel("Signature coefficient")
    ax.grid(axis="x")
    soften(ax)

    ax_cat.set_xlim(0, 1)
    ax_cat.set_xticks([])
    ax_cat.tick_params(left=False, labelleft=False)
    for spine in ax_cat.spines.values():
        spine.set_visible(False)
    ax_cat.set_ylim(-0.5, len(df) - 0.5)
    for idx, row in df.iterrows():
        ax_cat.text(0.14, idx, row["category"].capitalize(), va="center", fontsize=8.9, color=cat_colors[row["category"]])

    base.save_figure(fig, "Figure5_signature_characteristics")


def robust_orientation_panel(fig, spec, stat, title, tag):
    inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=spec, height_ratios=[0.9, 4.1], hspace=0.05)
    ax_head = fig.add_subplot(inner[0])
    ax = fig.add_subplot(inner[1])
    both = stat["orientation_all"]
    df = pd.DataFrame([
        OrderedDict([("label", "Original"), ("cindex", both["original"]["cindex"]), ("p", base.format_p(both["original"]["pval"]))]),
        OrderedDict([("label", "Flipped"), ("cindex", both["flipped"]["cindex"]), ("p", base.format_p(both["flipped"]["pval"]))]),
    ])

    ax_head.axis("off")
    ax_head.text(0.00, 0.84, "{}  {}".format(tag, title), fontsize=13.3, fontweight="bold", va="top")
    ax_head.text(0.00, 0.24, "Sign reversal is shown as a transportability sensitivity analysis, not confirmatory validation.", fontsize=9.3, color=COLORS["muted"], va="top")

    y = np.array([1, 0])
    colors = [COLORS["muted"], COLORS["gold"]]
    ax.hlines(y, 0.20, df["cindex"], color=colors, linewidth=2.0)
    ax.scatter(df["cindex"], y, s=72, color=colors, edgecolor="white", linewidth=0.9, zorder=3)
    ax.axvline(0.60, linestyle="--", color="#2D7A58", linewidth=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels(df["label"])
    ax.set_xlim(0.20, 0.88)
    ax.set_xlabel("C-index")
    ax.grid(axis="x")
    soften(ax)
    ax.text(0.865, 1, "p={}".format(df.iloc[0]["p"]), ha="right", va="center", fontsize=9.0, color=COLORS["muted"])
    ax.text(0.865, 0, "p={}".format(df.iloc[1]["p"]), ha="right", va="center", fontsize=9.0, color=COLORS["muted"])


def plot_robustness(gse_stats, cptac_stats, cutoff_df):
    fig = plt.figure(figsize=(15.0, 10.5))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.34, wspace=0.22)
    fig.suptitle("Figure 4. Robustness diagnostics", x=0.06, y=0.99, ha="left", fontsize=18, fontweight="bold")

    robust_orientation_panel(fig, gs[0, 0], gse_stats, "GSE29609 orientation handling", "A")
    robust_orientation_panel(fig, gs[0, 1], cptac_stats, "CPTAC orientation handling", "B")

    inner_c = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1, 0], height_ratios=[0.9, 4.1], hspace=0.05)
    ax_c_head = fig.add_subplot(inner_c[0])
    ax_c = fig.add_subplot(inner_c[1])
    ax_c_head.axis("off")
    ax_c_head.text(0.00, 0.84, "C  Sensitivity to alternative cutoffs", fontsize=13.3, fontweight="bold", va="top")
    ax_c_head.text(0.00, 0.24, "CPTAC remains stable, whereas GSE29609 varies more because of the small sample size.", fontsize=9.3, color=COLORS["muted"], va="top")
    for cohort, color in (("GSE29609", COLORS["coral"]), ("CPTAC", COLORS["navy"])):
        sub = cutoff_df[cutoff_df["cohort"] == cohort]
        ax_c.plot(sub["quantile"], sub["minus_log10_p"], marker="o", markersize=4.3, linewidth=2.0, color=color, label=cohort)
    ax_c.axhline(-np.log10(0.05), linestyle="--", color="#2D7A58", linewidth=1.0)
    ax_c.set_xlabel("Risk-score cutoff quantile")
    ax_c.set_ylabel("-log10(log-rank p)")
    ax_c.grid(True)
    soften(ax_c)
    ax_c.legend(frameon=False, loc="upper left")

    inner_d = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1, 1], height_ratios=[0.9, 4.1], hspace=0.05)
    ax_d_head = fig.add_subplot(inner_d[0])
    ax_d = fig.add_subplot(inner_d[1])
    ax_d_head.axis("off")
    ax_d_head.text(0.00, 0.84, "D  High-risk group balance", fontsize=13.3, fontweight="bold", va="top")
    ax_d_head.text(0.00, 0.24, "Risk-group proportion stays stable as the quantile cutoff is varied.", fontsize=9.3, color=COLORS["muted"], va="top")
    for cohort, color in (("GSE29609", COLORS["coral"]), ("CPTAC", COLORS["navy"])):
        sub = cutoff_df[cutoff_df["cohort"] == cohort]
        ax_d.plot(sub["quantile"], sub["high_risk_fraction"], marker="o", markersize=4.3, linewidth=2.0, color=color, label=cohort)
    ax_d.axhspan(0.35, 0.65, color="#E6F1EC", zorder=0)
    ax_d.set_xlabel("Risk-score cutoff quantile")
    ax_d.set_ylabel("High-risk fraction")
    ax_d.set_ylim(0, 1)
    ax_d.grid(True)
    soften(ax_d)
    ax_d.legend(frameon=False, loc="upper left")

    base.save_figure(fig, "Figure4_robustness_diagnostics")


def build_support_files(signature, gse_stats, cptac_stats):
    s1 = ["# S1. Optimization log excerpts", ""]
    log_files = [
        base.RESULTS_DIR / "best_signature.txt",
        base.RESULTS_DIR / "opt_final_best.txt",
        base.RESULTS_DIR / "run_summary.txt",
    ]
    existing_logs = [path for path in log_files if path.exists()]
    if existing_logs:
        s1.extend(["Original optimization logs are preserved below.", ""])
        for path in existing_logs:
            s1.extend([
                "## results/{}".format(path.name),
                "```",
                base.read_text_excerpt(path),
                "```",
                "",
            ])
    else:
        s1.append("No optional optimization log excerpts were found in `results/`; the package build remains reproducible without them.")
    (SUPP_DIR / "S1_optimization_log_excerpts.md").write_text("\n".join(s1) + "\n", encoding="utf-8")

    s2 = [
        "# S2. Reproducibility script list and runtime settings",
        "",
        "## Preserved code",
        "- `scripts/_figure_base.py`",
        "- `scripts/09_build_heliyon_revision_package.py`",
        "",
        "## Output directory",
        "- `output/heliyon_revision_package/`",
        "",
        "## Runtime settings",
        "- Random seed: {}".format(base.SEED),
        "- Bootstrap resamples: {}".format(base.BOOTSTRAP_N),
        "- Reporting stance: original orientation retained as primary; sign reversal reported as sensitivity",
        "- New Table 5 adds clinicopathologic model comparison",
    ]
    (SUPP_DIR / "S2_reproducibility_script_list_and_runtime_settings.md").write_text("\n".join(s2) + "\n", encoding="utf-8")

    version_df = pd.DataFrame([
        ("python", sys.version.split()[0]),
        ("pandas", pd.__version__),
        ("numpy", np.__version__),
        ("matplotlib", matplotlib.__version__),
    ], columns=["Software / dependency", "Version"])
    version_df.to_csv(str(SUPP_DIR / "S3_software_dependency_versions.tsv"), sep="\t", index=False)
    base.write_markdown_table(version_df, SUPP_DIR / "S3_software_dependency_versions.md")

    summary_df = pd.DataFrame([
        OrderedDict([("Cohort", "GSE29609"), ("Endpoint", gse_stats["endpoint"]), ("Original C-index", base.format_float(gse_stats["orientation_all"]["original"]["cindex"])), ("Sign-reversed C-index", base.format_float(gse_stats["orientation_all"]["flipped"]["cindex"])), ("Interpretation", "Sensitivity only")]),
        OrderedDict([("Cohort", "CPTAC"), ("Endpoint", cptac_stats["endpoint"]), ("Original C-index", base.format_float(cptac_stats["orientation_all"]["original"]["cindex"])), ("Sign-reversed C-index", base.format_float(cptac_stats["orientation_all"]["flipped"]["cindex"])), ("Interpretation", "Original direction validated")]),
    ])
    summary_df.to_csv(str(SUPP_DIR / "summary_metrics.tsv"), sep="\t", index=False)

    readme = [
        "# Heliyon Revision Submission Package",
        "",
        "This package revises the manuscript for Heliyon and addresses reviewer concerns about sign inversion, small-sample uncertainty, and missing clinicopathologic comparison.",
        "",
        "## Code",
        "- Shared base generator: `scripts/_figure_base.py`",
        "- Heliyon revision builder: `scripts/09_build_heliyon_revision_package.py`",
        "",
        "## Key external results",
        "- GSE29609 original direction: p={}, C-index={}".format(base.format_p(gse_stats["orientation_all"]["original"]["pval"]), base.format_float(gse_stats["orientation_all"]["original"]["cindex"])),
        "- GSE29609 sign-reversed sensitivity: p={}, C-index={} ({})".format(base.format_p(gse_stats["orientation_all"]["flipped"]["pval"]), base.format_float(gse_stats["orientation_all"]["flipped"]["cindex"]), base.format_ci(gse_stats["cindex_low"], gse_stats["cindex_high"])),
        "- CPTAC original direction: p={}, C-index={} ({})".format(base.format_p(cptac_stats["orientation_all"]["original"]["pval"]), base.format_float(cptac_stats["orientation_all"]["original"]["cindex"]), base.format_ci(cptac_stats["cindex_low"], cptac_stats["cindex_high"])),
        "",
        "## New main-table addition",
        "- Table 5 compares the molecular signature with clinicopathologic models in TCGA and GSE29609.",
    ]
    (OUT_DIR / "README_submission_package.md").write_text("\n".join(readme) + "\n", encoding="utf-8")


def main():
    ensure_dirs()
    set_style()

    signature = base.load_signature()
    tcga_primary = base.load_tcga_baseline()
    gse_data = base.load_gse29609(signature)
    cptac_data = base.load_cptac(signature)

    gse_stats = base.analyze_external_cohort(gse_data, "time", "death_cancer", "GSE29609", "death from cancer")
    cptac_stats = base.analyze_external_cohort(cptac_data, "OS_MONTHS", "OS_STATUS", "CPTAC", "overall survival")

    cutoff_df = pd.concat([
        base.compute_cutoff_sensitivity(gse_stats["time"], gse_stats["event"], gse_stats["risk_used"], "GSE29609"),
        base.compute_cutoff_sensitivity(cptac_stats["time"], cptac_stats["event"], cptac_stats["risk_used"], "CPTAC"),
    ], ignore_index=True)
    cutoff_df.to_csv(str(SUPP_DIR / "cutoff_sensitivity.tsv"), sep="\t", index=False)

    base.build_table_1(tcga_primary, gse_stats, cptac_stats, gse_data, cptac_data)
    build_table_2_revised(gse_stats, cptac_stats)
    base.build_table_3(signature)
    build_table_4_revised(gse_data, gse_stats, cptac_data, cptac_stats, cutoff_df)
    build_table_5_clinical_comparison(signature, gse_data)

    plot_workflow(signature, tcga_primary, gse_stats, cptac_stats)
    plot_external_km(gse_stats, cptac_stats)
    plot_discrimination_summary(gse_stats, cptac_stats)
    plot_signature_characteristics(signature)
    plot_robustness(gse_stats, cptac_stats, cutoff_df)
    build_support_files(signature, gse_stats, cptac_stats)

    print("Heliyon revision package generated under:", str(OUT_DIR))


if __name__ == "__main__":
    main()
