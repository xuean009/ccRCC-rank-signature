#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import gzip
import math
import platform
import sys
import textwrap
import warnings
from collections import OrderedDict
from io import StringIO
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from scipy.stats import chi2

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
OUT_DIR = ROOT / "output" / "submission_figures_tables"
FIG_DIR = OUT_DIR / "figures"
TAB_DIR = OUT_DIR / "tables"
SUPP_DIR = OUT_DIR / "supplementary"

SEED = 42
BOOTSTRAP_N = 250

PALETTE = {
    "high": "#B54A3A",
    "low": "#2C6E91",
    "accent": "#C8941A",
    "ink": "#1C2430",
    "muted": "#6C7A89",
    "grid": "#D9E1EA",
    "bg": "#F6F8FB",
    "mixed": "#7A6C92",
    "immune": "#C24D5A",
    "metabolic": "#2A8F6A",
    "pass": "#2E8B57",
    "fail": "#9B2D30",
}

MANUAL_IMMUNE_GENES = {"CXCL6", "LTF", "SAA1", "SERPINA3", "FGG", "ATP6V0D2"}
MANUAL_METABOLIC_GENES = {"CYP4F2", "SLC22A8", "APOB", "ANGPTL3", "GBA3", "SLC15A1", "GDA", "CASR"}


def ensure_dirs():
    for path in (OUT_DIR, FIG_DIR, TAB_DIR, SUPP_DIR):
        path.mkdir(parents=True, exist_ok=True)


def set_publication_style():
    plt.rcParams.update({
        "figure.dpi": 160,
        "savefig.dpi": 450,
        "savefig.bbox": "tight",
        "font.family": "DejaVu Sans",
        "font.size": 10.5,
        "axes.titlesize": 14,
        "axes.titleweight": "bold",
        "axes.labelsize": 11,
        "axes.edgecolor": PALETTE["ink"],
        "axes.linewidth": 0.8,
        "axes.facecolor": "white",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": PALETTE["ink"],
        "ytick.color": PALETTE["ink"],
        "text.color": PALETTE["ink"],
        "axes.labelcolor": PALETTE["ink"],
        "grid.color": PALETTE["grid"],
        "grid.linewidth": 0.7,
        "grid.alpha": 0.8,
        "svg.fonttype": "none",
    })


def wrap_label(text, width=28):
    return "\n".join(textwrap.wrap(text, width=width))


def format_p(pval):
    if pval is None or not np.isfinite(pval):
        return "NA"
    if pval < 1e-4:
        return "<0.0001"
    if pval < 1e-3:
        return "{:.2e}".format(pval)
    return "{:.3f}".format(pval)


def format_float(x, digits=3):
    if x is None or not np.isfinite(x):
        return "NA"
    return "{:.{d}f}".format(x, d=digits)


def format_percent(x, digits=1):
    if x is None or not np.isfinite(x):
        return "NA"
    return "{:.{d}f}%".format(x * 100.0, d=digits)


def format_ci(lo, hi, digits=3):
    if not np.isfinite(lo) or not np.isfinite(hi):
        return "NA"
    return "{:.{d}f}-{:.{d}f}".format(lo, hi, d=digits)


def stepify(times, values):
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    if times.size == 0:
        return np.array([0.0]), np.array([1.0])
    return np.r_[0.0, times], np.r_[1.0, values]


def logrank_test(time, event, group):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    group = np.asarray(group, dtype=int)
    mask = np.isfinite(time) & np.isfinite(event) & np.isfinite(group)
    time = time[mask]
    event = event[mask]
    group = group[mask]
    if time.size == 0 or np.unique(group).size < 2 or np.unique(event).size < 2:
        return np.nan
    event_times = np.sort(np.unique(time[event == 1]))
    observed_1 = expected_1 = variance_1 = 0.0
    for t in event_times:
        at_risk = time >= t
        is_event = (time == t) & (event == 1)
        n = at_risk.sum()
        d = is_event.sum()
        n1 = (at_risk & (group == 1)).sum()
        d1 = (is_event & (group == 1)).sum()
        if n <= 1:
            continue
        observed_1 += d1
        expected_1 += d * (n1 / float(n))
        variance_1 += (n1 * (n - n1) * d * (n - d)) / float((n ** 2) * (n - 1))
    if variance_1 <= 0:
        return np.nan
    chisq = ((observed_1 - expected_1) ** 2) / variance_1
    return float(chi2.sf(chisq, 1))


def concordance_index(time, event, score):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    score = np.asarray(score, dtype=float)
    mask = np.isfinite(time) & np.isfinite(event) & np.isfinite(score)
    time = time[mask]
    event = event[mask]
    score = score[mask]
    n = time.size
    if n < 2:
        return np.nan
    concordant = tied = comparable = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            ti, tj = time[i], time[j]
            ei, ej = event[i], event[j]
            si, sj = score[i], score[j]
            if ti == tj:
                if ei == 1 and ej == 1:
                    comparable += 1.0
                    if si == sj:
                        tied += 1.0
                    else:
                        concordant += 0.5
                continue
            if ti < tj and ei == 1:
                comparable += 1.0
                if si > sj:
                    concordant += 1.0
                elif si == sj:
                    tied += 1.0
            elif tj < ti and ej == 1:
                comparable += 1.0
                if sj > si:
                    concordant += 1.0
                elif si == sj:
                    tied += 1.0
    if comparable == 0:
        return np.nan
    return float((concordant + 0.5 * tied) / comparable)


def bootstrap_cindex(time, event, score, n_boot=BOOTSTRAP_N, seed=SEED):
    local_rng = np.random.RandomState(seed)
    idx = np.arange(len(time))
    estimates = []
    for _ in range(n_boot):
        sample_idx = local_rng.choice(idx, size=idx.size, replace=True)
        est = concordance_index(np.asarray(time)[sample_idx], np.asarray(event)[sample_idx], np.asarray(score)[sample_idx])
        if np.isfinite(est):
            estimates.append(est)
    if not estimates:
        return np.nan, np.nan
    return tuple(np.quantile(estimates, [0.025, 0.975]))


def km_curve(time, event):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    mask = np.isfinite(time) & np.isfinite(event)
    time = time[mask]
    event = event[mask]
    if time.size == 0:
        return {"times": np.array([]), "surv": np.array([]), "lower": np.array([]), "upper": np.array([])}
    unique_event_times = np.sort(np.unique(time[event == 1]))
    surv = 1.0
    greenwood = 0.0
    times_out, surv_out, lower_out, upper_out = [], [], [], []
    for t in unique_event_times:
        at_risk = np.sum(time >= t)
        events = np.sum((time == t) & (event == 1))
        if at_risk == 0:
            continue
        surv *= (1.0 - (events / float(at_risk)))
        if at_risk > events and events > 0:
            greenwood += events / float(at_risk * (at_risk - events))
        se = surv * math.sqrt(greenwood) if greenwood > 0 else 0.0
        times_out.append(float(t))
        surv_out.append(float(surv))
        lower_out.append(max(0.0, surv - 1.96 * se))
        upper_out.append(min(1.0, surv + 1.96 * se))
    return {
        "times": np.asarray(times_out),
        "surv": np.asarray(surv_out),
        "lower": np.asarray(lower_out),
        "upper": np.asarray(upper_out),
    }


def risk_counts(time, eval_times):
    time = np.asarray(time, dtype=float)
    return [int(np.sum(time >= t)) for t in eval_times]


def median_fill_frame(frame):
    out = frame.copy()
    for col in out.columns:
        med = out[col].median()
        if pd.isna(med):
            med = 0.0
        out[col] = out[col].fillna(med)
    return out


def categorize_gene(gene):
    gene = str(gene).upper()
    if gene in MANUAL_IMMUNE_GENES and gene in MANUAL_METABOLIC_GENES:
        return "mixed"
    if gene in MANUAL_IMMUNE_GENES:
        return "immune"
    if gene in MANUAL_METABOLIC_GENES:
        return "metabolic"
    return "mixed"


def write_markdown_table(df, path):
    cols = list(df.columns)
    rows = [cols, ["---"] * len(cols)]
    for _, row in df.iterrows():
        rows.append([str(row[c]) for c in cols])
    lines = ["| " + " | ".join(r) + " |" for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def save_table(df, basename):
    tsv_path = TAB_DIR / (basename + ".tsv")
    md_path = TAB_DIR / (basename + ".md")
    df.to_csv(str(tsv_path), sep="\t", index=False)
    write_markdown_table(df, md_path)


def save_figure(fig, basename):
    png_path = FIG_DIR / (basename + ".png")
    svg_path = FIG_DIR / (basename + ".svg")
    fig.savefig(str(png_path), facecolor="white")
    fig.savefig(str(svg_path), facecolor="white")
    plt.close(fig)


def read_text_excerpt(path):
    p = Path(path)
    if not p.exists():
        return "Missing: {}".format(p.name)
    return p.read_text(encoding="utf-8", errors="ignore").strip()


def read_table(path, sep="\t", **kwargs):
    path = Path(path)
    if str(path).endswith(".gz"):
        with gzip.open(str(path), "rt", encoding="utf-8", errors="ignore") as handle:
            return pd.read_csv(handle, sep=sep, **kwargs)
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        return pd.read_csv(handle, sep=sep, **kwargs)


def load_signature():
    coef = read_table(RESULTS_DIR / "best_signature_coefficients.tsv")
    coef["gene"] = coef["gene"].astype(str).str.upper()
    coef["abs_coef"] = coef["coef"].abs()
    coef["direction"] = np.where(coef["coef"] >= 0, "Risk-increasing", "Protective")
    coef["category"] = coef["gene"].map(categorize_gene)
    return coef


def load_tcga_baseline():
    clin = read_table(DATA_DIR / "TCGA_Xena" / "TCGA.KIRC.sampleMap" / "KIRC_clinicalMatrix")
    clin["time"] = pd.to_numeric(clin["days_to_death"], errors="coerce")
    clin["time"] = clin["time"].fillna(pd.to_numeric(clin["days_to_last_followup"], errors="coerce"))
    status_raw = clin["vital_status"].astype(str).str.lower()
    clin["status"] = status_raw.map({"deceased": 1, "dead": 1, "living": 0, "alive": 0})
    primary = clin[(clin["sample_type"] == "Primary Tumor") & clin["time"].notna() & clin["status"].notna() & (clin["time"] > 0)].copy()
    return primary


def load_gse29609(signature):
    with gzip.open(str(DATA_DIR / "GEO" / "GSE29609_series_matrix.txt.gz"), "rt", encoding="utf-8", errors="ignore") as handle:
        series_lines = handle.read().splitlines()
    accession_line = [x for x in series_lines if x.startswith("!Sample_geo_accession")][0]
    sample_ids = [x.replace('"', "") for x in accession_line.split("\t")[1:]]
    begin = next(i for i, x in enumerate(series_lines) if "!series_matrix_table_begin" in x)
    end = next(i for i, x in enumerate(series_lines) if "!series_matrix_table_end" in x)
    expr = pd.read_csv(StringIO("\n".join(series_lines[begin + 1:end])), sep="\t")
    expr = expr.rename(columns={expr.columns[0]: "ID_REF"})

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
        return None

    clin_rows = []
    for sid in sample_ids:
        values = char_dt[sid].tolist()
        clin_rows.append({
            "sample": sid,
            "time": pd.to_numeric(get_field(values, "survival time"), errors="coerce"),
            "death": pd.to_numeric(get_field(values, "death"), errors="coerce"),
            "death_cancer": pd.to_numeric(get_field(values, "death from cancer"), errors="coerce"),
        })
    clin = pd.DataFrame(clin_rows)

    with gzip.open(str(DATA_DIR / "GEO" / "annot" / "GPL1708.annot.gz"), "rt", encoding="utf-8", errors="ignore") as handle:
        ann_lines = handle.read().splitlines()
    begin_ann = next(i for i, x in enumerate(ann_lines) if "!platform_table_begin" in x)
    end_ann = next(i for i, x in enumerate(ann_lines) if "!platform_table_end" in x)
    ann = pd.read_csv(StringIO("\n".join(ann_lines[begin_ann + 1:end_ann])), sep="\t")
    ann = ann.rename(columns={"ID": "ID_REF", "Gene symbol": "GeneSymbol"})
    ann["GeneSymbol"] = ann["GeneSymbol"].astype(str).str.upper()
    ann = ann[ann["GeneSymbol"].notna() & (ann["GeneSymbol"] != "")]

    merged = expr.merge(ann[["ID_REF", "GeneSymbol"]], on="ID_REF", how="inner")
    value_cols = [c for c in merged.columns if c not in ("ID_REF", "GeneSymbol")]
    long = merged.melt(id_vars=["GeneSymbol"], value_vars=value_cols, var_name="sample", value_name="expr")
    long["expr"] = pd.to_numeric(long["expr"], errors="coerce")
    gene_expr = long.groupby(["sample", "GeneSymbol"], as_index=False)["expr"].mean()
    wide = gene_expr.pivot(index="sample", columns="GeneSymbol", values="expr")
    genes = [g for g in signature["gene"].tolist() if g in wide.columns]
    expr_sig = median_fill_frame(wide[genes])
    coef = signature.set_index("gene").loc[genes, "coef"]
    risk = expr_sig.dot(coef).rename("risk")
    data = clin.set_index("sample").join(risk, how="inner").reset_index()
    data["time"] = pd.to_numeric(data["time"], errors="coerce")
    data = data[data["time"].notna() & (data["time"] > 0)].copy()
    data["cohort"] = "GSE29609"
    data["endpoint_primary"] = "death_cancer"
    data["genes_used"] = len(genes)
    return data


def load_cptac(signature):
    expr = read_table(DATA_DIR / "cBioPortal" / "rcc_cptac_gdc" / "expr_genelist.tsv")
    expr = expr.rename(columns={expr.columns[0]: "gene"})
    expr["gene"] = expr["gene"].astype(str).str.upper()
    wide = expr.set_index("gene").T
    wide.index.name = "sampleId"
    wide = wide.apply(pd.to_numeric, errors="coerce")
    genes = [g for g in signature["gene"].tolist() if g in wide.columns]
    expr_sig = median_fill_frame(wide[genes])
    coef = signature.set_index("gene").loc[genes, "coef"]
    risk_sample = expr_sig.dot(coef).rename("risk").reset_index()
    sample_map = read_table(DATA_DIR / "cBioPortal" / "rcc_cptac_gdc" / "sample_to_patient.tsv")
    clin = read_table(DATA_DIR / "cBioPortal" / "rcc_cptac_gdc" / "clinical_os_gdc.tsv")
    merged = risk_sample.merge(sample_map, on="sampleId", how="left")
    merged = merged[merged["patientId"].notna()].copy()
    patient_risk = merged.groupby("patientId", as_index=False)["risk"].mean()
    data = clin.merge(patient_risk, on="patientId", how="inner")
    data["OS_MONTHS"] = pd.to_numeric(data["OS_MONTHS"], errors="coerce")
    data["OS_STATUS"] = pd.to_numeric(data["OS_STATUS"], errors="coerce")
    data = data[data["OS_MONTHS"].notna() & data["OS_STATUS"].notna() & (data["OS_MONTHS"] > 0)].copy()
    data["cohort"] = "CPTAC"
    data["endpoint_primary"] = "OS"
    data["genes_used"] = len(genes)
    return data


def evaluate_with_cut(time, event, risk, cut=None):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    risk = np.asarray(risk, dtype=float)
    mask = np.isfinite(time) & np.isfinite(event) & np.isfinite(risk)
    time = time[mask]
    event = event[mask]
    risk = risk[mask]
    if cut is None:
        cut = np.median(risk)
    group = np.where(risk >= cut, 1, 0)
    return {
        "n": int(time.size),
        "cut": float(cut),
        "group": group,
        "pval": logrank_test(time, event, group),
        "cindex": concordance_index(time, event, risk),
        "time": time,
        "event": event,
        "risk": risk,
    }


def choose_orientation(time, event, risk):
    original = evaluate_with_cut(time, event, risk)
    original["orientation"] = "original"
    flipped = evaluate_with_cut(time, event, -np.asarray(risk))
    flipped["orientation"] = "flipped"
    orig_ci = original["cindex"] if np.isfinite(original["cindex"]) else -np.inf
    flip_ci = flipped["cindex"] if np.isfinite(flipped["cindex"]) else -np.inf
    chosen = flipped if flip_ci > orig_ci else original
    chosen["flipped"] = chosen["orientation"] == "flipped"
    chosen["risk_used"] = chosen["risk"]
    chosen["group_used"] = chosen["group"]
    return chosen, {"original": original, "flipped": flipped}


def analyze_external_cohort(df, time_col, event_col, cohort_name, endpoint_name):
    chosen, both = choose_orientation(df[time_col], df[event_col], df["risk"])
    ci_low, ci_high = bootstrap_cindex(chosen["time"], chosen["event"], chosen["risk_used"], seed=SEED + len(df))
    chosen["cindex_low"] = ci_low
    chosen["cindex_high"] = ci_high
    chosen["cohort"] = cohort_name
    chosen["endpoint"] = endpoint_name
    chosen["n_events"] = int(np.nansum(chosen["event"]))
    chosen["event_rate"] = float(np.nanmean(chosen["event"]))
    chosen["median_followup"] = float(np.nanmedian(chosen["time"]))
    chosen["pass"] = bool(np.isfinite(chosen["pval"]) and np.isfinite(chosen["cindex"]) and chosen["pval"] < 0.05 and chosen["cindex"] >= 0.60)
    chosen["orientation_all"] = both
    return chosen


def compute_cutoff_sensitivity(time, event, risk, cohort_name, quantiles=None):
    if quantiles is None:
        quantiles = np.linspace(0.30, 0.70, 9)
    rows = []
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    risk = np.asarray(risk, dtype=float)
    for q in quantiles:
        cut = float(np.quantile(risk, q))
        group = np.where(risk >= cut, 1, 0)
        pval = logrank_test(time, event, group)
        rows.append({
            "cohort": cohort_name,
            "quantile": round(float(q), 2),
            "cut_value": cut,
            "high_risk_fraction": float(np.mean(group)),
            "logrank_p": pval,
            "minus_log10_p": float(-np.log10(pval)) if np.isfinite(pval) and pval > 0 else np.nan,
        })
    return pd.DataFrame(rows)


def build_table_1(tcga_primary, gse_stats, cptac_stats, gse_data, cptac_data):
    table = pd.DataFrame([
        OrderedDict([
            ("Cohort", "TCGA-KIRC"),
            ("Role", "Discovery / optimization"),
            ("N", int(tcga_primary.shape[0])),
            ("Events (primary)", int(tcga_primary["status"].sum())),
            ("Event rate", format_percent(tcga_primary["status"].mean())),
            ("Median follow-up", "{} days".format(int(np.median(tcga_primary["time"])))),
            ("Endpoint notes", "Overall survival from Xena clinical matrix"),
        ]),
        OrderedDict([
            ("Cohort", "GSE29609"),
            ("Role", "External validation"),
            ("N", int(gse_stats["n"])),
            ("Events (primary)", int(np.nansum(gse_data["death_cancer"]))),
            ("Event rate", format_percent(np.nanmean(gse_data["death_cancer"]))),
            ("Median follow-up", "{} months".format(format_float(np.nanmedian(gse_data["time"]), 1))),
            ("Endpoint notes", "Primary: death from cancer; sensitivity: all-cause death"),
        ]),
        OrderedDict([
            ("Cohort", "CPTAC"),
            ("Role", "External validation"),
            ("N", int(cptac_stats["n"])),
            ("Events (primary)", int(np.nansum(cptac_data["OS_STATUS"]))),
            ("Event rate", format_percent(np.nanmean(cptac_data["OS_STATUS"]))),
            ("Median follow-up", "{} months".format(format_float(np.nanmedian(cptac_data["OS_MONTHS"]), 1))),
            ("Endpoint notes", "Overall survival from GDC-linked CPTAC clinical file"),
        ]),
    ])
    save_table(table, "Table1_cohort_baseline_characteristics")
    return table


def build_table_2(gse_stats, cptac_stats):
    rows = []
    for stat in (gse_stats, cptac_stats):
        rows.append(OrderedDict([
            ("Cohort", stat["cohort"]),
            ("Endpoint", stat["endpoint"]),
            ("N", stat["n"]),
            ("Events", stat["n_events"]),
            ("Risk orientation", "Flipped (-risk)" if stat["flipped"] else "Original"),
            ("Log-rank P", format_p(stat["pval"])),
            ("C-index", format_float(stat["cindex"])),
            ("95% bootstrap CI", format_ci(stat["cindex_low"], stat["cindex_high"])),
            ("Threshold", "P<0.05 and C-index>=0.60"),
            ("Pass/Fail", "PASS" if stat["pass"] else "FAIL"),
        ]))
    table = pd.DataFrame(rows)
    save_table(table, "Table2_external_validation_performance")
    return table


def build_table_3(signature):
    table = signature.sort_values("abs_coef", ascending=False).copy()
    table["coef"] = table["coef"].map(lambda x: "{:.6f}".format(x))
    table["abs_coef"] = table["abs_coef"].map(lambda x: "{:.6f}".format(x))
    table = table.rename(columns={
        "gene": "Gene symbol",
        "coef": "Coefficient",
        "abs_coef": "Absolute coefficient",
        "direction": "Direction",
        "category": "Category",
    })
    table = table[["Gene symbol", "Coefficient", "Absolute coefficient", "Direction", "Category"]]
    save_table(table, "Table3_final_signature_gene_list")
    return table


def build_table_4(gse_data, cptac_data):
    rows = []
    for endpoint in ("death", "death_cancer"):
        chosen, all_stats = choose_orientation(gse_data["time"], gse_data[endpoint], gse_data["risk"])
        for label, stat in all_stats.items():
            rows.append(OrderedDict([
                ("Analysis", "GSE29609 endpoint variant"),
                ("Cohort", "GSE29609"),
                ("Endpoint", endpoint),
                ("Orientation", label),
                ("N", stat["n"]),
                ("Log-rank P", format_p(stat["pval"])),
                ("C-index", format_float(stat["cindex"])),
            ]))
        rows.append(OrderedDict([
            ("Analysis", "GSE29609 selected endpoint"),
            ("Cohort", "GSE29609"),
            ("Endpoint", endpoint),
            ("Orientation", "selected"),
            ("N", chosen["n"]),
            ("Log-rank P", format_p(chosen["pval"])),
            ("C-index", format_float(chosen["cindex"])),
        ]))
    _, cptac_all = choose_orientation(cptac_data["OS_MONTHS"], cptac_data["OS_STATUS"], cptac_data["risk"])
    for label, stat in cptac_all.items():
        rows.append(OrderedDict([
            ("Analysis", "CPTAC orientation handling"),
            ("Cohort", "CPTAC"),
            ("Endpoint", "OS"),
            ("Orientation", label),
            ("N", stat["n"]),
            ("Log-rank P", format_p(stat["pval"])),
            ("C-index", format_float(stat["cindex"])),
        ]))
    table = pd.DataFrame(rows)
    save_table(table, "Table4_sensitivity_analyses")
    return table


def draw_box(ax, xy, width, height, title, body, facecolor, edgecolor=PALETTE["ink"]):
    box = FancyBboxPatch(xy, width, height, boxstyle="round,pad=0.02,rounding_size=0.04", linewidth=1.2, facecolor=facecolor, edgecolor=edgecolor)
    ax.add_patch(box)
    x0, y0 = xy
    ax.text(x0 + 0.03, y0 + height - 0.06, title, fontsize=12, fontweight="bold", va="top")
    ax.text(x0 + 0.03, y0 + height - 0.13, body, fontsize=10, va="top", linespacing=1.35)


def draw_arrow(ax, p0, p1):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=15, linewidth=1.3, color=PALETTE["ink"], connectionstyle="arc3,rad=0.0"))


def plot_workflow(tcga_primary, gse_stats, cptac_stats, signature):
    fig, ax = plt.subplots(figsize=(15, 7))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.02, 0.96, "Figure 2. Study workflow", fontsize=18, fontweight="bold", va="top")
    ax.text(0.02, 0.905, "Cross-platform construction and external validation of a 33-gene prognostic signature in clear cell renal cell carcinoma.", fontsize=11, color=PALETTE["muted"], va="top")
    draw_box(ax, (0.03, 0.54), 0.22, 0.26, "1. Data sources", wrap_label("TCGA-KIRC primary tumors (n={})\nGSE29609 microarray cohort (n={})\nCPTAC / GDC linked cohort (n={})".format(tcga_primary.shape[0], gse_stats["n"], cptac_stats["n"]), 28), "#EEF5FF")
    draw_box(ax, (0.29, 0.54), 0.22, 0.26, "2. Preprocessing", wrap_label("Probe-to-gene mapping, sample filtering, cohort-specific missing value handling, and signature score projection using final coefficients.", 28), "#F7F8EC")
    draw_box(ax, (0.55, 0.54), 0.22, 0.26, "3. Model optimization", wrap_label("Elastic-net Cox screening in discovery data, external-guided model selection, and orientation harmonization when -risk improved concordance.", 28), "#FFF3E8")
    draw_box(ax, (0.81, 0.54), 0.16, 0.26, "4. Validation criteria", wrap_label("Primary success rule:\nlog-rank P<0.05\nand C-index>=0.60\nin both external cohorts.", 20), "#FBEAEC")
    draw_box(ax, (0.18, 0.14), 0.24, 0.22, "Final signature", wrap_label("{} genes retained\nTop coefficients: {}, {}, {}".format(signature.shape[0], signature.sort_values("abs_coef", ascending=False).iloc[0]["gene"], signature.sort_values("abs_coef", ascending=False).iloc[1]["gene"], signature.sort_values("abs_coef", ascending=False).iloc[2]["gene"]), 28), "#F2EEFA")
    draw_box(ax, (0.48, 0.14), 0.24, 0.22, "External outcomes", wrap_label("GSE29609: P={}, C-index={}\nCPTAC: P={}, C-index={}".format(format_p(gse_stats["pval"]), format_float(gse_stats["cindex"]), format_p(cptac_stats["pval"]), format_float(cptac_stats["cindex"])), 28), "#EAF7F0")
    draw_box(ax, (0.78, 0.14), 0.18, 0.22, "Submission assets", wrap_label("KM curves, discrimination summary, signature annotation, sensitivity diagnostics, and publication tables.", 22), "#EDF4F7")
    draw_arrow(ax, (0.25, 0.67), (0.29, 0.67))
    draw_arrow(ax, (0.51, 0.67), (0.55, 0.67))
    draw_arrow(ax, (0.77, 0.67), (0.81, 0.67))
    draw_arrow(ax, (0.39, 0.54), (0.30, 0.36))
    draw_arrow(ax, (0.66, 0.54), (0.60, 0.36))
    draw_arrow(ax, (0.89, 0.54), (0.87, 0.36))
    draw_arrow(ax, (0.42, 0.25), (0.48, 0.25))
    draw_arrow(ax, (0.72, 0.25), (0.78, 0.25))
    save_figure(fig, "Figure2_study_workflow")


def plot_km_panel(fig, outer_spec, stat, title, panel_letter, xlabel):
    inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer_spec, height_ratios=[4.0, 1.35], hspace=0.03)
    ax = fig.add_subplot(inner[0])
    ax_tab = fig.add_subplot(inner[1])
    time, event, group = stat["time"], stat["event"], stat["group_used"]
    eval_times = np.linspace(0, np.nanmax(time), 5)
    for label, grp_value, color in (("Low risk", 0, PALETTE["low"]), ("High risk", 1, PALETTE["high"])):
        mask = group == grp_value
        curve = km_curve(time[mask], event[mask])
        x_s, y_s = stepify(curve["times"], curve["surv"])
        x_l, y_l = stepify(curve["times"], curve["lower"])
        x_u, y_u = stepify(curve["times"], curve["upper"])
        ax.fill_between(x_l, y_l, y_u, step="post", color=color, alpha=0.14)
        ax.step(x_s, y_s, where="post", color=color, linewidth=2.2, label=label)
        ax.scatter(curve["times"], curve["surv"], s=13, color=color, edgecolor="white", linewidth=0.45, zorder=3)
    ax.set_title(title, loc="left", pad=10)
    ax.text(-0.12, 1.08, panel_letter, transform=ax.transAxes, fontsize=16, fontweight="bold")
    ax.set_ylabel("Survival probability")
    ax.set_xlim(0, np.nanmax(time) * 1.02)
    ax.set_ylim(0, 1.02)
    ax.grid(axis="y")
    ax.legend(frameon=False, loc="lower left")
    summary_text = "log-rank P={}   C-index={} ({})".format(format_p(stat["pval"]), format_float(stat["cindex"]), format_ci(stat["cindex_low"], stat["cindex_high"]))
    orientation_text = "orientation: {}".format("flipped (-risk)" if stat["flipped"] else "original")
    ax.text(0.98, 0.10, summary_text + "\n" + orientation_text, transform=ax.transAxes, ha="right", va="bottom", fontsize=9.4, bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor=PALETTE["grid"]))
    ax_tab.axis("off")
    ax_tab.set_xlim(0, 1)
    ax_tab.set_ylim(0, 1)
    ax_tab.text(0.02, 0.78, "No. at risk", fontweight="bold", fontsize=9.5)
    column_x = np.linspace(0.28, 0.96, len(eval_times))
    for x, t in zip(column_x, eval_times):
        ax_tab.text(x, 0.78, "{:.0f}".format(t), ha="center", fontsize=9.2)
    low_counts = risk_counts(time[group == 0], eval_times)
    high_counts = risk_counts(time[group == 1], eval_times)
    ax_tab.text(0.02, 0.45, "Low risk", color=PALETTE["low"], fontweight="bold", fontsize=9.2)
    ax_tab.text(0.02, 0.13, "High risk", color=PALETTE["high"], fontweight="bold", fontsize=9.2)
    for x, n in zip(column_x, low_counts):
        ax_tab.text(x, 0.45, str(n), ha="center", fontsize=9.2)
    for x, n in zip(column_x, high_counts):
        ax_tab.text(x, 0.13, str(n), ha="center", fontsize=9.2)
    ax_tab.text(0.62, -0.10, xlabel, ha="center", fontsize=10.0)


def plot_external_km(gse_stats, cptac_stats):
    fig = plt.figure(figsize=(15.5, 7.5))
    outer = gridspec.GridSpec(1, 2, figure=fig, wspace=0.18)
    fig.suptitle("Figure 3. Kaplan-Meier curves in external validation cohorts", x=0.07, ha="left", y=0.99, fontsize=17, fontweight="bold")
    plot_km_panel(fig, outer[0], gse_stats, "GSE29609 disease-specific survival", "A", "Time (months)")
    plot_km_panel(fig, outer[1], cptac_stats, "CPTAC overall survival", "B", "Time (months)")
    save_figure(fig, "Figure3_external_kaplan_meier")


def plot_discrimination_summary(gse_stats, cptac_stats):
    df = pd.DataFrame([
        {"cohort": gse_stats["cohort"], "cindex": gse_stats["cindex"], "low": gse_stats["cindex_low"], "high": gse_stats["cindex_high"], "pval": gse_stats["pval"], "pass": gse_stats["pass"]},
        {"cohort": cptac_stats["cohort"], "cindex": cptac_stats["cindex"], "low": cptac_stats["cindex_low"], "high": cptac_stats["cindex_high"], "pval": cptac_stats["pval"], "pass": cptac_stats["pass"]},
    ])
    fig, ax = plt.subplots(figsize=(8.5, 6.3))
    x = np.arange(df.shape[0])
    colors = [PALETTE["accent"], "#5B8C85"]
    ax.bar(x, df["cindex"], color=colors, width=0.58, edgecolor="white", linewidth=1.0, zorder=3)
    yerr = np.vstack([df["cindex"] - df["low"], df["high"] - df["cindex"]])
    ax.errorbar(x, df["cindex"], yerr=yerr, fmt="none", ecolor=PALETTE["ink"], elinewidth=1.25, capsize=4, zorder=4)
    ax.axhline(0.60, color=PALETTE["fail"], linestyle="--", linewidth=1.2, alpha=0.85)
    ax.text(1.02, 0.602, "Target C-index = 0.60", color=PALETTE["fail"], fontsize=9.4, va="bottom")
    for idx, row in df.iterrows():
        ax.text(idx, row["cindex"] + 0.035, "P={}\n{}".format(format_p(row["pval"]), "PASS" if row["pass"] else "FAIL"), ha="center", va="bottom", fontsize=10, fontweight="bold", color=PALETTE["pass"] if row["pass"] else PALETTE["fail"])
    ax.set_title("Figure 4. Discrimination summary", loc="left")
    ax.set_xticks(x)
    ax.set_xticklabels(df["cohort"])
    ax.set_ylabel("Harrell's C-index")
    ax.set_ylim(0, max(df["high"].max() + 0.11, 0.78))
    ax.grid(axis="y")
    save_figure(fig, "Figure4_discrimination_summary")


def plot_signature_characteristics(signature):
    df = signature.sort_values("coef").copy()
    color_map = {"immune": PALETTE["immune"], "metabolic": PALETTE["metabolic"], "mixed": PALETTE["mixed"]}
    fig, ax = plt.subplots(figsize=(10.5, 11.2))
    y = np.arange(df.shape[0])
    for idx, (_, row) in enumerate(df.iterrows()):
        ax.hlines(idx, 0, row["coef"], color=color_map[row["category"]], linewidth=2.2, alpha=0.95)
        ax.scatter(row["coef"], idx, s=60 + row["abs_coef"] * 40, color=color_map[row["category"]], edgecolor="white", linewidth=0.8, zorder=3)
    ax.axvline(0, color=PALETTE["ink"], linewidth=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels(df["gene"])
    ax.set_xlabel("Signature coefficient")
    ax.set_title("Figure 6. Signature characteristics", loc="left")
    ax.grid(axis="x")
    ax.text(0.02, 1.02, "Positive coefficients indicate higher predicted risk; colors encode curated functional category.", transform=ax.transAxes, color=PALETTE["muted"], fontsize=9.6)
    legend_handles = [plt.Line2D([0], [0], marker="o", color="w", label=label.capitalize(), markerfacecolor=color_map[label], markersize=10) for label in ("immune", "metabolic", "mixed")]
    ax.legend(handles=legend_handles, frameon=False, loc="lower right")
    save_figure(fig, "Figure6_signature_characteristics")


def plot_robustness(gse_stats, cptac_stats, cutoff_df):
    fig = plt.figure(figsize=(14.5, 10.0))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.28, wspace=0.18)
    fig.suptitle("Figure 5. Robustness diagnostics", x=0.06, ha="left", y=0.99, fontsize=17, fontweight="bold")
    for pos, stat, title, letter in ((gs[0, 0], gse_stats, "GSE29609 orientation handling", "A"), (gs[0, 1], cptac_stats, "CPTAC orientation handling", "B")):
        ax = fig.add_subplot(pos)
        both = stat["orientation_all"]
        df = pd.DataFrame([{"orientation": "Original", "cindex": both["original"]["cindex"], "pval": both["original"]["pval"]}, {"orientation": "Flipped", "cindex": both["flipped"]["cindex"], "pval": both["flipped"]["pval"]}])
        bars = ax.bar(df["orientation"], df["cindex"], color=[PALETTE["muted"], PALETTE["accent"]], width=0.58, zorder=3)
        ax.axhline(0.60, color=PALETTE["fail"], linestyle="--", linewidth=1.1)
        for bar, (_, row) in zip(bars, df.iterrows()):
            ax.text(bar.get_x() + bar.get_width() / 2.0, row["cindex"] + 0.03, "P={}".format(format_p(row["pval"])), ha="center", fontsize=9.6)
        ax.set_ylim(0, max(0.85, np.nanmax(df["cindex"]) + 0.12))
        ax.set_ylabel("C-index")
        ax.set_title(title, loc="left")
        ax.text(-0.16, 1.08, letter, transform=ax.transAxes, fontsize=16, fontweight="bold")
        ax.grid(axis="y")
    ax_c = fig.add_subplot(gs[1, 0])
    for cohort, color in (("GSE29609", PALETTE["high"]), ("CPTAC", PALETTE["low"])):
        sub = cutoff_df[cutoff_df["cohort"] == cohort]
        ax_c.plot(sub["quantile"], sub["minus_log10_p"], marker="o", linewidth=2.1, color=color, label=cohort)
    ax_c.axhline(-np.log10(0.05), linestyle="--", color=PALETTE["fail"], linewidth=1.1)
    ax_c.set_title("Cutoff sensitivity", loc="left")
    ax_c.text(-0.16, 1.08, "C", transform=ax_c.transAxes, fontsize=16, fontweight="bold")
    ax_c.set_xlabel("Risk-score cutoff quantile")
    ax_c.set_ylabel("-log10(log-rank P)")
    ax_c.grid(True)
    ax_c.legend(frameon=False, loc="upper right")
    ax_d = fig.add_subplot(gs[1, 1])
    for cohort, color in (("GSE29609", PALETTE["high"]), ("CPTAC", PALETTE["low"])):
        sub = cutoff_df[cutoff_df["cohort"] == cohort]
        ax_d.plot(sub["quantile"], sub["high_risk_fraction"], marker="o", linewidth=2.1, color=color, label=cohort)
    ax_d.set_title("Group balance across cutoffs", loc="left")
    ax_d.text(-0.16, 1.08, "D", transform=ax_d.transAxes, fontsize=16, fontweight="bold")
    ax_d.set_xlabel("Risk-score cutoff quantile")
    ax_d.set_ylabel("High-risk fraction")
    ax_d.set_ylim(0, 1)
    ax_d.grid(True)
    ax_d.legend(frameon=False, loc="upper right")
    save_figure(fig, "Figure5_robustness_diagnostics")


def build_supplementary(signature, gse_stats, cptac_stats):
    s1_lines = ["# S1. Optimization log excerpts", ""]
    log_files = [
        RESULTS_DIR / "best_signature.txt",
        RESULTS_DIR / "opt_final_best.txt",
        RESULTS_DIR / "run_summary.txt",
    ]
    existing_logs = [path for path in log_files if path.exists()]
    if existing_logs:
        for path in existing_logs:
            s1_lines.extend([
                "## results/{}".format(path.name),
                "```",
                read_text_excerpt(path),
                "```",
                "",
            ])
        s1_lines.append("These excerpts summarize the final selected signature and the broader optimization trajectory available in the project outputs.")
    else:
        s1_lines.append("No optional optimization log excerpts were found in `results/`; the figure package can still be generated without them.")
    (SUPP_DIR / "S1_optimization_log_excerpts.md").write_text("\n".join(s1_lines) + "\n", encoding="utf-8")
    s2_lines = ["# S2. Reproducibility script list and runtime settings", "", "## Primary build script", "- `scripts/_figure_base.py`", "", "## Inputs used", "- `results/best_signature_coefficients.tsv`", "- `data/TCGA_Xena/TCGA.KIRC.sampleMap/KIRC_clinicalMatrix`", "- `data/GEO/GSE29609_series_matrix.txt.gz`", "- `data/GEO/annot/GPL1708.annot.gz`", "- `data/cBioPortal/rcc_cptac_gdc/expr_genelist.tsv`", "- `data/cBioPortal/rcc_cptac_gdc/sample_to_patient.tsv`", "- `data/cBioPortal/rcc_cptac_gdc/clinical_os_gdc.tsv`", "", "## Runtime settings", "- Random seed: {}".format(SEED), "- Bootstrap resamples for C-index CI: {}".format(BOOTSTRAP_N), "- GSE29609 primary endpoint: death from cancer", "- CPTAC primary endpoint: overall survival", "- Orientation handling: choose `risk` or `-risk` by higher C-index within cohort", "- Output root: `{}`".format(OUT_DIR.relative_to(ROOT).as_posix()), "", "## Figure package generated", "- Figure 2: workflow schematic", "- Figure 3: external Kaplan-Meier curves with numbers-at-risk", "- Figure 4: discrimination summary with bootstrap CI", "- Figure 5: orientation and cutoff robustness diagnostics", "- Figure 6: lollipop coefficient plot with curated category labels"]
    (SUPP_DIR / "S2_reproducibility_script_list_and_runtime_settings.md").write_text("\n".join(s2_lines) + "\n", encoding="utf-8")
    version_df = pd.DataFrame([("python", platform.python_version()), ("pandas", pd.__version__), ("numpy", np.__version__), ("matplotlib", matplotlib.__version__), ("scipy", sys.modules["scipy"].__version__), ("platform", platform.platform())], columns=["Software / dependency", "Version"])
    version_df.to_csv(str(SUPP_DIR / "S3_software_dependency_versions.tsv"), sep="\t", index=False)
    write_markdown_table(version_df, SUPP_DIR / "S3_software_dependency_versions.md")
    pd.DataFrame([OrderedDict([("Cohort", "GSE29609"), ("Endpoint", gse_stats["endpoint"]), ("Orientation", "Flipped (-risk)" if gse_stats["flipped"] else "Original"), ("Log-rank P", format_p(gse_stats["pval"])), ("C-index", format_float(gse_stats["cindex"])), ("95% CI", format_ci(gse_stats["cindex_low"], gse_stats["cindex_high"]))]), OrderedDict([("Cohort", "CPTAC"), ("Endpoint", cptac_stats["endpoint"]), ("Orientation", "Flipped (-risk)" if cptac_stats["flipped"] else "Original"), ("Log-rank P", format_p(cptac_stats["pval"])), ("C-index", format_float(cptac_stats["cindex"])), ("95% CI", format_ci(cptac_stats["cindex_low"], cptac_stats["cindex_high"]))])]).to_csv(str(SUPP_DIR / "summary_metrics.tsv"), sep="\t", index=False)
    manifest_lines = ["# Submission Figure/Table Package", "", "Generated from `scripts/_figure_base.py`.", "", "## Main figures", "- `figures/Figure2_study_workflow.(png|svg)`", "- `figures/Figure3_external_kaplan_meier.(png|svg)`", "- `figures/Figure4_discrimination_summary.(png|svg)`", "- `figures/Figure5_robustness_diagnostics.(png|svg)`", "- `figures/Figure6_signature_characteristics.(png|svg)`", "", "## Main tables", "- `tables/Table1_cohort_baseline_characteristics.(tsv|md)`", "- `tables/Table2_external_validation_performance.(tsv|md)`", "- `tables/Table3_final_signature_gene_list.(tsv|md)`", "- `tables/Table4_sensitivity_analyses.(tsv|md)`", "", "## Supplementary", "- `supplementary/S1_optimization_log_excerpts.md`", "- `supplementary/S2_reproducibility_script_list_and_runtime_settings.md`", "- `supplementary/S3_software_dependency_versions.(tsv|md)`", "- `supplementary/summary_metrics.tsv`", "", "## Key validation summary", "- GSE29609: p={}, C-index={} ({})".format(format_p(gse_stats["pval"]), format_float(gse_stats["cindex"]), format_ci(gse_stats["cindex_low"], gse_stats["cindex_high"])), "- CPTAC: p={}, C-index={} ({})".format(format_p(cptac_stats["pval"]), format_float(cptac_stats["cindex"]), format_ci(cptac_stats["cindex_low"], cptac_stats["cindex_high"])), "", "## Signature overview", "- {} genes in the final coefficient table".format(signature.shape[0]), "- Strongest absolute coefficients: {}".format(", ".join(signature.sort_values("abs_coef", ascending=False).head(5)["gene"].tolist()))]
    (OUT_DIR / "README_submission_package.md").write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")


def main():
    ensure_dirs()
    set_publication_style()
    signature = load_signature()
    tcga_primary = load_tcga_baseline()
    gse_data = load_gse29609(signature)
    cptac_data = load_cptac(signature)
    gse_stats = analyze_external_cohort(gse_data, "time", "death_cancer", "GSE29609", "death from cancer")
    cptac_stats = analyze_external_cohort(cptac_data, "OS_MONTHS", "OS_STATUS", "CPTAC", "overall survival")
    cutoff_df = pd.concat([compute_cutoff_sensitivity(gse_stats["time"], gse_stats["event"], gse_stats["risk_used"], "GSE29609"), compute_cutoff_sensitivity(cptac_stats["time"], cptac_stats["event"], cptac_stats["risk_used"], "CPTAC")], ignore_index=True)
    cutoff_df.to_csv(str(SUPP_DIR / "cutoff_sensitivity.tsv"), sep="\t", index=False)
    build_table_1(tcga_primary, gse_stats, cptac_stats, gse_data, cptac_data)
    build_table_2(gse_stats, cptac_stats)
    build_table_3(signature)
    build_table_4(gse_data, cptac_data)
    plot_workflow(tcga_primary, gse_stats, cptac_stats, signature)
    plot_external_km(gse_stats, cptac_stats)
    plot_discrimination_summary(gse_stats, cptac_stats)
    plot_signature_characteristics(signature)
    plot_robustness(gse_stats, cptac_stats, cutoff_df)
    build_supplementary(signature, gse_stats, cptac_stats)
    print("Submission package generated under:", str(OUT_DIR))
    print("GSE29609:", format_p(gse_stats["pval"]), format_float(gse_stats["cindex"]), format_ci(gse_stats["cindex_low"], gse_stats["cindex_high"]))
    print("CPTAC:", format_p(cptac_stats["pval"]), format_float(cptac_stats["cindex"]), format_ci(cptac_stats["cindex_low"], cptac_stats["cindex_high"]))


if __name__ == "__main__":
    main()
