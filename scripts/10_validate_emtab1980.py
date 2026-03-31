#!/usr/bin/env python

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index


ROOT = Path(__file__).resolve().parents[1]
DIR_DATA = ROOT / "data"
DIR_RES = ROOT / "results"
DIR_EXT = DIR_DATA / "external_clinical" / "E-MTAB-1980"


def rank_rows(df):
    ranks = df.rank(axis=1, method="average", na_option="keep")
    max_rank = ranks.max(axis=1)
    out = ranks.div(max_rank, axis=0)
    return out.fillna(0.5)


def km_stats(time, status, risk):
    mask = np.isfinite(time) & pd.notna(status) & np.isfinite(risk)
    time = np.asarray(time[mask], dtype=float)
    status = np.asarray(status[mask], dtype=int)
    risk = np.asarray(risk[mask], dtype=float)
    n0 = len(risk)
    if n0 < 10 or len(set(status)) < 2:
        return {"pval": np.nan, "cindex": np.nan, "n": n0}
    cut0 = np.nanmedian(risk)
    grp = risk >= cut0
    if len(set(grp)) < 2:
        return {"pval": np.nan, "cindex": np.nan, "n": n0}
    try:
        lr = logrank_test(time[grp], time[~grp], event_observed_A=status[grp], event_observed_B=status[~grp])
        pval = float(lr.p_value)
    except Exception:
        pval = np.nan
    try:
        cidx = float(concordance_index(time, -risk, status))
    except Exception:
        cidx = np.nan
    return {"pval": pval, "cindex": cidx, "n": n0}


def orientation_metrics(time, status, risk):
    original = km_stats(time, status, risk)
    sign_reversed = km_stats(time, status, -risk)
    c1 = original["cindex"] if np.isfinite(original["cindex"]) else -np.inf
    c2 = sign_reversed["cindex"] if np.isfinite(sign_reversed["cindex"]) else -np.inf
    return {
        "original": original,
        "sign_reversed": sign_reversed,
        "best_discrimination_orientation": "sign_reversed" if c2 > c1 else "original",
    }


def main():
    coef = pd.read_csv(DIR_RES / "best_signature_coefficients.tsv", sep="\t")
    coef["gene"] = coef["gene"].astype(str).str.upper()
    coef = coef.set_index("gene")["coef"]

    expr = pd.read_csv(DIR_EXT / "expression_signature_gene33.tsv", sep="\t")
    expr["gene"] = expr["gene"].astype(str).str.upper()
    expr = expr.set_index("gene").T
    expr = expr.apply(pd.to_numeric, errors="coerce")

    clin = pd.read_csv(DIR_EXT / "clinical_matched_101.tsv", sep="\t")
    common = expr.index.intersection(clin["sample_id"])
    expr = expr.loc[common]
    clin = clin.set_index("sample_id").loc[common].reset_index().rename(columns={"index": "sample_id"})

    expr_r = rank_rows(expr[coef.index])
    risk = expr_r.dot(coef)
    metrics = orientation_metrics(clin["time_months"].to_numpy(), clin["status_os"].to_numpy(), risk.to_numpy())
    original = metrics["original"]
    sign_reversed = metrics["sign_reversed"]

    summary = {
        "cohort": "E-MTAB-1980",
        "n_samples": int(len(clin)),
        "n_events": int(clin["status_os"].sum()),
        "candidate_genes_for_ranking": int(expr_r.shape[1]),
        "signature_genes_used": int(sum(g in expr_r.columns for g in coef.index)),
        "manuscript_orientation": "original",
        "original_logrank_p": original["pval"],
        "original_cindex": original["cindex"],
        "sign_reversed_logrank_p": sign_reversed["pval"],
        "sign_reversed_cindex": sign_reversed["cindex"],
        "best_discrimination_orientation": metrics["best_discrimination_orientation"],
        "logrank_p": original["pval"],
        "cindex": original["cindex"],
        "median_cut_high_n": int((risk >= np.nanmedian(risk)).sum()),
        "median_cut_low_n": int((risk < np.nanmedian(risk)).sum()),
    }

    scores = pd.DataFrame(
        {
            "sample_id": clin["sample_id"],
            "time_months": clin["time_months"],
            "status_os": clin["status_os"],
            "risk_score": risk,
            "risk_score_sign_reversed": -risk,
        }
    )

    DIR_RES.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([summary]).to_csv(DIR_RES / "external_validation_emtab1980_summary.tsv", sep="\t", index=False)
    scores.to_csv(DIR_RES / "external_validation_emtab1980_scores.tsv", sep="\t", index=False)
    (DIR_RES / "external_validation_emtab1980_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
