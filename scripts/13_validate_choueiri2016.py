#!/usr/bin/env python

import csv
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index


ROOT = Path(__file__).resolve().parents[1]
DIR_DATA = ROOT / "data"
DIR_RES = ROOT / "results"
DIR_EXT = DIR_DATA / "external_clinical" / "Choueiri2016"

BASE = "https://www.cbioportal.org/api"
STUDY = "ccrcc_iatlas_choueiri_2016"
UA = "Mozilla/5.0"
H_JSON = {"Accept": "application/json", "User-Agent": UA}
H_POST = {"Accept": "application/json", "Content-Type": "application/json", "User-Agent": UA}


def cbio_get(session, path, **params):
    resp = session.get(f"{BASE}{path}", headers=H_JSON, params=params, timeout=180)
    resp.raise_for_status()
    return resp.json()


def cbio_post(session, path, payload):
    resp = session.post(f"{BASE}{path}", headers=H_POST, data=json.dumps(payload), timeout=240)
    resp.raise_for_status()
    return resp.json()


def list_all(session, path, **params):
    out = []
    page = 0
    while True:
        arr = cbio_get(session, path, pageSize=5000, pageNumber=page, **params)
        if not arr:
            break
        out.extend(arr)
        page += 1
        time.sleep(0.1)
    return out


def load_gene_map(session):
    rows = list_all(session, "/genes")
    out = {}
    for rec in rows:
        hugo = str(rec.get("hugoGeneSymbol", "")).upper().strip()
        eid = rec.get("entrezGeneId")
        if hugo and eid is not None and hugo not in out:
            out[hugo] = int(eid)
    return out


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


def parse_status(x):
    x = str(x).strip()
    if not x or x.lower() == "nan":
        return np.nan
    if x.startswith("1:"):
        return 1
    if x.startswith("0:"):
        return 0
    return np.nan


def endpoint_summary(cohort_name, study_id, study_name, endpoint_name, time, status, risk, n_genes_requested, n_genes_used, missing_genes):
    metrics = orientation_metrics(time.to_numpy(), status.to_numpy(), risk.to_numpy())
    original = metrics["original"]
    return {
        "cohort": cohort_name,
        "study_id": study_id,
        "study_name": study_name,
        "endpoint": endpoint_name,
        "n_samples": int(np.isfinite(pd.to_numeric(time, errors="coerce")).sum()),
        "n_events": int(np.nansum(status)),
        "signature_genes_requested": int(n_genes_requested),
        "signature_genes_used": int(n_genes_used),
        "missing_signature_genes": ";".join(missing_genes),
        "manuscript_orientation": "original",
        "original_logrank_p": original["pval"],
        "original_cindex": original["cindex"],
        "sign_reversed_logrank_p": metrics["sign_reversed"]["pval"],
        "sign_reversed_cindex": metrics["sign_reversed"]["cindex"],
        "best_discrimination_orientation": metrics["best_discrimination_orientation"],
        "logrank_p": original["pval"],
        "cindex": original["cindex"],
        "median_cut_high_n": int((risk >= np.nanmedian(risk)).sum()),
        "median_cut_low_n": int((risk < np.nanmedian(risk)).sum()),
    }


def main():
    DIR_EXT.mkdir(parents=True, exist_ok=True)
    DIR_RES.mkdir(parents=True, exist_ok=True)

    coef = pd.read_csv(DIR_RES / "best_signature_coefficients.tsv", sep="\t")
    coef["gene"] = coef["gene"].astype(str).str.upper()
    coef = coef.set_index("gene")["coef"]
    signature_genes = list(coef.index)

    session = requests.Session()

    study = cbio_get(session, f"/studies/{STUDY}")
    profiles = cbio_get(session, f"/studies/{STUDY}/molecular-profiles")
    expr_profile = None
    for prof in profiles:
        if prof.get("molecularAlterationType") == "MRNA_EXPRESSION" and prof.get("datatype") == "CONTINUOUS":
            expr_profile = str(prof.get("molecularProfileId", ""))
            break
    if expr_profile is None:
        raise SystemExit("No continuous RNA expression profile found")

    samples = cbio_get(session, f"/studies/{STUDY}/samples")
    sample_df = pd.DataFrame(samples)[["sampleId", "patientId", "sampleType"]].copy()
    sample_df = sample_df.rename(columns={"sampleType": "sample_type_api"})

    clinical_patient = pd.DataFrame(list_all(session, f"/studies/{STUDY}/clinical-data", projection="SUMMARY", clinicalDataType="PATIENT"))
    patient_wide = clinical_patient.pivot_table(index="patientId", columns="clinicalAttributeId", values="value", aggfunc="first")
    clin = sample_df.merge(patient_wide, left_on="patientId", right_index=True, how="left")

    clin["os_months"] = pd.to_numeric(clin["OS_MONTHS"], errors="coerce")
    clin["status_os"] = clin["OS_STATUS"].map(parse_status)
    clin["pfs_months"] = pd.to_numeric(clin["PFS_MONTHS"], errors="coerce")
    clin["status_pfs"] = clin["PFS_STATUS"].map(parse_status)

    gene_map = load_gene_map(session)
    kept_genes = [g for g in signature_genes if g in gene_map]
    missing_genes = [g for g in signature_genes if g not in gene_map]
    entrez_ids = [gene_map[g] for g in kept_genes]

    sample_ids = sample_df["sampleId"].tolist()
    payload = {"sampleIds": sample_ids, "entrezGeneIds": entrez_ids}
    expr_data = cbio_post(session, f"/molecular-profiles/{expr_profile}/molecular-data/fetch", payload)

    expr_mat = {g: [""] * len(sample_ids) for g in kept_genes}
    sidx = {sid: i for i, sid in enumerate(sample_ids)}
    eid_to_gene = {eid: g for eid, g in zip(entrez_ids, kept_genes)}
    for rec in expr_data:
        sid = rec.get("sampleId")
        eid = rec.get("entrezGeneId")
        val = rec.get("value")
        i = sidx.get(sid)
        g = eid_to_gene.get(eid)
        if i is not None and g is not None:
            expr_mat[g][i] = "" if val is None else str(val)

    expr_export = pd.DataFrame({"gene": kept_genes})
    for sid in sample_ids:
        expr_export[sid] = [expr_mat[g][sidx[sid]] for g in kept_genes]
    expr_export.to_csv(DIR_EXT / "expression_signature_gene33.tsv", sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)

    keep_cols = [
        "sampleId",
        "patientId",
        "sample_type_api",
        "AGE_AT_DIAGNOSIS",
        "SEX",
        "CLINICAL_STAGE",
        "ICI_PATHWAY",
        "ICI_RX",
        "ICI_TARGET",
        "CLINICAL_BENEFIT",
        "RESPONDER",
        "PROGRESSION",
        "OS_MONTHS",
        "OS_STATUS",
        "os_months",
        "status_os",
        "PFS_MONTHS",
        "PFS_STATUS",
        "pfs_months",
        "status_pfs",
    ]
    clin_out = clin[keep_cols].copy().rename(columns={"sampleId": "sample_id", "patientId": "patient_id"})
    clin_out.to_csv(DIR_EXT / "clinical_matched_16.tsv", sep="\t", index=False)

    expr = expr_export.copy()
    expr["gene"] = expr["gene"].astype(str).str.upper()
    expr = expr.set_index("gene").T
    expr = expr.apply(pd.to_numeric, errors="coerce")

    common = expr.index.intersection(clin_out["sample_id"])
    expr = expr.loc[common]
    clin_use = clin_out.set_index("sample_id").loc[common].reset_index()

    expr_r = rank_rows(expr[kept_genes])
    risk_raw = expr_r.dot(coef.loc[kept_genes])

    os_summary = endpoint_summary(
        "Choueiri2016",
        STUDY,
        study["name"],
        "OS",
        clin_use["os_months"],
        clin_use["status_os"],
        risk_raw,
        len(signature_genes),
        len(kept_genes),
        missing_genes,
    )
    pfs_summary = endpoint_summary(
        "Choueiri2016",
        STUDY,
        study["name"],
        "PFS",
        clin_use["pfs_months"],
        clin_use["status_pfs"],
        risk_raw,
        len(signature_genes),
        len(kept_genes),
        missing_genes,
    )

    score_df = clin_use.copy()
    score_df["risk_score_raw"] = risk_raw
    score_df["risk_score_sign_reversed"] = -risk_raw
    score_df.to_csv(DIR_RES / "external_validation_choueiri2016_scores.tsv", sep="\t", index=False)

    summary_df = pd.DataFrame([os_summary, pfs_summary])
    summary_df.to_csv(DIR_RES / "external_validation_choueiri2016_summary.tsv", sep="\t", index=False)
    (DIR_RES / "external_validation_choueiri2016_summary.json").write_text(
        json.dumps({"endpoints": [os_summary, pfs_summary]}, indent=2),
        encoding="utf-8",
    )

    readme = f"""# Choueiri 2016 external validation cache

Recovered on: 2026-03-25

## Public sources

- cBioPortal study page: https://www.cbioportal.org/study/summary?id={STUDY}
- cBioPortal API study endpoint: https://www.cbioportal.org/api/studies/{STUDY}
- PMID: https://pubmed.ncbi.nlm.nih.gov/27169994/

## What is stored here

- `clinical_matched_16.tsv`: sample-matched clinical table exported from cBioPortal.
- `expression_signature_gene33.tsv`: public RNA-seq expression for the frozen 33-gene signature.

## Notes

- This cohort is a very small metastatic clear-cell RCC nivolumab trial cohort.
- Both `OS` and `PFS` endpoints are exposed by cBioPortal and are evaluated from the same frozen signature score.
- Cancer Medicine summaries retain the original score direction; sign-reversed metrics are supplementary robustness diagnostics only.
- Because `n = 16`, this cohort should be treated as a supplementary sensitivity validation only.
"""
    (DIR_EXT / "README.md").write_text(readme, encoding="utf-8")

    print(json.dumps({"endpoints": [os_summary, pfs_summary]}, indent=2))


if __name__ == "__main__":
    main()
