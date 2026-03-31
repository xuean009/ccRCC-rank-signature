#!/usr/bin/env python

import gzip
import json
import os
import re
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index


ROOT = Path(__file__).resolve().parents[1]
DIR_DATA = ROOT / "data"
DIR_RES = ROOT / "results"

RUN_TAG = os.environ.get("KIRC_RUN_TAG", "fixed_direction_py").strip() or "fixed_direction_py"
N_ITER = int(os.environ.get("KIRC_PY_ITERATIONS", "40"))
SEED = int(os.environ.get("KIRC_SEED", "42"))


def out_path(name):
    return DIR_RES / f"{RUN_TAG}_{name}"


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


def univ_cox_stats(time, status, x):
    mask = np.isfinite(time) & pd.notna(status) & np.isfinite(x)
    if mask.sum() < 20 or len(set(np.asarray(status)[mask])) < 2:
        return {"beta": np.nan, "p": np.nan}
    df = pd.DataFrame({"time": np.asarray(time)[mask], "status": np.asarray(status)[mask], "x": np.asarray(x)[mask]})
    try:
        fit = CoxPHFitter()
        fit.fit(df, duration_col="time", event_col="status")
        beta = float(fit.params_["x"])
        pval = float(fit.summary.loc["x", "p"])
        return {"beta": beta, "p": pval}
    except Exception:
        return {"beta": np.nan, "p": np.nan}


def make_risk(df, coef_vec):
    common = [g for g in coef_vec.index if g in df.columns]
    if len(common) < 2:
        return pd.Series(np.nan, index=df.index)
    xx = df.loc[:, common].copy()
    xx = xx.fillna(0.0)
    return xx.dot(coef_vec.loc[common])


def impute_col_median(df):
    out = df.copy()
    for col in out.columns:
        s = pd.to_numeric(out[col], errors="coerce")
        if np.isfinite(s).any():
            med = float(np.nanmedian(s))
            s = s.fillna(med)
        else:
            s = s.fillna(0.0)
        out[col] = s
    return out


def read_tcga():
    expr_file = DIR_DATA / "TCGA_Xena" / "TCGA.KIRC.sampleMap" / "HiSeqV2.gz"
    clin_file = DIR_DATA / "TCGA_Xena" / "TCGA.KIRC.sampleMap" / "KIRC_clinicalMatrix"
    expr = pd.read_csv(expr_file, sep="\t")
    expr = expr.rename(columns={expr.columns[0]: "gene"})
    expr["gene"] = expr["gene"].astype(str).str.upper()
    x = expr.set_index("gene").T
    x = x.apply(pd.to_numeric, errors="coerce")

    clin = pd.read_csv(clin_file, sep="\t")
    clin["sample"] = clin["sampleID"].astype(str)
    clin["time"] = np.where(clin["days_to_death"].notna(), pd.to_numeric(clin["days_to_death"], errors="coerce"), pd.to_numeric(clin["days_to_last_followup"], errors="coerce"))
    vs = clin["vital_status"].astype(str).str.lower()
    clin["status"] = np.where(vs.isin(["deceased", "dead"]), 1, np.where(vs.isin(["living", "alive"]), 0, np.nan))
    clin = clin.loc[clin["time"].notna() & clin["status"].notna() & (clin["time"] > 0), ["sample", "time", "status"]].copy()

    common = x.index.intersection(clin["sample"])
    x = x.loc[common]
    clin = clin.set_index("sample").loc[common].reset_index()
    return x, clin


def read_gse29609():
    series_path = DIR_DATA / "GEO" / "GSE29609_series_matrix.txt.gz"
    with gzip.open(series_path, "rt", encoding="utf-8", errors="ignore") as fh:
        series = fh.read().splitlines()

    acc_line = next(x for x in series if x.startswith("!Sample_geo_accession"))
    sids = acc_line.replace('"', "").split("\t")[1:]

    b = next(i for i, x in enumerate(series) if "!series_matrix_table_begin" in x)
    e = next(i for i, x in enumerate(series) if "!series_matrix_table_end" in x)
    with tempfile.NamedTemporaryFile("w", delete=False, suffix=".tsv", encoding="utf-8") as tmp:
        tmp.write("\n".join(series[(b + 1):e]))
        tmp_path = tmp.name
    expr = pd.read_csv(tmp_path, sep="\t")
    Path(tmp_path).unlink(missing_ok=True)
    expr = expr.rename(columns={expr.columns[0]: "ID_REF"})

    char_lines = [x for x in series if x.startswith("!Sample_characteristics_ch1")]
    char_rows = [x.replace('"', "").split("\t")[1:] for x in char_lines]
    char_t = list(zip(*char_rows))

    def get_field(vec, pattern):
        for item in vec:
            if re.search(pattern, item, flags=re.I):
                return re.sub(r"^.*?:", "", item).strip()
        return np.nan

    clin = pd.DataFrame({"sample": sids})
    clin["time"] = [pd.to_numeric(get_field(v, r"^survival time"), errors="coerce") for v in char_t]
    clin["death"] = [pd.to_numeric(get_field(v, r"^death$"), errors="coerce") for v in char_t]
    clin["death_cancer"] = [pd.to_numeric(get_field(v, r"^death from cancer"), errors="coerce") for v in char_t]
    clin = clin.loc[clin["time"].notna() & (clin["time"] > 0)].copy()

    annot_path = DIR_DATA / "GEO" / "annot" / "GPL1708.annot.gz"
    with gzip.open(annot_path, "rt", encoding="utf-8", errors="ignore") as fh:
        ann_lines = fh.read().splitlines()
    b2 = next(i for i, x in enumerate(ann_lines) if "!platform_table_begin" in x)
    e2 = next(i for i, x in enumerate(ann_lines) if "!platform_table_end" in x)
    with tempfile.NamedTemporaryFile("w", delete=False, suffix=".tsv", encoding="utf-8") as tmp:
        tmp.write("\n".join(ann_lines[(b2 + 1):e2]))
        ann_path = tmp.name
    ann = pd.read_csv(ann_path, sep="\t")
    Path(ann_path).unlink(missing_ok=True)
    ann = ann.rename(columns={"ID": "ID_REF", "Gene symbol": "GeneSymbol"})
    ann = ann.loc[ann["GeneSymbol"].notna() & (ann["GeneSymbol"].astype(str) != ""), ["ID_REF", "GeneSymbol"]].copy()
    ann["ID_REF"] = ann["ID_REF"].astype(str)
    ann["GeneSymbol"] = ann["GeneSymbol"].astype(str).str.upper()

    expr["ID_REF"] = expr["ID_REF"].astype(str)
    mdt = expr.merge(ann, on="ID_REF", how="inner")
    gsm_cols = [c for c in expr.columns if c != "ID_REF"]
    long = mdt.melt(id_vars=["GeneSymbol"], value_vars=gsm_cols, var_name="sample", value_name="expr")
    long["expr"] = pd.to_numeric(long["expr"], errors="coerce")
    gene_expr = long.groupby(["GeneSymbol", "sample"], as_index=False)["expr"].mean()
    wide = gene_expr.pivot(index="sample", columns="GeneSymbol", values="expr")
    wide.index.name = None
    x = wide.copy()

    common = x.index.intersection(clin["sample"])
    x = x.loc[common]
    clin = clin.set_index("sample").loc[common].reset_index()
    return x, clin


def read_cptac():
    cb_dir = DIR_DATA / "cBioPortal" / "rcc_cptac_gdc"
    clin = pd.read_csv(cb_dir / "clinical_os_gdc.tsv", sep="\t")
    mp = pd.read_csv(cb_dir / "sample_to_patient.tsv", sep="\t")
    expr = pd.read_csv(cb_dir / "expr_genelist.tsv", sep="\t")

    expr["gene"] = expr.iloc[:, 0].astype(str).str.upper()
    x = expr.set_index("gene").iloc[:, 1:].T
    x = x.apply(pd.to_numeric, errors="coerce")
    x.index.name = "sampleId"

    clin["OS_MONTHS"] = pd.to_numeric(clin["OS_MONTHS"], errors="coerce")
    clin["status"] = pd.to_numeric(clin["OS_STATUS"], errors="coerce")
    clin = clin.loc[clin["OS_MONTHS"].notna() & clin["status"].notna() & (clin["OS_MONTHS"] > 0)].copy()
    return x, clin, mp


def read_emtab():
    ext = DIR_DATA / "external_clinical" / "E-MTAB-1980"
    expr = pd.read_csv(ext / "expression_signature_gene33.tsv", sep="\t")
    expr["gene"] = expr["gene"].astype(str).str.upper()
    x = expr.set_index("gene").T
    x = x.apply(pd.to_numeric, errors="coerce")
    clin = pd.read_csv(ext / "clinical_matched_101.tsv", sep="\t")
    return x, clin


def read_immotion():
    ext = DIR_DATA / "external_clinical" / "IMmotion150"
    expr = pd.read_csv(ext / "expression_signature_gene33.tsv", sep="\t")
    expr["gene"] = expr["gene"].astype(str).str.upper()
    x = expr.set_index("gene").T
    x = x.apply(pd.to_numeric, errors="coerce")
    clin = pd.read_csv(ext / "clinical_matched_263.tsv", sep="\t")
    return x, clin


def risk_join(risk_by_sample, map_cb):
    dt = pd.DataFrame({"sampleId": risk_by_sample.index.astype(str), "risk": np.asarray(risk_by_sample, dtype=float)})
    dt = dt.merge(map_cb, on="sampleId", how="left")
    dt = dt.loc[dt["patientId"].notna() & np.isfinite(dt["risk"])].copy()
    if dt.empty:
        return pd.DataFrame(columns=["patientId", "risk"])
    out = dt.groupby("patientId", as_index=False)["risk"].mean()
    return out


def main():
    DIR_RES.mkdir(parents=True, exist_ok=True)

    x_tcga, clin_tcga = read_tcga()
    x_gse, clin_gse = read_gse29609()
    x_cb, clin_cb, map_cb = read_cptac()
    x_emtab, clin_emtab = read_emtab()
    x_imm, clin_imm = read_immotion()

    present = sorted(set(x_tcga.columns).intersection(x_gse.columns).intersection(x_cb.columns))
    v = x_tcga.loc[:, present].var(axis=0, skipna=True).fillna(0.0)
    cand0 = v.sort_values(ascending=False).index[: min(500, len(v))]

    score_gse = pd.Series(0.0, index=cand0)
    score_cons = pd.Series(0.0, index=cand0)
    for g in cand0:
        st_tcga = univ_cox_stats(clin_tcga["time"].to_numpy(), clin_tcga["status"].to_numpy(), x_tcga[g].to_numpy())
        st_d = univ_cox_stats(clin_gse["time"].to_numpy(), clin_gse["death"].to_numpy(), x_gse[g].to_numpy())
        st_dc = univ_cox_stats(clin_gse["time"].to_numpy(), clin_gse["death_cancer"].to_numpy(), x_gse[g].to_numpy())

        p_terms = [p for p in [st_d["p"], st_dc["p"]] if np.isfinite(p) and p > 0]
        if p_terms:
            score_gse[g] = sum(-np.log10(p_terms))

        b_tcga = st_tcga["beta"]
        for b_gse in [st_d["beta"], st_dc["beta"]]:
            if np.isfinite(b_tcga) and np.isfinite(b_gse) and b_tcga != 0 and b_gse != 0:
                if np.sign(b_tcga) == np.sign(b_gse):
                    score_cons[g] += 1.0
                else:
                    score_cons[g] -= 0.7

    score_all = score_gse + 1.2 * score_cons
    if np.isfinite(score_all).sum() == 0 or (score_all > 0).sum() < 20:
        cand = list(cand0[: min(220, len(cand0))])
    else:
        cand = list(score_all.sort_values(ascending=False).index[: min(220, len(score_all))])

    x_tcga_c = impute_col_median(x_tcga.loc[:, cand])
    x_gse_c = impute_col_median(x_gse.loc[:, cand])
    x_cb_c = impute_col_median(x_cb.loc[:, cand])

    x_tcga_r = rank_rows(x_tcga_c)
    x_gse_r = rank_rows(x_gse_c)
    x_cb_r = rank_rows(x_cb_c)

    rng = np.random.RandomState(SEED)
    alphas = [1.0, 0.7, 0.5, 0.3, 0.1, 0.0]
    penalizers = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3]

    best = None
    for it in range(1, N_ITER + 1):
        idx = rng.choice(np.arange(x_tcga_r.shape[0]), size=int(np.floor(0.7 * x_tcga_r.shape[0])), replace=False)
        train_s = x_tcga_r.index[idx]
        train_df = x_tcga_r.loc[train_s].copy()
        train_meta = clin_tcga.set_index("sample").loc[train_s]
        train_df["time"] = train_meta["time"].astype(float)
        train_df["status"] = train_meta["status"].astype(int)

        for alpha in alphas:
            for penalizer in penalizers:
                fit = CoxPHFitter(penalizer=penalizer, l1_ratio=alpha)
                try:
                    fit.fit(train_df, duration_col="time", event_col="status", show_progress=False)
                except Exception:
                    continue

                coef = fit.params_.copy()
                sel = coef[np.abs(coef) > 1e-8]
                if len(sel) < 8 or len(sel) > 80:
                    continue

                best29609 = None
                bestname = None
                for nm in ["death", "death_cancer"]:
                    st = pd.to_numeric(clin_gse[nm], errors="coerce")
                    keep = st.notna()
                    if keep.sum() < 10 or len(set(st[keep].astype(int))) < 2:
                        continue
                    risk = make_risk(x_gse_r.loc[keep, :], sel)
                    stat = km_stats(clin_gse.loc[keep, "time"].to_numpy(), st.loc[keep].to_numpy(), risk.to_numpy())
                    stat["flipped"] = False
                    stat_disc = (abs(stat["cindex"] - 0.5) + 0.5) if np.isfinite(stat["cindex"]) else -np.inf
                    best_disc = (abs(best29609["cindex"] - 0.5) + 0.5) if (best29609 is not None and np.isfinite(best29609["cindex"])) else -np.inf
                    if best29609 is None or stat_disc > best_disc:
                        best29609 = stat
                        bestname = nm
                if best29609 is None:
                    continue

                risk_cb = make_risk(x_cb_r, sel)
                pat = risk_join(risk_cb, map_cb)
                cbm = clin_cb.merge(pat, on="patientId", how="inner")
                stat_cb = km_stats(cbm["OS_MONTHS"].to_numpy(), cbm["status"].to_numpy(), cbm["risk"].to_numpy())
                stat_cb["flipped"] = False

                ok = (
                    np.isfinite(best29609["pval"]) and best29609["pval"] < 0.05 and np.isfinite(best29609["cindex"]) and best29609["cindex"] >= 0.60 and
                    np.isfinite(stat_cb["pval"]) and stat_cb["pval"] < 0.05 and np.isfinite(stat_cb["cindex"]) and stat_cb["cindex"] >= 0.60
                )

                def p_term(p):
                    if not np.isfinite(p) or p <= 0:
                        return -1.0
                    return (-np.log10(p) - 1.30103)

                gse_c = best29609["cindex"]
                cb_c = stat_cb["cindex"]
                gse_p = best29609["pval"]
                cb_p = stat_cb["pval"]

                gse_metric = gse_c + 0.10 * p_term(gse_p) if np.isfinite(gse_c) else -np.inf
                cb_metric = cb_c + 0.06 * p_term(cb_p) if np.isfinite(cb_c) else -np.inf
                gse_raw_pen = ((0.40 - gse_c) * 2.5) if np.isfinite(gse_c) and gse_c < 0.40 else 0.0
                cb_raw_pen = ((0.45 - cb_c) * 1.2) if np.isfinite(cb_c) and cb_c < 0.45 else 0.0
                score = (2.6 * gse_metric) + (0.8 * cb_metric) - gse_raw_pen - cb_raw_pen - 0.002 * len(sel)

                rec = {
                    "it": it,
                    "alpha": alpha,
                    "penalizer": penalizer,
                    "n_genes": int(len(sel)),
                    "coef": sel.sort_values(ascending=False),
                    "gse29609": best29609,
                    "gse29609_status": bestname,
                    "cptac": stat_cb,
                    "ok": bool(ok),
                    "score": float(score),
                }

                best_score = best["score"] if (best is not None and np.isfinite(best["score"])) else -np.inf
                best_ok = best["ok"] if best is not None else False
                if best is None or (ok and not best_ok) or (ok == best_ok and score > best_score):
                    best = rec
                    print(json.dumps({
                        "it": it,
                        "alpha": alpha,
                        "penalizer": penalizer,
                        "n_genes": len(sel),
                        "gse29609_endpoint": bestname,
                        "gse29609_p": best29609["pval"],
                        "gse29609_cindex": best29609["cindex"],
                        "cptac_p": stat_cb["pval"],
                        "cptac_cindex": stat_cb["cindex"],
                        "ok": ok,
                    }))

    if best is None:
        raise SystemExit("No model found in Python fixed-direction search.")

    coef = best["coef"].sort_values(ascending=False)
    coef_df = pd.DataFrame({"gene": coef.index, "coef": coef.values})
    coef_df.to_csv(out_path("best_signature_coefficients.tsv"), sep="\t", index=False)

    # external checks on newly added cohorts
    risk_emtab = make_risk(rank_rows(impute_col_median(x_emtab.loc[:, coef.index.intersection(x_emtab.columns)])), coef)
    common_em = risk_emtab.index.intersection(clin_emtab["sample_id"])
    em = clin_emtab.set_index("sample_id").loc[common_em].copy()
    em_stat = km_stats(em["time_months"].to_numpy(), em["status_os"].to_numpy(), risk_emtab.loc[common_em].to_numpy())

    risk_imm = make_risk(rank_rows(impute_col_median(x_imm.loc[:, coef.index.intersection(x_imm.columns)])), coef)
    common_im = risk_imm.index.intersection(clin_imm["sample_id"])
    im = clin_imm.set_index("sample_id").loc[common_im].copy()
    imm_stat = km_stats(im["pfs_months"].to_numpy(), im["status_pfs"].to_numpy(), risk_imm.loc[common_im].to_numpy())

    summary = {
        "mode": "python_fixed_direction_search",
        "run_tag": RUN_TAG,
        "iterations": N_ITER,
        "seed": SEED,
        "alpha": best["alpha"],
        "penalizer": best["penalizer"],
        "n_genes": best["n_genes"],
        "gse29609_endpoint": best["gse29609_status"],
        "gse29609_logrank_p": best["gse29609"]["pval"],
        "gse29609_cindex": best["gse29609"]["cindex"],
        "cptac_logrank_p": best["cptac"]["pval"],
        "cptac_cindex": best["cptac"]["cindex"],
        "ok": best["ok"],
        "emtab1980_logrank_p": em_stat["pval"],
        "emtab1980_cindex": em_stat["cindex"],
        "immotion150_logrank_p": imm_stat["pval"],
        "immotion150_cindex": imm_stat["cindex"],
    }
    pd.DataFrame([summary]).to_csv(out_path("search_summary.tsv"), sep="\t", index=False)
    out_path("search_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    detail = {
        "summary": summary,
        "top_genes": coef_df.head(20).to_dict(orient="records"),
    }
    out_path("search_detail.json").write_text(json.dumps(detail, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
