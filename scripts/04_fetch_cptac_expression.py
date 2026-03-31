#!/usr/bin/env python3
"""Fetch CPTAC expression for the shared candidate gene set.

Input: ``data/gene_lists/common_topvar_1500.tsv`` (column: ``gene``)
Output: ``data/cBioPortal/rcc_cptac_gdc/expr_genelist.tsv`` (gene x sample)

This expands CPTAC expression beyond the 10-gene signature so optimization can use a large shared gene set.
"""

import csv
import json
import time
from pathlib import Path

import requests

BASE = "https://www.cbioportal.org/api"
UA = "Mozilla/5.0"
H_JSON = {"Accept": "application/json", "User-Agent": UA}
H_POST = {"Accept": "application/json", "Content-Type": "application/json", "User-Agent": UA}

STUDY = "rcc_cptac_gdc"
ROOT = Path(__file__).resolve().parents[1]
GENELIST = ROOT / "data" / "gene_lists" / "common_topvar_1500.tsv"
OUTDIR = ROOT / "data" / "cBioPortal" / STUDY
OUTDIR.mkdir(parents=True, exist_ok=True)
OUT_EXPR = OUTDIR / "expr_genelist.tsv"

sess = requests.Session()

# expression profile: TPM
profiles_resp = sess.get(f"{BASE}/studies/{STUDY}/molecular-profiles", headers=H_JSON, timeout=60)
profiles_resp.raise_for_status()
profiles = profiles_resp.json()
expr_profile = None
for p in profiles:
    pid = p.get("molecularProfileId", "")
    if pid.endswith("mrna_seq_tpm") and p.get("molecularAlterationType") == "MRNA_EXPRESSION":
        expr_profile = pid
        break
if expr_profile is None:
    raise SystemExit("No TPM profile found")

# samples (use all)
samples_resp = sess.get(f"{BASE}/studies/{STUDY}/samples", headers=H_JSON, timeout=60)
samples_resp.raise_for_status()
samples = samples_resp.json()
sample_ids = [x["sampleId"] for x in samples]

# load gene list
with GENELIST.open("r", encoding="utf-8") as f:
    rows = [r.strip().split("\t")[0] for r in f.read().splitlines() if r.strip()]
# drop header
genes = [g.strip().upper() for g in rows if g.strip().upper() != "GENE"]

# build hugo->entrez map once by downloading genes table (paged)
# This avoids per-gene keyword lookups.
hugo_to_entrez = {}
page = 0
page_size = 20000
r = sess.get(f"{BASE}/genes", headers=H_JSON, params={"pageSize": page_size, "pageNumber": 0}, timeout=180)
r.raise_for_status()
arr = r.json()
for rec in arr:
    hugo = str(rec.get("hugoGeneSymbol", "")).upper()
    eid = rec.get("entrezGeneId")
    if hugo and eid is not None:
        hugo_to_entrez[hugo] = int(eid)

entrez = []
kept_genes = []
for g in genes:
    if g in hugo_to_entrez:
        kept_genes.append(g)
        entrez.append(hugo_to_entrez[g])

print(f"genes requested={len(genes)} mapped={len(kept_genes)} samples={len(sample_ids)}")

# init matrix as dict of gene->list of values aligned to sample_ids
mat = {g: [""] * len(sample_ids) for g in kept_genes}

# fetch by chunks of entrez ids to keep payload reasonable
chunk_genes = 200
for i in range(0, len(entrez), chunk_genes):
    eids = entrez[i:i+chunk_genes]
    gsub = kept_genes[i:i+chunk_genes]
    payload = {"sampleIds": sample_ids, "entrezGeneIds": eids}
    resp = sess.post(f"{BASE}/molecular-profiles/{expr_profile}/molecular-data/fetch", headers=H_POST, data=json.dumps(payload), timeout=240)
    if resp.status_code != 200:
        raise SystemExit(f"fetch failed {resp.status_code}: {resp.text[:200]}")
    data = resp.json()
    # index sample
    sidx = {sid: j for j, sid in enumerate(sample_ids)}
    # map eid->gene
    eid2g = {eid: g for eid, g in zip(eids, gsub)}
    for rec in data:
        eid = rec.get("entrezGeneId")
        sid = rec.get("sampleId")
        val = rec.get("value")
        g = eid2g.get(eid)
        j = sidx.get(sid)
        if g is not None and j is not None:
            mat[g][j] = "" if val is None else str(val)
    time.sleep(0.15)

with OUT_EXPR.open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gene"] + sample_ids)
    for g in kept_genes:
        w.writerow([g] + mat[g])

print("Wrote", OUT_EXPR)
