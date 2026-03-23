#!/usr/bin/env python3
"""Download complete OS follow-up for CPTAC RCC from GDC.

Goal: fix the issue that cBioPortal rcc_cptac_gdc exposes OS_MONTHS only for deceased patients.
We query GDC /cases by submitter_id (C3L-xxxx) to obtain days_to_death OR days_to_last_follow_up
plus vital_status, for BOTH deceased and living.

Outputs:
  ``data/cBioPortal/rcc_cptac_gdc/clinical_os_gdc.tsv``

Input:
  ``data/cBioPortal/rcc_cptac_gdc/sample_to_patient.tsv``
"""

import csv
import math
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
CB_DIR = ROOT / "data" / "cBioPortal" / "rcc_cptac_gdc"
MAP_PATH = CB_DIR / "sample_to_patient.tsv"
OUT_PATH = CB_DIR / "clinical_os_gdc.tsv"

GDC = "https://api.gdc.cancer.gov"
UA = "Mozilla/5.0"

# load patient ids
pids = set()
with MAP_PATH.open("r", encoding="utf-8") as f:
    for i, line in enumerate(f):
        if i == 0:
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2 and parts[1].strip():
            pids.add(parts[1].strip())

pids = sorted(pids)
print("patientIds", len(pids))

sess = requests.Session()
headers = {"Content-Type": "application/json", "Accept": "application/json", "User-Agent": UA}

FIELDS = [
    "submitter_id",
    # vital status is under demographic for CPTAC-3
    "demographic.vital_status",
    "demographic.days_to_death",
    # follow-up time may appear under diagnoses
    "diagnoses.days_to_last_follow_up",
]

# chunk query
chunk = 100
rows = []
for i in range(0, len(pids), chunk):
    sub = pids[i:i+chunk]
    flt = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "submitter_id", "value": sub}},
        ],
    }
    payload = {
        "filters": flt,
        "fields": ",".join(FIELDS),
        "format": "JSON",
        "size": len(sub),
    }
    r = sess.post(f"{GDC}/cases", headers=headers, json=payload, timeout=120)
    r.raise_for_status()
    hits = r.json().get("data", {}).get("hits", [])

    for h in hits:
        sid = h.get("submitter_id")
        diags = h.get("diagnoses") or []
        d0 = diags[0] if diags else {}
        demo = h.get("demographic") or {}
        vital = (demo.get("vital_status") or "").lower()
        dtd = demo.get("days_to_death")
        dtl = d0.get("days_to_last_follow_up")
        # choose time
        t = None
        if dtd is not None:
            try:
                t = float(dtd)
            except Exception:
                t = None
        if (t is None or not math.isfinite(t)) and dtl is not None:
            try:
                t = float(dtl)
            except Exception:
                t = None
        status = None
        if vital in ("dead", "deceased"):
            status = 1
        elif vital in ("alive", "living"):
            status = 0
        if sid and t is not None and math.isfinite(t) and t > 0 and status is not None:
            rows.append((sid, t / 30.44, status))

    time.sleep(0.15)

# write
rows.sort(key=lambda x: x[0])
with OUT_PATH.open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["patientId", "OS_MONTHS", "OS_STATUS"])
    for sid, months, status in rows:
        w.writerow([sid, f"{months:.6f}", status])

print("wrote", OUT_PATH, "rows", len(rows))
