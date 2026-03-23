#!/usr/bin/env python3
"""Download CPTAC ccRCC sample metadata from cBioPortal.

Outputs under ``data/cBioPortal/rcc_cptac_gdc/``:
  - ``sample_to_patient.tsv``
  - ``clinical_os.tsv`` (cBioPortal-reported OS, used only as a fallback)

The final model uses ``clinical_os_gdc.tsv`` from ``05_download_cptac_gdc_clinical.py``
for survival time because cBioPortal does not expose complete follow-up for living patients.
"""

import csv
import time
from pathlib import Path

import requests

BASE = "https://www.cbioportal.org/api"
UA = "Mozilla/5.0"
HEADERS_JSON = {"Accept": "application/json", "User-Agent": UA}

STUDY = "rcc_cptac_gdc"
ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "data" / "cBioPortal" / STUDY
OUTDIR.mkdir(parents=True, exist_ok=True)

s = requests.Session()

# samples & mapping
samples = s.get(f"{BASE}/studies/{STUDY}/samples", headers=HEADERS_JSON, timeout=60).json()
sample_ids = [x["sampleId"] for x in samples]
sample_to_patient = {x["sampleId"]: x.get("patientId") for x in samples}

map_path = OUTDIR / "sample_to_patient.tsv"
with map_path.open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["sampleId", "patientId"])
    for sid in sample_ids:
        w.writerow([sid, sample_to_patient.get(sid, "")])

# patient-level clinical (OS)
os_months = {}
os_status = {}
page = 0
page_size = 5000
while True:
    r = s.get(
        f"{BASE}/studies/{STUDY}/clinical-data",
        headers=HEADERS_JSON,
        params={"projection": "SUMMARY", "clinicalDataType": "PATIENT", "pageSize": page_size, "pageNumber": page},
        timeout=180,
    )
    r.raise_for_status()
    arr = r.json()
    if not arr:
        break
    for rec in arr:
        pid = rec.get("patientId")
        attr = rec.get("clinicalAttributeId")
        val = rec.get("value")
        if attr == "OS_MONTHS":
            os_months[pid] = val
        elif attr == "OS_STATUS":
            os_status[pid] = val
    page += 1
    time.sleep(0.15)

keep_patients = sorted(set(os_months) & set(os_status))

clin_path = OUTDIR / "clinical_os.tsv"
with clin_path.open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["patientId", "OS_MONTHS", "OS_STATUS"])
    for pid in keep_patients:
        w.writerow([pid, os_months.get(pid, ""), os_status.get(pid, "")])

print("OK", STUDY, "patients_os", len(keep_patients), "samples", len(sample_ids))
