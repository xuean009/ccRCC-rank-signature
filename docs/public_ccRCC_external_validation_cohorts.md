# Public ccRCC cohorts with survival endpoints for external validation

Date checked: 2026-03-25

This note summarizes public cohorts that are plausible additions to the current external-validation design beyond `GSE29609` and `CPTAC`.

## Recommended priority

### Tier 1: strongest additions for prognosis-oriented external validation

1. `GSE73731` (GEO)
   - Source: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73731>
   - Expression platform: Affymetrix HG-U133 Plus 2.0 (`GPL570`)
   - Tumor sample count: 265 ccRCC tumors
   - Why it is attractive:
     - Much larger than `GSE29609`
     - Pure ccRCC cohort
     - The associated publication explicitly reports cancer-specific-survival stratification in the 265-sample cohort: <https://pubmed.ncbi.nlm.nih.gov/28779136/>
   - Caveat:
     - Sample-level survival fields were not obvious in the GEO series-matrix header during spot checking, so survival time/status may need to be reconstructed from the linked publication or supplementary clinical table rather than the plain GEO matrix alone.
   - Practical recommendation:
     - Best candidate to replace `GSE29609` as the main microarray-based external validation set if the clinical supplement can be recovered cleanly.

2. `E-MTAB-1980` (ArrayExpress / BioStudies)
   - Study page: <https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-1980>
   - IDF file: <https://www.ebi.ac.uk/biostudies/files/E-MTAB-1980/E-MTAB-1980.idf.txt>
   - Processed expression file: <https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-1980/E-MTAB-1980.processed.1.zip>
   - Expression platform: Agilent (`A-MEXP-2183`)
   - Tumor sample count: 101 ccRCC samples
   - Why it is attractive:
     - Independent ccRCC cohort
     - Public processed expression matrix is directly downloadable
     - ArrayExpress metadata links it to the UTokyo Nat Genet 2013 ccRCC molecular study: <https://pubmed.ncbi.nlm.nih.gov/23797736/>
   - Caveat:
     - The ArrayExpress `sdrf` exposed obvious sex/disease/sample identifiers, but not survival columns in the quick metadata check. In many later ccRCC prognostic studies this cohort is used as an OS validation cohort, so the survival table likely exists in linked clinical supplements rather than in the simple SDRF itself.
   - Practical recommendation:
     - Strong OS-oriented external-validation candidate, but use only after the clinical follow-up table is pinned down and archived locally.

### Tier 2: usable, but endpoint or clinical context differs from classic untreated OS validation

3. `E-MTAB-3267` (ArrayExpress / BioStudies)
   - Study page: <https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-3267>
   - SDRF file: <https://www.ebi.ac.uk/biostudies/files/E-MTAB-3267/E-MTAB-3267.sdrf.txt>
   - Expression platform: Affymetrix HuGene (`A-AFFY-141`)
   - Sample count in the associated publication: 53
   - Endpoint available from SDRF:
     - `progression free survival`
     - `progression`
     - treatment-response annotations under sunitinib
   - Supporting publication: <https://pubmed.ncbi.nlm.nih.gov/28779136/>
   - Caveat:
     - This is a treated metastatic setting and the endpoint is PFS rather than OS.
   - Practical recommendation:
     - Good as a secondary transportability/sensitivity cohort, not ideal as the only external survival validation set for a general ccRCC prognosis signature.

4. `ccrcc_iatlas_choueiri_2016` (cBioPortal / CRI iAtlas harmonized)
   - Study page: <https://www.cbioportal.org/study/summary?id=ccrcc_iatlas_choueiri_2016>
   - cBioPortal API study id: `ccrcc_iatlas_choueiri_2016`
   - RNA-seq samples: 16
   - Endpoints exposed by cBioPortal API:
     - `OS_MONTHS`, `OS_STATUS`
     - `PFS_MONTHS`, `PFS_STATUS`
   - Clinical context:
     - Metastatic clear-cell RCC, nivolumab phase 1 biomarker trial
   - Supporting study metadata:
     - cBioPortal study summary/API
     - citation PMID: <https://pubmed.ncbi.nlm.nih.gov/27169994/>
   - Practical recommendation:
     - Can be merged into a "small public metastatic ccRCC RNA-seq validation panel", but sample size is very small and therapy confounding is substantial.

5. `rcc_iatlas_immotion150_2018` (cBioPortal / CRI iAtlas harmonized)
   - Study page: <https://www.cbioportal.org/study/summary?id=rcc_iatlas_immotion150_2018>
   - cBioPortal API study id: `rcc_iatlas_immotion150_2018`
   - RNA-seq samples: 263
   - Endpoints exposed by cBioPortal API:
     - `PFS_MONTHS`, `PFS_STATUS`
   - Clinical context:
     - Treatment-naive metastatic RCC trial cohort
   - Caveat:
     - The study is labeled RCC rather than explicitly ccRCC in the portal title, so histology filtering should be checked before use as a strict ccRCC external set.
   - Practical recommendation:
     - Valuable as a large therapy-response/PFS sensitivity cohort after confirming histology labels; not a direct substitute for untreated localized-disease OS validation.
   - Current project status:
     - Successfully recovered and validated as a public RNA-seq/PFS sensitivity cohort.
     - Current frozen-signature result in this repository:
       - `n = 263`
       - events `= 164`
       - `log-rank P = 0.0003`
       - `C-index = 0.608`

## Legacy / access-challenging option

6. `ICGC RECA-EU`
   - Historically used in many ccRCC signature papers as an RNA-seq + clinical validation cohort, often around 91 patients after filtering.
   - What is currently recoverable from the official post-portal-retirement route:
     - The public `icgc25k-open` bucket still exposes `RECA-EU` donor/sample/specimen shards.
     - Official open-data headers confirm that donor rows include `donor_vital_status`, `donor_survival_time`, and `donor_interval_of_last_followup`.
     - Specimen rows include pathology fields such as histology, grade, and stage, so a ccRCC-like subset may still be filterable.
   - Caveat:
     - During this pass, the current public route did not expose a directly downloadable `exp_seq.RECA-EU` matrix or donor-level `exp_seq` shards.
     - The official mapping-TSV link returned `AccessDenied`, so the current alternate host for the expression files could not be resolved from public metadata alone.
   - Practical recommendation:
     - `RECA-EU` remains biologically attractive, but do not count it as a reproducible public RNA-seq validation cohort unless the expression download path can be pinned down cleanly.

## What this means for the current project

### Best near-term upgrade path

1. Keep `CPTAC` as the proteogenomic external cohort.
2. Add `GSE73731` as the main microarray external validation cohort.
3. Add `E-MTAB-1980` if the linked survival table can be recovered cleanly.

This would give a much stronger validation story than `GSE29609` alone because:

- `GSE73731` is substantially larger than `GSE29609`
- `E-MTAB-1980` is an independent non-TCGA public ccRCC cohort often reused in prognostic studies
- both remain closer to a general prognosis setting than therapy-trial PFS datasets

### If survival tables cannot be cleanly recovered for Tier 1

Fallback plan:

1. Treat `E-MTAB-3267` and `rcc_iatlas_immotion150_2018` as endpoint-shifted validation cohorts.
2. Present them explicitly as PFS-based transportability analyses, not as direct OS validation.
3. Keep `GSE29609` only as a small historical sensitivity cohort.

## Pooling strategy if multiple small public cohorts are merged

Do not directly pool all cohorts into one Cox model with a single endpoint unless they share:

- the same endpoint definition
- comparable treatment setting
- compatible censoring conventions

Safer alternatives:

1. Rank-normalize within each cohort, compute the frozen signature score, then meta-analyze cohort-specific HRs/C-indices/log-rank P values.
2. Pool only OS-like cohorts together (`GSE29609`, `GSE73731`, `E-MTAB-1980`, possibly `ICGC RECA-EU` if accessible) and stratify Cox models by cohort.
3. Keep PFS cohorts (`E-MTAB-3267`, `IMmotion150`, `Choueiri 2016`) in a separate sensitivity-validation section.

## Bottom line

The most realistic public replacements or upgrades for `GSE29609` are:

- `GSE73731`
- `E-MTAB-1980`

The most realistic mergeable small public survival/PFS cohorts are:

- `E-MTAB-3267`
- `ccrcc_iatlas_choueiri_2016`

The largest additional public therapy-trial cohort worth checking is:

- `rcc_iatlas_immotion150_2018`
