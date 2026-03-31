# Cross-platform rank-based prognostic signature for clear cell renal cell carcinoma

This repository contains the cleaned core code used to build and visualize a rank-based prognostic signature for clear cell renal cell carcinoma (ccRCC), including bulk transcriptome modeling, public-cohort frozen-score projection, single-cell RNA-seq visualization, and spatial transcriptomics visualization.

## Public summary

This repository packages the research code behind a kidney-cancer prognosis study. In simple terms, the project asks whether a gene-expression risk score for clear cell renal cell carcinoma can still carry useful prognostic information after it is moved across different public datasets, platforms, and clinical settings.

Instead of relying on absolute expression values, the retained model uses within-sample gene ranks, which makes the score less sensitive to platform-specific scaling differences. The current public-facing workflow starts from a frozen 33-gene coefficient table and projects that locked score into multiple public cohorts rather than retraining the model in each dataset.

For most readers, the main takeaways are:

- this is a reproducibility-oriented research repository, not a clinical software package
- the strongest independent overall-survival support in the current analysis comes from `E-MTAB-1980`
- `CPTAC`, `IMmotion150`, and `E-MTAB-3267` are included to test transportability across additional public settings, but they play different inferential roles
- the current manuscript-oriented build entry point is `scripts/17_build_cancer_medicine_package.py`

## Abstract

Transcriptome-derived prognostic signatures for clear cell renal cell carcinoma (ccRCC) often lose performance after transfer across platforms and cohorts. The current Cancer Medicine revision is framed around a legacy frozen 33-gene rank-based score projected across four public cohorts spanning RNA-seq and microarray platforms and both resected and therapy-specific settings: CPTAC (n=237, overall survival), E-MTAB-1980 (n=101, overall survival), IMmotion150 (n=263, progression-free survival), and E-MTAB-3267 (n=53, progression-free survival). Archived search-stage scripts indicate that GSE29609 informed early feature filtering and CPTAC contributed to historical cross-platform model-definition constraints, so E-MTAB-1980 is treated as the primary independent overall-survival validation cohort, CPTAC as a development-linked overall-survival projection cohort, and IMmotion150 plus E-MTAB-3267 as secondary endpoint/context sensitivity cohorts. The original score direction remains associated with outcome in E-MTAB-1980 (C-index 0.746, 95% CI 0.650-0.838) and CPTAC (C-index 0.656, 95% CI 0.549-0.749), with additional directionally concordant support in IMmotion150 and E-MTAB-3267. In TCGA-KIRC, the score adds prognostic information beyond age, pathologic stage, and grade, improving the clinicopathologic-model C-index from 0.762 to 0.814 and the 1000-bootstrap optimism-corrected C-index to 0.810; in an exploratory E-MTAB-1980 benchmark, adding the score increases the clinicopathologic-model C-index from 0.804 to 0.824. These findings support cross-platform reproducibility while underscoring the need for prospective validation and external calibration before clinical deployment.

## Repository layout

- `scripts/01_download_tcga_xena.R`: download TCGA-KIRC clinical and expression matrices from UCSC Xena.
- `scripts/02_download_cptac_cbioportal.py`: download CPTAC sample-to-patient mapping and cBioPortal OS metadata.
- `scripts/03_build_common_gene_list.R`: build the shared high-variance candidate gene list across TCGA and GSE29609.
- `scripts/04_fetch_cptac_expression.py`: fetch CPTAC expression for the shared candidate genes through cBioPortal.
- `scripts/05_download_cptac_gdc_clinical.py`: retrieve complete CPTAC follow-up from the GDC API.
- `scripts/06_build_rank_signature.R`: archival exploratory discovery script retained for traceability. It reflects the earlier search-stage optimizer in which GSE29609 and CPTAC influenced model-definition constraints and should not be treated as the current Cancer Medicine projection pipeline.
- `scripts/07_plot_scrna_signature.R`: generate the core scRNA-seq visualizations for GSE159115.
- `scripts/08_plot_spatial_signature.R`: generate the core spatial transcriptomics visualizations for GSE250163.
- `scripts/_figure_base.py`: helper code for manuscript figure/table assembly.
- `scripts/09_build_heliyon_revision_package.py`: build the archived Heliyon revision package.
- `scripts/10_validate_emtab1980.py`: score the frozen signature in E-MTAB-1980.
- `scripts/11_validate_immotion150.py`: score the frozen signature in IMmotion150.
- `scripts/14_validate_emtab3267.R`: score the frozen signature in E-MTAB-3267.
- `scripts/15_build_frontiers_submission_package.py`: archived Frontiers package builder retained only for provenance.
- `scripts/16_check_frontiers_consistency.py`: archived Frontiers consistency checker retained only for provenance.
- `scripts/17_build_cancer_medicine_package.py`: build the current Cancer Medicine submission package aligned to the study hierarchy defined for this analysis.
- `results/best_signature_coefficients.tsv`: frozen legacy 33-gene coefficient table used for downstream projection and visualization.

## Tested environment

The original analysis/build logs captured the following software stack:

- R 4.3.3
- Bioconductor 3.18
- Seurat 5.4.0
- SeuratObject 5.3.0
- glmnet 4.1-10
- data.table 1.18.2.1
- ggplot2 4.0.2
- patchwork 1.3.2
- cowplot 1.2.0
- UCSCXenaTools 1.6.1
- Python 3.10.9
- pandas 2.2.3
- numpy 1.26.4
- matplotlib 3.7.0

Additional Python packages required by the manuscript figure builder:

- scipy
- statsmodels
- requests

## Public data sources

Bulk modeling and prespecified Cancer Medicine frozen-score projection:

- TCGA-KIRC from UCSC Xena
  - Clinical information and IlluminaHiSeq RNASeqV2 expression
  - Downloaded by `scripts/01_download_tcga_xena.R`
  - Portal: <https://xena.ucsc.edu/>
- GEO GSE29609
  - Series matrix: <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE29609/matrix/GSE29609_series_matrix.txt.gz>
  - Platform annotation GPL1708: <https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL1nnn/GPL1708/annot/GPL1708.annot.gz>
- CPTAC ccRCC cohort
  - cBioPortal study: `rcc_cptac_gdc`
  - cBioPortal: <https://www.cbioportal.org/study/summary?id=rcc_cptac_gdc>
  - GDC API for complete follow-up: <https://api.gdc.cancer.gov/>

Single-cell and spatial analyses:

- GEO GSE159115 (scRNA-seq)
  - Supplementary raw archive: <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159115/suppl/GSE159115_RAW.tar>
  - Annotation file: <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159115/suppl/GSE159115_ccRCC_anno.csv.gz>
- GEO GSE250163 (spatial transcriptomics)
  - Supplementary raw archive: <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250163/suppl/GSE250163_RAW.tar>

Raw public datasets are not included in this repository. See [data/README.md](data/README.md) for expected file locations.

## Minimal run order

1. Download or place public bulk data files under `data/` as described in `data/README.md`.
2. Treat `results/best_signature_coefficients.tsv` as the frozen coefficient table used in the current Cancer Medicine manuscript.
3. Run `scripts/10_validate_emtab1980.py`, `scripts/11_validate_immotion150.py`, and `scripts/14_validate_emtab3267.R` to refresh the retained projection/validation summaries used in the current Cancer Medicine manuscript if needed.
4. Run `scripts/17_build_cancer_medicine_package.py` to assemble the current Cancer Medicine submission package under `output/cancer_medicine_submission_package/`.
5. Use `scripts/09_build_heliyon_revision_package.py`, `scripts/15_build_frontiers_submission_package.py`, and `scripts/16_check_frontiers_consistency.py` only when you intentionally need the archived journal-specific packages retained for provenance.

## Notes

- All code paths have been converted to project-relative paths.
- Large raw data files, temporary outputs, and intermediate analysis products are intentionally excluded from version control.
- `r_lib/` and `.mplconfig/` are optional local runtime caches and are intentionally excluded from version control; the scripts still fall back to the system R/Python environments when those directories are absent.
- `results/best_signature_coefficients.tsv` is included because it is a small derived artifact required by the downstream validation scripts and the final manuscript package.
- Ancillary exploratory files for GSE29609 and Choueiri2016 remain in the repository for provenance, but they are not part of the retained external-validation hierarchy described in the current Cancer Medicine analysis.
