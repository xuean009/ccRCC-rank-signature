# Cross-platform rank-based prognostic signature for clear cell renal cell carcinoma

This repository contains the cleaned core code used to build and visualize a rank-based prognostic signature for clear cell renal cell carcinoma (ccRCC), including bulk transcriptome modeling, external validation, single-cell RNA-seq visualization, and spatial transcriptomics visualization.

## Abstract

Transcriptome-derived prognostic signatures for clear cell renal cell carcinoma (ccRCC) often lose performance after transfer across platforms and cohorts. We developed a 33-gene rank-based Cox signature in TCGA-KIRC and examined its transportability in GSE29609 and CPTAC using only public data. Shared genes were prefiltered by variance and inter-cohort relevance, expression values were converted to within-sample ranks during model development, and the finalized coefficient set was projected to external cohorts. In CPTAC, the original score direction showed significant overall-survival stratification with a C-index of 0.722 (95% bootstrap confidence interval [CI] 0.652-0.791; log-rank p<1e-4). In GSE29609, however, the original score direction produced inverse concordance (C-index 0.262, p=0.060), whereas a sign-reversed sensitivity analysis yielded a C-index of 0.738 (95% CI 0.599-0.849; p=0.021), indicating preserved ranking but unstable score direction. In TCGA, the signature remained independently associated with survival beyond age, pathologic stage, and grade, improving the clinicopathologic-model C-index from 0.762 to 0.812. In GSE29609, an adapted SSIGN-like clinicopathologic score based on available T/N/M stage, Fuhrman grade, and necrosis achieved a C-index of 0.794, and adding the sign-reversed signature produced only limited incremental improvement (C-index 0.803; likelihood-ratio p=0.275). These results suggest that the model captures biologically relevant prognostic structure, but score direction is not fully transportable across cohorts. The signature should therefore be interpreted as a reproducible, hypothesis-generating biomarker that requires external calibration and prospective validation before clinical deployment.

## Repository layout

- `scripts/01_download_tcga_xena.R`: download TCGA-KIRC clinical and expression matrices from UCSC Xena.
- `scripts/02_download_cptac_cbioportal.py`: download CPTAC sample-to-patient mapping and cBioPortal OS metadata.
- `scripts/03_build_common_gene_list.R`: build the shared high-variance candidate gene list across TCGA and GSE29609.
- `scripts/04_fetch_cptac_expression.py`: fetch CPTAC expression for the shared candidate genes through cBioPortal.
- `scripts/05_download_cptac_gdc_clinical.py`: retrieve complete CPTAC follow-up from the GDC API.
- `scripts/06_build_rank_signature.R`: preprocess bulk data, perform within-sample rank transformation, fit the elastic-net Cox model, and validate it in GSE29609 and CPTAC.
- `scripts/07_plot_scrna_signature.R`: generate the core scRNA-seq visualizations for GSE159115.
- `scripts/08_plot_spatial_signature.R`: generate the core spatial transcriptomics visualizations for GSE250163.
- `scripts/_figure_base.py`: helper code for manuscript figure/table assembly.
- `scripts/09_build_heliyon_revision_package.py`: build the final Heliyon revision figure/table package.
- `results/best_signature_coefficients.tsv`: frozen final 33-gene coefficient table used for downstream visualization.

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

Bulk modeling and external validation:

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
2. Run `scripts/01_download_tcga_xena.R`, `scripts/02_download_cptac_cbioportal.py`, `scripts/03_build_common_gene_list.R`, `scripts/04_fetch_cptac_expression.py`, and `scripts/05_download_cptac_gdc_clinical.py` as needed.
3. Run `scripts/06_build_rank_signature.R` to train the rank-based elastic-net Cox model and write `results/best_signature_coefficients.tsv`.
4. Run `scripts/07_plot_scrna_signature.R` and `scripts/08_plot_spatial_signature.R` for the single-cell and spatial plots.
5. Run `scripts/09_build_heliyon_revision_package.py` to assemble the final manuscript figure/table package under `output/heliyon_revision_package/`.

## Notes

- All code paths have been converted to project-relative paths.
- Large raw data files, temporary outputs, and intermediate analysis products are intentionally excluded from version control.
- `results/best_signature_coefficients.tsv` is included because it is a small derived artifact required by the downstream visualization scripts and the final manuscript package.
