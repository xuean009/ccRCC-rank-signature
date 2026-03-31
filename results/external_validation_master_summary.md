# External Validation Master Summary

This summary now follows the prespecified Cancer Medicine study hierarchy and lists only the four public cohorts retained in the cautious projection/validation framework.

| Cohort | Endpoint | Platform | n | Events | Log-rank P | C-index | Flipped | Suggested role |
| --- | --- | --- | ---: | ---: | ---: | ---: | :---: | --- |
| CPTAC | overall survival | RNA-seq | 237 | 36 | 0.016 | 0.656 | No | Development-linked OS projection cohort |
| E-MTAB-1980 | overall survival | microarray (Agilent A-MEXP-2183) | 101 | 23 | 0.0007 | 0.746 | No | Primary independent OS validation cohort |
| IMmotion150 | progression-free survival | RNA-seq (TPM) | 263 | 164 | 0.0003 | 0.608 | No | Secondary endpoint/context sensitivity cohort |
| E-MTAB-3267 | progression-free survival | microarray (Affymetrix Human Gene 1.0 ST) | 53 | 39 | 0.080 | 0.644 | No | Secondary endpoint/context sensitivity cohort |

## Notes

- `GSE29609` and `Choueiri2016` remain in the repository only as ancillary exploratory artifacts and are not part of the prespecified validation/projection hierarchy used in the manuscript.
- `E-MTAB-1980` is the primary independent overall-survival validation cohort.
- `CPTAC` is reported as a development-linked overall-survival projection cohort because archived search-stage scripts used CPTAC during model-definition ranking.
- `IMmotion150` and `E-MTAB-3267` are secondary endpoint/context sensitivity analyses because they use progression-free survival in therapy-specific settings.
