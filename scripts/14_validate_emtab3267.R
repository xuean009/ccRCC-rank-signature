#!/usr/bin/env Rscript

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

ROOT <- normalizePath(file.path(get_script_dir(), ".."), winslash = "/", mustWork = FALSE)
local_lib <- file.path(ROOT, "r_lib")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(data.table)
  library(oligo)
  library(pd.hugene.1.0.st.v1)
  library(hugene10sttranscriptcluster.db)
  library(AnnotationDbi)
  library(survival)
})

options(timeout = max(7200, getOption("timeout")))

DIR_DATA <- file.path(ROOT, "data")
DIR_RES <- file.path(ROOT, "results")
DIR_EXT <- file.path(DIR_DATA, "external_clinical", "E-MTAB-3267")
DIR_CEL <- file.path(DIR_EXT, "raw_cel")
dir.create(DIR_EXT, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_CEL, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_RES, showWarnings = FALSE, recursive = TRUE)

download_if_missing <- function(url, dest) {
  if (file.exists(dest) && file.info(dest)$size > 0) return(invisible(dest))
  download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
  dest
}

download_missing_cels_from_fire <- function(cel_names, exdir) {
  base_url <- "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/267/E-MTAB-3267/Files"
  for (nm in cel_names) {
    dest <- file.path(exdir, nm)
    if (file.exists(dest) && file.info(dest)$size > 0) next
    url <- sprintf("%s/%s", base_url, nm)
    message("Downloading missing CEL from BioStudies fire: ", nm)
    download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
  }
  invisible(file.path(exdir, cel_names))
}

extract_zip_if_needed <- function(zip_path, exdir) {
  listed <- utils::unzip(zip_path, list = TRUE)
  targets <- listed$Name[grepl("\\.CEL$", listed$Name, ignore.case = TRUE)]
  missing <- targets[!file.exists(file.path(exdir, targets))]
  if (length(missing)) {
    utils::unzip(zip_path, exdir = exdir)
  }
  invisible(file.path(exdir, targets))
}

rank_rows <- function(x) {
  t(apply(x, 1, function(v) {
    if (all(is.na(v))) return(rep(0.5, length(v)))
    r <- rank(v, ties.method = "average", na.last = "keep")
    m <- suppressWarnings(max(r, na.rm = TRUE))
    if (!is.finite(m) || m <= 0) return(rep(0.5, length(v)))
    out <- r / m
    out[is.na(out)] <- 0.5
    out
  }))
}

km_stats <- function(time, status, risk, cut0 = median(risk, na.rm = TRUE)) {
  keep <- which(is.finite(time) & !is.na(status) & is.finite(risk))
  time <- time[keep]; status <- status[keep]; risk <- risk[keep]
  n0 <- length(risk)
  if (n0 < 10 || length(unique(status)) < 2) {
    return(list(pval = NA_real_, cindex = NA_real_, n = n0))
  }
  cut0 <- as.numeric(median(risk, na.rm = TRUE))
  grp <- factor(ifelse(risk >= cut0, "High", "Low"))
  if (length(unique(grp)) < 2) {
    return(list(pval = NA_real_, cindex = NA_real_, n = n0))
  }
  p <- tryCatch(survdiff(Surv(time, status) ~ grp), error = function(e) NULL)
  if (is.null(p)) return(list(pval = NA_real_, cindex = NA_real_, n = n0))
  pval <- 1 - pchisq(p$chisq, df = 1)
  cidx <- tryCatch(as.numeric(concordance(Surv(time, status) ~ risk, reverse = TRUE)$concordance), error = function(e) NA_real_)
  list(pval = pval, cindex = cidx, n = n0)
}

orientation_metrics <- function(time, status, risk) {
  s1 <- km_stats(time, status, risk)
  s2 <- km_stats(time, status, -risk)
  c1 <- ifelse(is.finite(s1$cindex), s1$cindex, -Inf)
  c2 <- ifelse(is.finite(s2$cindex), s2$cindex, -Inf)
  list(
    original = s1,
    sign_reversed = s2,
    best_discrimination_orientation = ifelse(c2 > c1, "sign_reversed", "original")
  )
}

sdrf_path <- file.path(DIR_EXT, "E-MTAB-3267.sdrf.txt")
idf_path <- file.path(DIR_EXT, "E-MTAB-3267.idf.txt")
download_if_missing("https://www.ebi.ac.uk/biostudies/files/E-MTAB-3267/E-MTAB-3267.sdrf.txt", sdrf_path)
download_if_missing("https://www.ebi.ac.uk/biostudies/files/E-MTAB-3267/E-MTAB-3267.idf.txt", idf_path)

zip1 <- file.path(DIR_EXT, "E-MTAB-3267.raw.1.zip")
zip2 <- file.path(DIR_EXT, "E-MTAB-3267.raw.2.zip")
download_if_missing("https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3267/E-MTAB-3267.raw.1.zip", zip1)
download_if_missing("https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3267/E-MTAB-3267.raw.2.zip", zip2)

extract_zip_if_needed(zip1, DIR_CEL)
extract_zip_if_needed(zip2, DIR_CEL)

sdrf <- fread(sdrf_path)
setnames(
  sdrf,
  old = c(
    "Source Name",
    "Assay Name",
    "Array Data File",
    "Factor Value[individual]",
    "Characteristics[disease]",
    "Characteristics[histology type]",
    "Characteristics[sex]",
    "Characteristics[age]",
    "Characteristics[sunitinib response]",
    "Characteristics[progression]",
    "Characteristics[progression free survival]"
  ),
  new = c(
    "source_name",
    "assay_name",
    "array_data_file",
    "individual",
    "disease",
    "histology_type",
    "sex",
    "age",
    "sunitinib_response",
    "progression",
    "progression_free_survival"
  ),
  skip_absent = TRUE
)

tumor <- sdrf[
  disease == "Tumor" &
    histology_type == "Clear Cell" &
    !is.na(progression_free_survival),
  .(
    sample_id = sub("\\.CEL$", "", array_data_file, ignore.case = TRUE),
    source_name = source_name,
    assay_name = assay_name,
    individual = individual,
    sex = sex,
    age = suppressWarnings(as.numeric(age)),
    sunitinib_response = sunitinib_response,
    status_pfs = suppressWarnings(as.integer(progression)),
    pfs_months = suppressWarnings(as.numeric(progression_free_survival))
  )
]

cel_files <- file.path(DIR_CEL, paste0(tumor$sample_id, ".CEL"))
missing <- basename(cel_files[!file.exists(cel_files)])
if (length(missing)) {
  download_missing_cels_from_fire(missing, DIR_CEL)
}
if (!all(file.exists(cel_files))) {
  missing <- cel_files[!file.exists(cel_files)]
  stop("Missing CEL files after extraction: ", paste(basename(missing), collapse = ", "))
}

eset_cache <- file.path(DIR_EXT, "rma_expression_all_samples.rds")
if (file.exists(eset_cache)) {
  eset <- readRDS(eset_cache)
} else {
  raw <- read.celfiles(cel_files, pkgname = "pd.hugene.1.0.st.v1")
  sampleNames(raw) <- tumor$sample_id
  eset <- rma(raw, target = "core")
  saveRDS(eset, eset_cache)
}

expr <- exprs(eset)
annot <- AnnotationDbi::select(
  hugene10sttranscriptcluster.db,
  keys = rownames(expr),
  columns = c("SYMBOL"),
  keytype = "PROBEID"
)
annot <- as.data.table(annot)
annot <- annot[!is.na(SYMBOL) & SYMBOL != ""]
annot[, SYMBOL := toupper(SYMBOL)]
annot <- unique(annot, by = c("PROBEID", "SYMBOL"))

expr_dt <- as.data.table(expr, keep.rownames = "PROBEID")
expr_anno <- merge(annot, expr_dt, by = "PROBEID", all = FALSE)
sample_cols <- setdiff(colnames(expr_anno), c("PROBEID", "SYMBOL"))
expr_gene <- expr_anno[, lapply(.SD, mean, na.rm = TRUE), by = SYMBOL, .SDcols = sample_cols]
setnames(expr_gene, "SYMBOL", "gene")

coef <- fread(file.path(DIR_RES, "best_signature_coefficients.tsv"))
coef[, gene := toupper(gene)]
sig_genes <- coef$gene

expr_sig <- expr_gene[gene %in% sig_genes]
expr_sig <- expr_sig[match(sig_genes, gene)]
expr_sig <- expr_sig[!is.na(gene)]
fwrite(expr_sig, file.path(DIR_EXT, "expression_signature_gene33.tsv"), sep = "\t")

clin_out <- copy(tumor)
fwrite(clin_out, file.path(DIR_EXT, "clinical_matched_53.tsv"), sep = "\t")

mat <- as.matrix(expr_sig[, ..sample_cols])
rownames(mat) <- expr_sig$gene
x <- t(mat)
mode(x) <- "numeric"
x <- x[clin_out$sample_id, , drop = FALSE]
xr <- rank_rows(x[, coef$gene, drop = FALSE])
risk <- as.numeric(xr %*% coef$coef)
metrics <- orientation_metrics(clin_out$pfs_months, clin_out$status_pfs, risk)
stat <- metrics$original

summary_dt <- data.table(
  cohort = "E-MTAB-3267",
  endpoint = "PFS",
  n_samples = nrow(clin_out),
  n_events = sum(clin_out$status_pfs, na.rm = TRUE),
  candidate_genes_for_ranking = ncol(xr),
  signature_genes_used = nrow(coef),
  manuscript_orientation = "original",
  original_logrank_p = stat$pval,
  original_cindex = stat$cindex,
  sign_reversed_logrank_p = metrics$sign_reversed$pval,
  sign_reversed_cindex = metrics$sign_reversed$cindex,
  best_discrimination_orientation = metrics$best_discrimination_orientation,
  logrank_p = stat$pval,
  cindex = stat$cindex,
  median_cut_high_n = sum(risk >= median(risk, na.rm = TRUE)),
  median_cut_low_n = sum(risk < median(risk, na.rm = TRUE))
)

score_dt <- copy(clin_out)
score_dt[, risk_score := risk]
score_dt[, risk_score_sign_reversed := -risk]

fwrite(summary_dt, file.path(DIR_RES, "external_validation_emtab3267_summary.tsv"), sep = "\t")
fwrite(score_dt, file.path(DIR_RES, "external_validation_emtab3267_scores.tsv"), sep = "\t")
writeLines(jsonlite::toJSON(as.list(summary_dt[1]), pretty = TRUE, auto_unbox = TRUE), file.path(DIR_RES, "external_validation_emtab3267_summary.json"))

readme <- c(
  "# E-MTAB-3267 external validation cache",
  "",
  "Recovered on: 2026-03-25",
  "",
  "## Public sources",
  "",
  "- ArrayExpress/BioStudies study page: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-3267",
  "- SDRF: https://www.ebi.ac.uk/biostudies/files/E-MTAB-3267/E-MTAB-3267.sdrf.txt",
  "- PMID: https://pubmed.ncbi.nlm.nih.gov/25583177/",
  "",
  "## What is stored here",
  "",
  "- clinical_matched_53.tsv: 53 tumor samples with sunitinib response, progression status, and progression-free survival.",
  "- expression_signature_gene33.tsv: RMA-normalized Human Gene 1.0 ST expression collapsed to the frozen 33-gene signature.",
  "- rma_expression_all_samples.rds: cached normalized expression object generated from raw CEL files.",
  "",
  "## Notes",
  "",
  "- This is a metastatic clear-cell RCC sunitinib-treated cohort.",
  "- The appropriate endpoint is PFS/progression, not OS.",
  "- Raw CEL files were normalized locally with oligo::rma(target = 'core').",
  "- Cancer Medicine summaries retain the original score direction; sign-reversed metrics are supplementary robustness diagnostics only.",
  "- Missing CEL files not present in the legacy raw zip bundles were recovered from the BioStudies fire Files endpoint."
)
writeLines(readme, file.path(DIR_EXT, "README.md"))

cat(as.character(jsonlite::toJSON(as.list(summary_dt[1]), pretty = TRUE, auto_unbox = TRUE)), "\n")
