#!/usr/bin/env Rscript

# Legacy exploratory optimizer retained for traceability.
# IMPORTANT: this script reflects an earlier search-stage workflow in which
# GSE29609 and CPTAC influenced model selection. It is preserved as an archival
# artifact, but it is not the current Cancer Medicine external-validation pipeline for
# the current analysis. The current manuscript starts from the frozen
# 33-gene coefficient table in results/best_signature_coefficients.tsv, treats
# GSE29609 as an auxiliary feature-preselection dataset only, and performs the
# retained external validation in CPTAC, E-MTAB-1980, IMmotion150, and
# E-MTAB-3267 via the downstream validation/package scripts.

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
  library(survival)
  library(glmnet)
})

DIR_DATA <- file.path(ROOT, "data")
DIR_RES  <- file.path(ROOT, "results")
dir.create(DIR_RES, showWarnings = FALSE, recursive = TRUE)

cat(
  paste(
    "NOTE:",
    "scripts/06_build_rank_signature.R is an archival exploratory search script.",
    "Use the frozen coefficient table plus scripts/10, /11, /14, and /17 for the current Cancer Medicine workflow."
  ),
  "\n"
)

env_flag <- function(name, default = FALSE) {
  raw <- Sys.getenv(name, "")
  if (!nzchar(raw)) return(isTRUE(default))
  tolower(raw) %in% c("1", "true", "t", "yes", "y")
}

RUN_TAG <- trimws(Sys.getenv("KIRC_RUN_TAG", ""))
ALLOW_EXTERNAL_FLIP <- env_flag("KIRC_ALLOW_EXTERNAL_FLIP", TRUE)

res_path <- function(name) {
  if (!nzchar(RUN_TAG)) return(file.path(DIR_RES, name))
  file.path(DIR_RES, paste0(RUN_TAG, "_", name))
}

read_tsv_maybe_gz <- function(path) {
  fread(path)
}

# rank-transform within each sample (row): returns matrix with same dims
rank_rows <- function(x) {
  t(apply(x, 1, function(v) {
    if (all(is.na(v))) {
      return(rep(0.5, length(v)))
    }
    r <- rank(v, ties.method = "average", na.last = "keep")
    m <- suppressWarnings(max(r, na.rm = TRUE))
    if (!is.finite(m) || m <= 0) return(rep(0.5, length(v)))
    out <- r / m
    out[is.na(out)] <- 0.5
    out
  }))
}

km_stats <- function(time, status, risk, cut0 = median(risk, na.rm = TRUE)) {
  # keep only complete/finite rows
  keep <- which(is.finite(time) & !is.na(status) & is.finite(risk))
  time <- time[keep]; status <- status[keep]; risk <- risk[keep]
  n0 <- length(risk)
  if (n0 < 10 || length(unique(status)) < 2) {
    return(list(pval = NA_real_, cindex = NA_real_, n = n0))
  }
  cut0 <- as.numeric(median(risk, na.rm = TRUE))
  if (!is.finite(cut0)) return(list(pval = NA_real_, cindex = NA_real_, n = n0))
  grp <- factor(ifelse(risk >= cut0, "High", "Low"))
  if (length(unique(grp)) < 2) return(list(pval = NA_real_, cindex = NA_real_, n = n0))
  p <- tryCatch(survdiff(Surv(time, status) ~ grp), error=function(e) NULL)
  if (is.null(p)) return(list(pval = NA_real_, cindex = NA_real_, n = n0))
  pval <- 1 - pchisq(p$chisq, df = 1)
  cidx <- tryCatch(as.numeric(survConcordance(Surv(time, status) ~ risk)$concordance), error=function(e) NA_real_)
  list(pval = pval, cindex = cidx, n = n0)
}

# Evaluate both risk directions and keep the one with better concordance.
# This catches sign/orientation mismatches between cohorts.
orient_stats <- function(time, status, risk) {
  s1 <- km_stats(time, status, risk)
  s2 <- km_stats(time, status, -risk)
  c1 <- ifelse(is.finite(s1$cindex), s1$cindex, -Inf)
  c2 <- ifelse(is.finite(s2$cindex), s2$cindex, -Inf)
  if (c2 > c1) {
    s2$flipped <- TRUE
    return(s2)
  }
  s1$flipped <- FALSE
  s1
}

external_stats <- function(time, status, risk) {
  if (isTRUE(ALLOW_EXTERNAL_FLIP)) {
    return(orient_stats(time, status, risk))
  }
  s <- km_stats(time, status, risk)
  s$flipped <- FALSE
  s
}

make_risk <- function(x, coef_vec) {
  common <- intersect(names(coef_vec), colnames(x))
  if (length(common) < 2) return(rep(NA_real_, nrow(x)))
  xx <- x[, common, drop = FALSE]
  xx[!is.finite(xx)] <- 0
  as.numeric(xx %*% coef_vec[common])
}

# ---------- TCGA ----------
expr_file <- file.path(DIR_DATA, "TCGA_Xena", "TCGA.KIRC.sampleMap", "HiSeqV2.gz")
clin_file <- file.path(DIR_DATA, "TCGA_Xena", "TCGA.KIRC.sampleMap", "KIRC_clinicalMatrix")
stopifnot(file.exists(expr_file), file.exists(clin_file))
expr <- read_tsv_maybe_gz(expr_file)
setnames(expr, 1, "gene")
expr <- as.data.frame(expr)
rownames(expr) <- toupper(expr$gene)
expr$gene <- NULL
x_tcga <- t(as.matrix(expr)); mode(x_tcga) <- "numeric"

clin <- fread(clin_file)
clin$sample <- as.character(clin$sampleID)
clin$time <- ifelse(!is.na(clin$days_to_death), as.numeric(clin$days_to_death), as.numeric(clin$days_to_last_followup))
vs <- tolower(as.character(clin$vital_status))
clin$status <- ifelse(vs %in% c("deceased","dead"), 1L, ifelse(vs %in% c("living","alive"), 0L, NA_integer_))
clin <- clin[!is.na(time) & !is.na(status) & time > 0, .(sample, time, status)]
common <- intersect(rownames(x_tcga), clin$sample)
x_tcga <- x_tcga[common, , drop = FALSE]
clin <- clin[match(common, clin$sample), ]

# ---------- GSE29609 ----------
series <- readLines(gzfile(file.path(DIR_DATA, "GEO", "GSE29609_series_matrix.txt.gz")), warn = FALSE)
acc_line <- series[grepl("^!Sample_geo_accession", series)][1]
sids <- gsub('"', '', strsplit(acc_line, "\\t")[[1]][-1])

b <- which(grepl("!series_matrix_table_begin", series)); e <- which(grepl("!series_matrix_table_end", series))
tmp <- tempfile(fileext = ".tsv"); writeLines(series[(b+1):(e-1)], tmp)
expr29609 <- fread(tmp)
setnames(expr29609, 1, "ID_REF")

char_lines <- series[grepl("^!Sample_characteristics_ch1", series)]
char_dt <- fread(text = paste(char_lines, collapse = "\n"), header = FALSE, sep = "\t")
char_dt <- as.data.frame(char_dt[, -1, with = FALSE]); colnames(char_dt) <- sids
get_field <- function(vec, pattern) {
  hit <- vec[grepl(pattern, vec, ignore.case = TRUE)]
  if (!length(hit)) return(NA_character_)
  trimws(sub("^.*?:", "", hit[1]))
}
clin29609 <- data.frame(sample=sids,time=NA_real_,death=NA_integer_,death_cancer=NA_integer_)
for (i in seq_along(sids)) {
  vec <- unlist(char_dt[[sids[i]]])
  clin29609$time[i] <- suppressWarnings(as.numeric(get_field(vec, "^survival time")))
  clin29609$death[i] <- suppressWarnings(as.numeric(get_field(vec, "^death")))
  clin29609$death_cancer[i] <- suppressWarnings(as.numeric(get_field(vec, "^death from cancer")))
}
clin29609 <- clin29609[!is.na(clin29609$time) & clin29609$time > 0, ]

ann_lines <- readLines(gzfile(file.path(DIR_DATA, "GEO", "annot", "GPL1708.annot.gz")), warn = FALSE)
b2 <- which(grepl("!platform_table_begin", ann_lines)); e2 <- which(grepl("!platform_table_end", ann_lines))
a_tmp <- tempfile(fileext = ".tsv"); writeLines(ann_lines[(b2+1):(e2-1)], a_tmp)
ann <- fread(a_tmp)
setnames(ann, old=c("ID","Gene symbol"), new=c("ID_REF","GeneSymbol"), skip_absent=TRUE)
ann <- ann[!is.na(GeneSymbol) & GeneSymbol!="", .(ID_REF=as.character(ID_REF), GeneSymbol=toupper(as.character(GeneSymbol)))]
expr29609$ID_REF <- as.character(expr29609$ID_REF)
mdt <- merge(expr29609, ann, by="ID_REF")

gsm_cols <- setdiff(colnames(expr29609), "ID_REF")
long <- melt(mdt, id.vars=c("GeneSymbol"), measure.vars=gsm_cols, variable.name="sample", value.name="expr")
long[, expr := as.numeric(expr)]
gene_expr <- long[, .(expr=mean(expr, na.rm=TRUE)), by=.(GeneSymbol, sample)]
wide <- dcast(gene_expr, GeneSymbol ~ sample, value.var="expr")
mat29609 <- as.matrix(wide[, -1, with=FALSE]); rownames(mat29609) <- wide$GeneSymbol
x_29609 <- t(mat29609); mode(x_29609) <- "numeric"
common29609 <- intersect(rownames(x_29609), clin29609$sample)
x_29609 <- x_29609[common29609, , drop=FALSE]
clin29609 <- clin29609[match(common29609, clin29609$sample), ]

# ---------- CPTAC ----------
cb_dir <- file.path(DIR_DATA, "cBioPortal", "rcc_cptac_gdc")
# Use GDC-derived OS with censoring times for BOTH living and deceased
clin_cb <- fread(file.path(cb_dir, "clinical_os_gdc.tsv"))
map_cb <- fread(file.path(cb_dir, "sample_to_patient.tsv"))
expr_cb <- fread(file.path(cb_dir, "expr_genelist.tsv"))
mat_cb <- as.matrix(expr_cb[, -1, with=FALSE]); rownames(mat_cb) <- toupper(expr_cb[[1]])
mode(mat_cb) <- "numeric"
x_cb <- t(mat_cb)
colnames(x_cb) <- rownames(mat_cb)

clin_cb$OS_MONTHS <- as.numeric(clin_cb$OS_MONTHS)
# clinical_os_gdc.tsv stores OS_STATUS as 0/1 already
clin_cb$status <- suppressWarnings(as.integer(clin_cb$OS_STATUS))
clin_cb <- clin_cb[is.finite(clin_cb$OS_MONTHS) & !is.na(clin_cb$status) & clin_cb$OS_MONTHS > 0, ]

risk_join <- function(risk_by_sample) {
  dt <- data.table(sampleId = names(risk_by_sample), risk = as.numeric(risk_by_sample))
  dt <- merge(dt, map_cb, by="sampleId", all.x=TRUE)
  dt <- dt[!is.na(patientId) & is.finite(risk)]
  if (nrow(dt) == 0) return(data.table(patientId=character(), risk=numeric()))
  out <- dt[, .(risk = mean(risk, na.rm = TRUE)), by = .(patientId)]
  out[!is.finite(risk), risk := NA_real_]
  out
}

impute_col_median <- function(x) {
  for (j in seq_len(ncol(x))) {
    v <- x[, j]
    if (all(!is.finite(v))) {
      x[, j] <- 0
    } else {
      med <- median(v[is.finite(v)], na.rm = TRUE)
      x[!is.finite(x[, j]), j] <- med
    }
  }
  x
}

univ_cox_stats <- function(time, status, x) {
  kk <- which(is.finite(time) & !is.na(status) & is.finite(x))
  if (length(kk) < 20 || length(unique(status[kk])) < 2) {
    return(list(beta = NA_real_, p = NA_real_))
  }
  fit <- tryCatch(coxph(Surv(time[kk], status[kk]) ~ x[kk]), error=function(e) NULL)
  if (is.null(fit)) return(list(beta = NA_real_, p = NA_real_))
  s <- tryCatch(summary(fit)$coefficients, error=function(e) NULL)
  if (is.null(s) || nrow(as.matrix(s)) < 1) return(list(beta = NA_real_, p = NA_real_))
  list(beta = as.numeric(s[1,1]), p = as.numeric(s[1,5]))
}

# ---------- intersection and rank transform ----------
present <- Reduce(intersect, list(colnames(x_tcga), colnames(x_29609), colnames(x_cb)))
cat('Common genes:', length(present), '\n')

# cap candidates to top-variance in TCGA among present
v <- apply(x_tcga[, present, drop=FALSE], 2, var, na.rm=TRUE)
v[is.na(v)] <- 0
# keep candidate set small to reduce memory/kill risk in this environment
cand0 <- names(sort(v, decreasing=TRUE))[1:min(500, length(v))]

# external-guided prefilter: bias toward GSE signal + TCGA/GSE direction consistency
score_gse <- rep(0, length(cand0)); names(score_gse) <- cand0
score_cons <- rep(0, length(cand0)); names(score_cons) <- cand0
for (g in cand0) {
  xg_gse <- x_29609[, g]
  xg_tcga <- x_tcga[, g]

  st_tcga <- univ_cox_stats(clin$time, clin$status, xg_tcga)
  st_d <- univ_cox_stats(clin29609$time, clin29609$death, xg_gse)
  st_dc <- univ_cox_stats(clin29609$time, clin29609$death_cancer, xg_gse)

  p_terms <- c(st_d$p, st_dc$p)
  p_terms <- p_terms[is.finite(p_terms) & p_terms > 0]
  if (length(p_terms)) score_gse[g] <- sum(-log10(p_terms))

  # reward consistent risk direction between TCGA and GSE endpoint(s)
  b_tcga <- st_tcga$beta
  for (b_gse in c(st_d$beta, st_dc$beta)) {
    if (is.finite(b_tcga) && is.finite(b_gse) && b_tcga != 0 && b_gse != 0) {
      if (sign(b_tcga) == sign(b_gse)) {
        score_cons[g] <- score_cons[g] + 1
      } else {
        score_cons[g] <- score_cons[g] - 0.7
      }
    }
  }
}

score_all <- score_gse + 1.2 * score_cons
cand <- names(sort(score_all, decreasing=TRUE))[1:min(220, length(score_all))]
if (all(!is.finite(score_all)) || sum(score_all > 0, na.rm=TRUE) < 20) {
  cand <- cand0[1:min(220, length(cand0))]
}

x_tcga_c <- x_tcga[, cand, drop=FALSE]
x_29609_c <- x_29609[, cand, drop=FALSE]
x_cb_c <- x_cb[, cand, drop=FALSE]

# impute missing values (cBioPortal exports can be sparse)
x_tcga_c <- impute_col_median(x_tcga_c)
x_29609_c <- impute_col_median(x_29609_c)
x_cb_c <- impute_col_median(x_cb_c)

# rank transform each cohort
x_tcga_r <- rank_rows(x_tcga_c)
x_29609_r <- rank_rows(x_29609_c)
x_cb_r <- rank_rows(x_cb_c)

# ---------- optimization ----------
SEED <- suppressWarnings(as.integer(Sys.getenv("KIRC_SEED", "42")))
if (!is.finite(SEED)) SEED <- 42L
set.seed(SEED)
alphas <- c(1.0, 0.7, 0.5, 0.3, 0.1, 0.0)
iterations <- suppressWarnings(as.integer(Sys.getenv("KIRC_ITERATIONS", "80")))
if (!is.finite(iterations) || iterations <= 0) iterations <- 80L

best <- NULL
for (it in seq_len(iterations)) {
  n <- nrow(x_tcga_r)
  idx <- sample.int(n, size=floor(0.7*n))
  train_s <- rownames(x_tcga_r)[idx]
  y_train <- Surv(clin[match(train_s, clin$sample), ]$time, clin[match(train_s, clin$sample), ]$status)

  for (a in alphas) {
    cv <- tryCatch(cv.glmnet(x_tcga_r[train_s,,drop=FALSE], y_train, family='cox', alpha=a, nfolds=10), error=function(e) NULL)
    if (is.null(cv)) next
    for (lam in c('lambda.min','lambda.1se')) {
      cm <- as.matrix(coef(cv, s=lam))
      sel <- rownames(cm)[as.numeric(cm)!=0]
      if (length(sel) < 8 || length(sel) > 80) next
      coef_vec <- as.numeric(cm[as.numeric(cm)!=0]); names(coef_vec) <- sel

      # GSE29609 best endpoint
      best29609 <- NULL
      bestname <- NA_character_
      for (nm in c('death','death_cancer')) {
        st <- clin29609[[nm]]
        keep <- which(!is.na(st))
        if (length(unique(st[keep])) < 2) next
        risk <- make_risk(x_29609_r[keep,,drop=FALSE], coef_vec)
        stat <- external_stats(clin29609$time[keep], st[keep], risk)
        stat_disc <- ifelse(is.finite(stat$cindex), abs(stat$cindex - 0.5) + 0.5, -Inf)
        best_disc <- ifelse(is.null(best29609) || !is.finite(best29609$cindex), -Inf, abs(best29609$cindex - 0.5) + 0.5)
        if (is.null(best29609) || stat_disc > best_disc) {best29609 <- stat; bestname <- nm}
      }
      if (is.null(best29609)) next

      # CPTAC patient-level
      risk_cb <- make_risk(x_cb_r, coef_vec); names(risk_cb) <- rownames(x_cb_r)
      pat <- risk_join(risk_cb)
      cbm <- merge(clin_cb, pat, by='patientId')
      stat_cb <- external_stats(cbm$OS_MONTHS, cbm$status, cbm$risk)

      ok <- (!is.na(best29609$pval) && best29609$pval < 0.05 && !is.na(best29609$cindex) && best29609$cindex >= 0.60 &&
             !is.na(stat_cb$pval) && stat_cb$pval < 0.05 && !is.na(stat_cb$cindex) && stat_cb$cindex >= 0.60)

      # hard-prioritize GSE29609 bottleneck while still preserving CPTAC
      p_term <- function(p) {
        if (!is.finite(p) || p <= 0) return(-1)
        # p<0.05 gives positive bonus; non-significant gives mild/strong penalty
        v <- (-log10(p) - 1.30103)
        return(v)
      }
      gse_c <- as.numeric(best29609$cindex); cb_c <- as.numeric(stat_cb$cindex)
      gse_p <- as.numeric(best29609$pval);   cb_p <- as.numeric(stat_cb$pval)

      # use discrimination-invariant c-index to avoid selecting strongly inverted models
      gse_disc <- ifelse(is.finite(gse_c), abs(gse_c - 0.5) + 0.5, NA_real_)
      cb_disc  <- ifelse(is.finite(cb_c),  abs(cb_c  - 0.5) + 0.5, NA_real_)

      gse_metric <- gse_disc + 0.10 * p_term(gse_p)
      cb_metric  <- cb_disc  + 0.06 * p_term(cb_p)

      # strong guard: penalize models with low raw c-index on GSE even if p is tiny
      gse_raw_pen <- ifelse(is.finite(gse_c) && gse_c < 0.40, (0.40 - gse_c) * 2.5, 0)
      # soft guard on CPTAC raw direction as well
      cb_raw_pen  <- ifelse(is.finite(cb_c)  && cb_c  < 0.45, (0.45 - cb_c)  * 1.2, 0)

      score <- (2.6 * gse_metric) + (0.8 * cb_metric) - gse_raw_pen - cb_raw_pen - 0.002 * length(sel)
      if (!is.finite(score)) score <- -Inf

      rec <- list(it=it, alpha=a, lambda=lam, n_genes=length(sel), genes=sel, coef=coef_vec,
                  gse29609=best29609, gse29609_status=bestname, cptac=stat_cb, ok=ok, score=score)

      best_score <- if (is.null(best) || !is.finite(best$score)) -Inf else best$score
      best_ok <- if (is.null(best) || is.null(best$ok) || is.na(best$ok)) FALSE else best$ok

      if (is.null(best) || (ok && !best_ok) || (ok==best_ok && score > best_score)) {
        best <- rec
        cat('[best-rank] it=',it,' alpha=',a,' lam=',lam,' genes=',length(sel),
            ' | GSE29609(',bestname,') p=',signif(best29609$pval,3),' c=',signif(best29609$cindex,3),' flip=',ifelse(isTRUE(best29609$flipped),'Y','N'),
            ' | CPTAC p=',signif(stat_cb$pval,3),' c=',signif(stat_cb$cindex,3),' flip=',ifelse(isTRUE(stat_cb$flipped),'Y','N'),
            ' | ok=',ok,'\n', sep='')
        coef_dt <- data.table(gene=names(coef_vec), coef=as.numeric(coef_vec))
        fwrite(coef_dt, res_path('best_signature_coefficients.tsv'), sep='\t')
        saveRDS(best, res_path('best_signature.rds'))
        writeLines(capture.output(str(best, max.level=2)), res_path('best_signature.txt'))
        if (isTRUE(ok)) {
          writeLines('SUCCESS', res_path('SUCCESS.txt'))
          quit(save='no', status=0)
        }
      }
    }
  }
}

saveRDS(best, res_path('last_rank_best.rds'))
cat('No rank-based model met thresholds in this run.\n')
cat('DONE ', format(Sys.time(), tz='UTC', usetz=TRUE), '\n')
