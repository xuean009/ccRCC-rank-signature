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

suppressPackageStartupMessages(library(data.table))

DIR_DATA <- file.path(ROOT, "data")

# TCGA genes
expr_file <- file.path(DIR_DATA, "TCGA_Xena", "TCGA.KIRC.sampleMap", "HiSeqV2.gz")
expr <- fread(expr_file)
setnames(expr, 1, "gene")
genes_tcga <- unique(toupper(expr$gene))

# GSE29609 genes via annotation
series <- readLines(gzfile(file.path(DIR_DATA,'GEO','GSE29609_series_matrix.txt.gz')), warn=FALSE)
b <- which(grepl('!series_matrix_table_begin', series)); e <- which(grepl('!series_matrix_table_end', series))
tmp <- tempfile(fileext='.tsv'); writeLines(series[(b+1):(e-1)], tmp)
expr29609 <- fread(tmp)
setnames(expr29609,1,'ID_REF')
ann_lines <- readLines(gzfile(file.path(DIR_DATA,'GEO','annot','GPL1708.annot.gz')), warn=FALSE)
b2 <- which(grepl('!platform_table_begin', ann_lines)); e2 <- which(grepl('!platform_table_end', ann_lines))
a_tmp <- tempfile(fileext='.tsv'); writeLines(ann_lines[(b2+1):(e2-1)], a_tmp)
ann <- fread(a_tmp)
setnames(ann, old=c('ID','Gene symbol'), new=c('ID_REF','GeneSymbol'), skip_absent=TRUE)
ann <- ann[!is.na(GeneSymbol) & GeneSymbol!='', .(ID_REF=as.character(ID_REF), GeneSymbol=toupper(as.character(GeneSymbol)))]
expr29609$ID_REF <- as.character(expr29609$ID_REF)
mdt <- merge(expr29609[,.(ID_REF)], ann, by='ID_REF')
genes_29609 <- unique(mdt$GeneSymbol)

common <- intersect(genes_tcga, genes_29609)
cat('common genes:', length(common), '\n')

# Pick top variable genes in TCGA among common (variance across samples)
mat <- as.matrix(expr[match(common, toupper(expr$gene)), -1, with=FALSE])
rownames(mat) <- common
mode(mat) <- 'numeric'
# variance
v <- apply(mat, 1, var, na.rm=TRUE)
v[is.na(v)] <- 0
sel <- names(sort(v, decreasing=TRUE))[1:min(1500, length(v))]

out <- file.path(ROOT, "data", "gene_lists")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
fwrite(data.table(gene = sel), file.path(out, "common_topvar_1500.tsv"), sep = "\t")
cat("Wrote", file.path(out, "common_topvar_1500.tsv"), "\n")
