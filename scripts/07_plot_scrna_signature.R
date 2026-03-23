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
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

dir_data <- file.path(ROOT, "data", "scRNA", "GSE159115")
dir_raw <- file.path(dir_data, "raw")
dir_out <- file.path(ROOT, "output", "scst", "scRNA_GSE159115")
dir.create(dir_raw, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

raw_tar <- file.path(dir_data, "GSE159115_RAW.tar")
anno_gz <- file.path(dir_data, "GSE159115_ccRCC_anno.csv.gz")

if (!file.exists(raw_tar)) {
  download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159115/suppl/GSE159115_RAW.tar", raw_tar, mode = "wb", quiet = FALSE)
}
if (!file.exists(anno_gz)) {
  download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159115/suppl/GSE159115_ccRCC_anno.csv.gz", anno_gz, mode = "wb", quiet = FALSE)
}

h5_files <- list.files(dir_raw, pattern = "\\.h5$", full.names = TRUE)
if (length(h5_files) == 0) {
  untar(raw_tar, exdir = dir_raw)
  h5_files <- list.files(dir_raw, pattern = "\\.h5$", full.names = TRUE)
}
if (length(h5_files) == 0) stop("No .h5 files found after extraction")

sample_from_h5 <- function(x) {
  b <- basename(x)
  # e.g., GSM4819725_SI_18854_filtered_gene_bc_matrices_h5.h5 -> SI_18854
  sid <- sub("^GSM[0-9]+_(SI_[0-9]+)_.*$", "\\1", b)
  if (identical(sid, b)) sid <- sub("^(GSM[0-9]+_[^_]+_[^_]+)_.*$", "\\1", b)
  sid
}

max_samples <- suppressWarnings(as.integer(Sys.getenv("SC_MAX_SAMPLES", "8")))
if (!is.finite(max_samples) || max_samples <= 0) max_samples <- 8
max_cells_per_sample <- suppressWarnings(as.integer(Sys.getenv("SC_MAX_CELLS_PER_SAMPLE", "4000")))
if (!is.finite(max_cells_per_sample) || max_cells_per_sample <= 0) max_cells_per_sample <- 4000

h5_files <- h5_files[seq_len(min(length(h5_files), max_samples))]

seu_list <- list()
for (hf in h5_files) {
  sid <- sample_from_h5(hf)
  cat("Reading", sid, "from", basename(hf), "\n")
  m <- Read10X_h5(hf)
  if (is.list(m)) {
    # prefer Gene Expression if multiple modalities
    if ("Gene Expression" %in% names(m)) {
      m <- m[["Gene Expression"]]
    } else {
      m <- m[[1]]
    }
  }
  so <- CreateSeuratObject(counts = m, project = sid, min.cells = 3, min.features = 200)
  so <- RenameCells(so, new.names = paste0(sid, "_", colnames(so)))
  if (ncol(so) > max_cells_per_sample) {
    set.seed(42)
    keep <- sample(colnames(so), max_cells_per_sample)
    so <- subset(so, cells = keep)
  }
  so$sample <- sid
  seu_list[[sid]] <- so
  gc()
}

if (length(seu_list) == 1) {
  seu <- seu_list[[1]]
} else {
  seu <- Reduce(function(x, y) merge(x, y), seu_list)
}

# Seurat v5 merged objects may contain multiple layers; join for downstream scoring APIs.
if (exists("JoinLayers", mode = "function")) {
  seu <- tryCatch(JoinLayers(seu), error = function(e) seu)
}

# QC
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)

# Attach GEO-provided annotation if available
anno <- tryCatch(fread(anno_gz), error = function(e) NULL)
if (!is.null(anno) && all(c("cell", "anno") %in% colnames(anno))) {
  annomap <- anno[, .(cell, anno, patient, label)]
  annomap <- annomap[match(colnames(seu), cell)]
  seu$anno <- annomap$anno
  seu$patient <- annomap$patient
  seu$label <- annomap$label
}

# Fallback marker-based annotation if GEO anno missing
if (!"anno" %in% colnames(seu@meta.data) || all(is.na(seu$anno))) {
  marker_sets <- list(
    Tumor = c("CA9", "VIM", "EPCAM", "KRT8", "KRT18"),
    Tcell = c("CD3D", "CD3E", "CD8A"),
    Macro = c("CD68", "CD163", "LST1"),
    Endothelial = c("PECAM1", "VWF", "KDR")
  )
  seu <- NormalizeData(seu)
  seu <- AddModuleScore(seu, features = marker_sets, name = "CTYPE")
  sc_cols <- grep("^CTYPE[0-9]+$", colnames(seu@meta.data), value = TRUE)
  mm <- as.matrix(seu@meta.data[, sc_cols, drop = FALSE])
  colnames(mm) <- names(marker_sets)[seq_len(ncol(mm))]
  cls <- colnames(mm)[max.col(mm, ties.method = "first")]
  seu$anno <- cls
}

# Standard Seurat workflow
if (!"RNA" %in% Assays(seu)) DefaultAssay(seu) <- Assays(seu)[1]
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
seu <- ScaleData(seu, features = VariableFeatures(seu))
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 40)
seu <- RunUMAP(seu, dims = 1:30)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)

# Signature genes from current prognostic model
coef_file <- file.path(ROOT, "results", "best_signature_coefficients.tsv")
coef_dt <- fread(coef_file)
coef_dt[, abs_coef := abs(coef)]
core_genes <- coef_dt[order(-abs_coef)]$gene[1:min(5, nrow(coef_dt))]
core_genes <- core_genes[core_genes %in% rownames(seu)]
sig_genes <- unique(coef_dt$gene[coef_dt$gene %in% rownames(seu)])

if (length(sig_genes) >= 3) {
  seu <- AddModuleScore(seu, features = list(sig_genes), name = "IMRISK")
}

# Save core outputs first (so downstream steps can continue even if plotting fails)
saveRDS(seu, file.path(dir_out, "seurat_scRNA_gse159115.rds"))
fwrite(as.data.table(seu@meta.data, keep.rownames = "cell"), file.path(dir_out, "cell_metadata.tsv"), sep = "\t")

# Plots (best-effort)
try({
  p_umap <- DimPlot(seu, reduction = "umap", group.by = "anno", label = TRUE, repel = TRUE) + ggtitle("GSE159115 cell annotation")
  ggsave(file.path(dir_out, "umap_celltype.png"), p_umap, width = 9, height = 7, dpi = 180)
}, silent = TRUE)

if (length(core_genes) > 0) {
  try({
    pf <- FeaturePlot(seu, features = core_genes, reduction = "umap", ncol = min(3, length(core_genes)))
    ggsave(file.path(dir_out, "featureplot_core_genes.png"), pf, width = 12, height = 8, dpi = 180)
  }, silent = TRUE)

  try({
    pv <- VlnPlot(seu, features = core_genes, group.by = "anno", pt.size = 0, ncol = min(2, length(core_genes)))
    ggsave(file.path(dir_out, "vlnplot_core_genes_by_celltype.png"), pv, width = 12, height = 8, dpi = 180)
  }, silent = TRUE)
}

if ("IMRISK1" %in% colnames(seu@meta.data)) {
  try({
    p_risk <- VlnPlot(seu, features = "IMRISK1", group.by = "anno", pt.size = 0) + ggtitle("Immune-metabolic risk score by cell type")
    ggsave(file.path(dir_out, "vlnplot_signature_score_by_celltype.png"), p_risk, width = 10, height = 6, dpi = 180)
  }, silent = TRUE)

  if ("patient" %in% colnames(seu@meta.data)) {
    med <- aggregate(IMRISK1 ~ patient, data = seu@meta.data, FUN = median)
    fwrite(as.data.table(med), file.path(dir_out, "signature_score_by_patient.tsv"), sep = "\t")
  }
}

cat("Done scRNA analysis. Cells:", ncol(seu), " Genes:", nrow(seu), "\n")
