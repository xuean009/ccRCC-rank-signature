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
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

dir_data <- file.path(ROOT, "data", "spatial", "GSE250163")
dir_raw <- file.path(dir_data, "raw")
dir_vis <- file.path(dir_data, "visium")
dir_out <- file.path(ROOT, "output", "scst", "spatial_GSE250163")
dir.create(dir_raw, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_vis, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

raw_tar <- file.path(dir_data, "GSE250163_RAW.tar")
if (!file.exists(raw_tar)) {
  download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250163/suppl/GSE250163_RAW.tar", raw_tar, mode = "wb", quiet = FALSE)
}

# Extract once
if (length(list.files(dir_raw, pattern = "filtered_feature_bc_matrix\\.h5$", recursive = TRUE)) == 0) {
  untar(raw_tar, exdir = dir_raw)
}

all_h5 <- list.files(dir_raw, pattern = "filtered_feature_bc_matrix\\.h5$", full.names = TRUE, recursive = TRUE)
if (length(all_h5) == 0) stop("No spatial h5 files found in GSE250163 raw")

sample_ids <- sub("_filtered_feature_bc_matrix\\.h5$", "", basename(all_h5))
# keep first 2 slices for runtime control (can increase later)
sample_ids <- sample_ids[1:min(2, length(sample_ids))]

copy_or_gunzip <- function(src, dst) {
  if (grepl("\\.gz$", src)) {
    con_in <- gzfile(src, "rb")
    con_out <- file(dst, "wb")
    on.exit({close(con_in); close(con_out)}, add = TRUE)
    repeat {
      b <- readBin(con_in, what = "raw", n = 65536)
      if (length(b) == 0) break
      writeBin(b, con_out)
    }
  } else {
    file.copy(src, dst, overwrite = TRUE)
  }
}

spatial_objs <- list()
for (sid in sample_ids) {
  cat("Preparing spatial sample", sid, "\n")
  sdir <- file.path(dir_vis, sid)
  spdir <- file.path(sdir, "spatial")
  dir.create(spdir, recursive = TRUE, showWarnings = FALSE)

  # Source files in GEO raw extract
  src_h5 <- list.files(dir_raw, pattern = paste0("^", sid, "_filtered_feature_bc_matrix\\.h5$"), full.names = TRUE, recursive = TRUE)
  src_pos <- list.files(dir_raw, pattern = paste0("^", sid, "_tissue_positions_list\\.csv\\.gz$"), full.names = TRUE, recursive = TRUE)
  src_sf <- list.files(dir_raw, pattern = paste0("^", sid, "_scalefactors_json\\.json\\.gz$"), full.names = TRUE, recursive = TRUE)
  src_hi <- list.files(dir_raw, pattern = paste0("^", sid, "_tissue_hires_image\\.png\\.gz$"), full.names = TRUE, recursive = TRUE)
  src_lo <- list.files(dir_raw, pattern = paste0("^", sid, "_tissue_lowres_image\\.png\\.gz$"), full.names = TRUE, recursive = TRUE)

  if (length(src_h5) == 0) next

  copy_or_gunzip(src_h5[1], file.path(sdir, "filtered_feature_bc_matrix.h5"))
  if (length(src_pos)) copy_or_gunzip(src_pos[1], file.path(spdir, "tissue_positions_list.csv"))
  if (length(src_sf)) copy_or_gunzip(src_sf[1], file.path(spdir, "scalefactors_json.json"))
  if (length(src_hi)) copy_or_gunzip(src_hi[1], file.path(spdir, "tissue_hires_image.png"))
  if (length(src_lo)) copy_or_gunzip(src_lo[1], file.path(spdir, "tissue_lowres_image.png"))

  so <- Load10X_Spatial(data.dir = sdir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = sid)
  so <- NormalizeData(so)
  spatial_objs[[sid]] <- so
}

if (length(spatial_objs) == 0) stop("No spatial slices loaded")

# signature genes
coef_file <- file.path(ROOT, "results", "best_signature_coefficients.tsv")
coef_dt <- fread(coef_file)
coef_dt[, abs_coef := abs(coef)]
core_genes <- coef_dt[order(-abs_coef)]$gene[1:min(5, nrow(coef_dt))]

for (sid in names(spatial_objs)) {
  so <- spatial_objs[[sid]]
  sig_genes <- unique(coef_dt$gene[coef_dt$gene %in% rownames(so)])
  if (length(sig_genes) >= 3) {
    so <- AddModuleScore(so, features = list(sig_genes), name = "IMRISK")
  }

  # Save spot-level table
  md <- as.data.table(so@meta.data, keep.rownames = "spot")
  fwrite(md, file.path(dir_out, paste0(sid, "_spot_metadata.tsv")), sep = "\t")

  # Spatial plots
  gp <- list()
  if ("IMRISK1" %in% colnames(so@meta.data)) {
    gp[["risk"]] <- SpatialFeaturePlot(so, features = "IMRISK1") + ggtitle(paste0(sid, " immune-metabolic risk"))
  }

  cg <- core_genes[core_genes %in% rownames(so)]
  if (length(cg) > 0) {
    p_gene <- SpatialFeaturePlot(so, features = cg[1:min(4, length(cg))], ncol = 2)
    ggsave(file.path(dir_out, paste0(sid, "_spatial_core_genes.png")), p_gene, width = 12, height = 10, dpi = 180)
  }

  if (length(gp) > 0) {
    ggsave(file.path(dir_out, paste0(sid, "_spatial_risk.png")), gp[[1]], width = 8, height = 7, dpi = 180)
  }

  saveRDS(so, file.path(dir_out, paste0(sid, "_seurat_spatial.rds")))
}

cat("Done spatial analysis for", length(spatial_objs), "slices.\n")
