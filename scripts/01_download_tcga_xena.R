#!/usr/bin/env Rscript

# Download TCGA-KIRC data using UCSCXenaTools helper wrappers.
# This avoids hard-coding raw dataset filenames (which often change).

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

suppressPackageStartupMessages(library(UCSCXenaTools))

out_dir <- file.path(ROOT, "data", "TCGA_Xena")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Downloading TCGA-KIRC from Xena into: ", out_dir)

# 1) Clinical/phenotype (contains OS fields)
downloadTCGA(
  project = "KIRC",
  data_type = "Phenotype",
  file_type = "Clinical Information",
  destdir = out_dir,
  force = FALSE
)

# 2) Gene expression (RNASeq). We prefer the standard IlluminaHiSeq RNASeqV2 matrix
# which is widely used in legacy TCGA analyses on Xena.
# Note: this is not HTSeq counts; but it's sufficient for signature training.
downloadTCGA(
  project = "KIRC",
  data_type = "Gene Expression RNASeq",
  file_type = "IlluminaHiSeq RNASeqV2",
  destdir = out_dir,
  force = FALSE
)

message("Done.")
