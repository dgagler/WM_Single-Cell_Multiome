#!/usr/bin/env Rscript

# R script to run rPCA integration on HPC

# Load libraries
library(Seurat, lib.loc = "/gpfs/share/apps/R/4.1.2/lib64/R/library")
library(ggplot2, lib.loc = "/gpfs/share/apps/R/4.1.2/lib64/R/library")
library(dplyr, lib.loc = "/gpfs/share/apps/R/4.1.2/lib64/R/library")

# Load merged QC filtered data object
merged.obj <- readRDS("/gpfs/scratch/gagled01/single-cell/WM_working/scRNA_WM_Merged_Unfiltered_CellNamesFixed.rds")

# Split into list
lib.list <- SplitObject(merged.obj, split.by = "patient")

# Normalize and find variable features for each lib
lib.list <- lapply(X = lib.list, FUN = function(x) {
  x <- NormalizeData(x) 
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = lib.list)
lib.list <- lapply(X = lib.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Find anchors
lib.anchors <- FindIntegrationAnchors(object.list = lib.list, anchor.features = features, reduction = "rpca")

# Integrate data
lib.combined <- IntegrateData(anchorset = lib.anchors)

# Save out
saveRDS(lib.combined, "/gpfs/scratch/gagled01/single-cell/WM_working/5TGM1_IntegratedObject_StandardWorkflow.rds")