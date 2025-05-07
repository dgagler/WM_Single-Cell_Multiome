# Script to perform SCTransform integration on cluster cuz sucks nuts on local

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(SoupX)
library(celda)
library(DropletUtils)
library(glmGamPoi)

# Load data (need to fix paths)
wm.rna <- readRDS("/gpfs/scratch/gagled01/single-cell/WM_working/scRNA_WM_Merged_QCFiltered.rds")
ctrl.rna <- readRDS("/gpfs/scratch/gagled01/single-cell/WM_working/scRNA_Ctrl_Merged_QCFiltered.rds")

# SCtransform first on ctrls then on WM
ctrl.rna <- SCTransform(ctrl.rna, assay = "RNA", new.assay.name = "SCT", vst.flavor = "v2", ncells = 5000, vars.to.regress = "batch", method = "glmGamPoi") %>%
  RunPCA(npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.7)

wm.rna <- SCTransform(wm.rna, assay = "RNA", new.assay.name = "SCT", vst.flavor = "v2", ncells = 5000, vars.to.regress = "batch", method = "glmGamPoi", verbose = F) %>%
  RunPCA(npcs = 20)

# Integration set up
ifnb.list <- list(ctrl = ctrl.rna, stim = wm.rna)
features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
rm(wm.rna, ctrl.rna)

# Perform integration
anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT", reduction = "rpca",
                                  anchor.features = features, dims = 1:20)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:20)

# Standard workflow
integrated <- RunPCA(integrated, dims = 1:20)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:20)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:20)
integrated <- FindClusters(integrated, resolution = 0.7)
DimPlot(integrated, reduction = "umap")

# Save out
saveRDS(integrated, "/gpfs/scratch/gagled01/single-cell/WM_working/scRNA_SCTransformed_Integrated.rds")
