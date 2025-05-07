#!/usr/bin/env Rscript

# Read file
seurat <- readRDS("/gpfs/scratch/gagled01/single-cell/WM_working/scRNA_WM_NoControls_rPCA_SeuratIntegrated_CTAnnotated.rds")

# DefaultAssay(seurat.sub) <- "RNA"
seurat <- NormalizeData(seurat)
seurat <- ScaleData(seurat)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)

saveRDS(seurat, "/gpfs/scratch/gagled01/single-cell/WM_working/scRNA_WM_NoControls_rPCA_SeuratIntegrated_CTAnnotated_NormScaleFindVar.rds")