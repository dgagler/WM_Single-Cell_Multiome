library(Seurat)
library(Matrix)

cell_subset <- readRDS("/Users/gagled01/morganLab/Waldenstroms/singlecell/objects/scRNA_SCTransformed_Integrated.rds")

### CONVERTING TO MATRICES/LISTS FOR CELL TYPIST INPUT ###
# Read in Seurat object

# CellTypist requires a log-normalized (to 10,000) expression matrix, a list of cells, and a list of genes/features
raw.matrix <- cell_subset@assays$RNA@data 
cell.names <- Cells(cell_subset@assays$RNA)
gene.names <- rownames(cell_subset@assays$RNA)

# Write out sparse matrix
writeMM(raw.matrix, '/Users/gagled01/morganLab/Waldenstroms/singlecell/scRNA_SCTransformIntegrated_Matrix.mtx')

# Write out gene and cell lists
write.table(cell.names, '/Users/gagled01/morganLab/Waldenstroms/singlecell/scRNA_SCTransformIntegrated_cellNames.tsv', quote = F, sep = '\t', 
            col.names = F, row.names = F)
write.table(gene.names, '/Users/gagled01/morganLab/Waldenstroms/singlecell/scRNA_SCTransformIntegrated_geneNames.tsv', quote = F, sep = '\t',
            col.names = F, row.names = F)

length(cell.names)
length(gene.names)
dim(raw.matrix)
