# A multiomic analysis of Waldenström macroglobulinemia defines distinct disease subtypes
> This repository contains code used for the bioinformatic analysis of the scRNA-seq and scATAC-seq data in Gagler et al., 2025, Blood

The raw data for this study can be found on the NCBI's Gene Expression Omnibus (GEO) database under accession [GSE296167](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296167). The code included in this repository cover the major analytical workflow of the study, outlined as follows:

## Analytical workflow:

1. **scRNA_Analysis.Rmd** - scRNA preprocessing, QC, integration, and creation of input files for CellTypist
2. **CellTypist_Workflow.ipynb** - Running CellTypist in python
3. **PostCellTypist_h5_to_Seurat_Conversion.Rmd** - Converting CellTypist outputs back into Seurat objects
4. **scATAC_Analysis_PreprocessingQC_Clustering_scRNA_Integration.Rmd** - scATAC processing and integration with scRNA
5. **PostIntegration_AddingPeaksDeviations.Rmd** - Calling peaks and chromVAR deviations on integrated data
6. **run_TRUST4.sh** - Running TRUST4 on each patient
7. **TRUST4_Annotation.Rmd** - Adding TRUST4 annotations to integrated object

More details regarding analysis can be found in the comments of the associated R markdown scripts. In addition to many R packages including Seurat and ArchR, this analysis requires the installation and execution of [MACS2](https://github.com/macs3-project/MACS/wiki/Install-macs2) and [TRUST4](https://github.com/liulab-dfci/TRUST4). 

For questions, please reach out to dylangagler@gmail.com or gareth.morgan@nyulangone.org.
