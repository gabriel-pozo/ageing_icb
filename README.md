# ageing_icb
All code for reanalysis of Caushi dataset (GSE173351)

seurat_prepare.R --> Process and filtering with Seurat
scvi_integration.ipynb --> Import filtered data and run integration, batch correction and clustering
cd4_subset --> Subset and recluster CD4+ T cells in order to refine characterization
cd8_subset --> Subset and recluster CD8+ T cells in order to refine characterization
hotspot_cd8 --> Identify gene modules related to CD8+ T cells
pseudobulk_senescence_pca --> Convert anndata object into singlecellexperiment, build pseudobulk object with DESeq2 and plot PCA
