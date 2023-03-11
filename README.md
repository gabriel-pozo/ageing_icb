# ageing_icb
All code for reanalysis of Caushi dataset (GSE173351)

seurat_prepare.R --> Process and filtering with Seurat

icb_scvi_integration.py --> Import filtered data and run integration, batch correction and clustering

cd4_subset.py --> Subset and recluster CD4+ T cells in order to refine characterization

cd8_subset.py --> Subset and recluster CD8+ T cells in order to refine characterization

hotspot_cd8.py --> Identify gene modules related to CD8+ T cells

pseudobulk_senescence_pca.R --> Convert anndata object into singlecellexperiment, build pseudobulk object with DESeq2 and plot PCA
