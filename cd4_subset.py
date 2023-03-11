#!/usr/bin/env python
# coding: utf-8


#Load the libraries

import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import scanpy.external as sce
import seaborn as sns
import anndata as ad
import squidpy as sq
import scvi
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz
import warnings
warnings.filterwarnings("ignore")
from sccoda.util import comp_ana as mod
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import os
import anndata
import hotspot

adata = sc.read("D:Scanpy/raw_artigo_elderly_ICB.h5ad")
adata.obs
adata[adata.obs.cell_type=='CD4']
adata = adata[adata.obs.cell_type=='CD4']


adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

condition_key = "sample_id"
sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=2000, layer='counts', subset=True, batch_key=condition_key)

scvi.model.SCVI.setup_anndata(adata, layer = "counts", batch_key = condition_key,
 continuous_covariate_keys=['nCount_RNA', 'percent.mt'])

model = scvi.model.SCVI(adata, n_layers=2, n_latent=30)
model.train(max_epochs = 100)

adata.obsm['X_scVI'] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep = 'X_scVI')
sc.tl.umap(adata)
sc.pl.umap(adata, color = 'sample_id', frameon=False)
sc.pl.umap(adata, color = ['status', 'ICB_response', 'Stage', 'patient_id'], frameon=False, ncols=2)
sc.pl.umap(adata, color = ['IL7R', 'CCR7', 'SELL', 'TCF7'], frameon=False, ncols=2, cmap='viridis')
sc.pl.umap(adata, color = ['GZMA', 'CCL5', 'CXCL13', 'NKG7'], frameon=False, ncols=2, cmap='viridis')
sc.pl.umap(adata, color = ['CTLA4', 'PDCD1', 'HAVCR2', 'TIGIT'], frameon=False, ncols=2, cmap='viridis')
sc.pl.umap(adata, color = ['ITGAE', 'CXCR6', 'GZMK', 'KLRB1'], frameon=False, ncols=2, cmap='viridis')
sc.pl.umap(adata, color = ['HSPA1A', 'TGFB1', 'HLA-DPA1', 'CCR6'], frameon=False, ncols=2, cmap='viridis')

sc.tl.leiden(adata, resolution = 1.5)
sc.pl.umap(adata, color = 'leiden', frameon=False)
sc.pl.umap(adata, color = 'leiden', frameon=False, legend_loc='on data')

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
markers
markers.to_csv("D:Scanpy/cd4_markers_elderly_ICB.csv")

celltype_refined = {
"0":"CD4-HSPA1A+",
"1":"CD4-TGFB1+",
"2":"CD4-GZMK+",
"3":"CD4-CXCL13+",
"4":"CD4-HLA(II)",
"5":"CD4-naive",
"6":"CD4-CXCL13+",
"7":"CD4-HSPA1A+",
"8":"CD4-CCR6+",
"9":"CD4-GZMK+",
"10":"CD4-CCR6+",
"11":"CD4-TIGIT+",
"12":"CD4-ITGAE+",
"13":"CD4-CCR6+",
"14":"CD4-HLA(II)",
"15":"CD4-GZMK+"
}

adata.obs['celltype_refined'] = adata.obs.leiden.map(celltype_refined)

sc.pl.umap(adata, color = ['celltype_refined', 'status', 'ICB_response', 'patient_id'], ncols=2, frameon=False)

adata.write("D:Scanpy/cd4_elderly_ICB.h5ad")
model.save("D:Scanpy/model_cd4_elderly_ICB.model")


adata.obs.groupby(['patient_id']).count()
num_total_cells = adata.obs.groupby(['patient_id']).count()
num_total_cells = dict(zip(num_total_cells.index, num_total_cells.leiden))
num_total_cells
cell_type_counts = adata.obs.groupby(['patient_id', 'status', 'celltype_refined']).count()
cell_type_counts = cell_type_counts[cell_type_counts.sum(axis= 1 ) > 0].reset_index()
cell_type_counts = cell_type_counts[cell_type_counts.columns[0:5]]
cell_type_counts
cell_type_counts['total_cells'] = cell_type_counts.patient_id.map(num_total_cells).astype(int)
cell_type_counts['frequency'] = cell_type_counts.nCount_RNA / cell_type_counts.total_cells
cell_type_counts

import matplotlib.pyplot as plt
plt.figure(figsize = (10,4))
ax = sns.boxplot(data = cell_type_counts, x = 'celltype_refined', y = 'frequency', hue = 'status')
plt.xticks(rotation = 35, rotation_mode = 'anchor', ha = 'right')
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()


cov_df = pd.DataFrame({"status": ["Elderly", "Elderly", "Elderly", "Elderly", "Elderly", "Elderly", "Elderly", "Elderly",
                                 "Elderly", "Elderly", "Elderly", "Non-elderly", "Non-elderly", "Non-elderly", "Non-elderly",
                                 "Non-elderly"]}, index = ["MD01004", "MD043008", "MD01019", "MD043006", "MD01024", "MD01010",
                                                          "NY016007", "NY016016", "NY016021", "NY016022", "NY016025", "MD01005",
                                                          "MD043011", "MD043003", "NY016014", "NY016015"])


data_scanpy_1 = dat.from_scanpy(
    adata,
    cell_type_identifier="celltype_refined",
    sample_identifier="patient_id",
    covariate_df=cov_df
)
print(data_scanpy_1)

# Stacked barplot for each sample
viz.stacked_barplot(data_scanpy_1, feature_name="samples")
plt.xticks(rotation=90, fontsize=10)
plt.show()

# Stacked barplot for the levels of "Condition"
viz.stacked_barplot(data_scanpy_1, feature_name="status")
plt.xticks(rotation=0)
plt.show()


# Grouped boxplots. No facets, relative abundance, no dots.
viz.boxplots(
    data_scanpy_1,
    feature_name="status",
    plot_facets=False,
    y_scale="relative",
    add_dots=False,
    cmap={"Non-elderly": "darkorange", "Elderly": "royalblue"},
)
plt.show()

# Grouped boxplots. Facets, log scale, added dots and custom color palette.
viz.boxplots(
    data_scanpy_1,
    feature_name="status",
    plot_facets=True,
    y_scale="log",
    add_dots=True,
    cmap={"Non-elderly": "darkorange", "Elderly": "royalblue"},
)
plt.show()


sc.set_figure_params(dpi=300)
plt.figure(figsize = (24,15))
sc.pl.umap(adata, color = ['celltype_refined', 'status', 'ICB_response', 'patient_id'], frameon=False, ncols = 2, wspace=0.4,
          show=False, size=1)


sc.tl.rank_genes_groups(adata, 'celltype_refined', method = 'wilcoxon')

markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
markers

markers[markers.group=='CD4-TGFB1+'].head(20)
marker_genes_dict = {
    'CD4-CCR6+': ['CCR6', 'KLRB1'],
    'CD4-CXCL13+': ['CXCL13', 'CTLA4'],
    'CD4-GZMK+': ['GZMK', 'NKG7'],
    'CD4-HLA(II)': ['LGALS1', 'HLA-DRB1'],
    'CD4-HSPA1A+': ['HSPA1A', 'HSPA1B'],
    'CD4-ITGAE+': ['ITGAE', 'GZMB'],
    'CD4-TGFB1+': ['TGFB1', 'CDKN1A'],
    'CD4-TIGIT+':['TIGIT', 'IKZF2'],
    'CD4-naive':['SELL', 'CCR7']
}
sc.pl.dotplot(adata, marker_genes_dict, 'celltype_refined', dendrogram=True)
