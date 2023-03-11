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
adata[adata.obs.cell_type=='CD8']
adata = adata[adata.obs.cell_type=='CD8']

adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

condition_key = "sample_id"
sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=2000, layer='counts', subset=True, batch_key=condition_key)

scvi.model.SCVI.setup_anndata(adata, layer = "counts", batch_key = condition_key,
 continuous_covariate_keys=['nCount_RNA', 'percent.mt'])

model = scvi.model.SCVI(adata, n_layers=2, n_latent=30)
model.train(max_epochs=100)

adata.obsm['X_scVI'] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep = 'X_scVI')
sc.tl.umap(adata)
sc.pl.umap(adata, color = 'sample_id', frameon=False)
sc.pl.umap(adata, color = ['status', 'ICB_response', 'Stage', 'patient_id'], frameon=False, ncols=2)

sc.tl.leiden(adata, resolution = 1)
sc.pl.umap(adata, color = 'leiden', frameon=False)
sc.pl.umap(adata, color = ['FGFBP2', 'ITGAE', 'CXCR6', 'KLRB1'], ncols=2, cmap='viridis', frameon=False)
sc.pl.umap(adata, color = ['PDCD1', 'HAVCR2', 'LAG3', 'CTLA4'], ncols=2, cmap='viridis', frameon=False)
sc.pl.umap(adata, color = ['GZMK', 'TIGIT', 'ZNF683', 'GNLY'], ncols=2, cmap='viridis', frameon=False)
sc.pl.umap(adata, color = ['ENTPD1', 'NT5E', 'EOMES', 'TNFRSF9'], ncols=2, cmap='viridis', frameon=False)
sc.pl.umap(adata, color = ['HLA-DRA', 'CD74', 'CXCL13', 'MKI67'], ncols=2, cmap='viridis', frameon=False)

sc.pl.umap(adata, color = 'leiden', frameon=False, legend_loc='on data')
sc.pl.umap(adata, color = 'leiden', frameon=False, groups=('8'),)
sc.pl.umap(adata, color = ['CD8A', 'CD8B', 'CD4'], frameon=False, cmap='viridis')
sc.pl.umap(adata, color = ['GZMA', 'GZMB', 'GZMH', 'PRF1'], frameon=False, cmap='viridis', ncols=2)
sc.pl.umap(adata, color = ['TCF7', 'IL7R', 'CCR7', 'SELL'], frameon=False, cmap='viridis', ncols=2)
sc.pl.umap(adata, color = ['TOX', 'CD27', 'TIGIT', 'HAVCR2'], frameon=False, cmap='viridis', ncols=2)

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
markers
markers.to_csv("D:Scanpy/markers_cd8_elderly_ICB.csv")

celltype_refined = {"0":"CD8-memory",
"1":"CD8-effector",
"2":"CD8-effector",
"3":"CD8-effector",
"4":"CD8-IL7R+",
"5":"CD8-CXCL13+",
"6":"CD8-memory",
"7":"DP-CD4+/CD8+",
"8":"CD8-CXCL13+",
"9":"CD8-KLRB1+",
"10":"CD8-effector",
"11":"CD8-cycling",
"12":"CD8-AREG+",
"13":"CD8-effector",
"14":"CD8-memory",
"15":"CD8-CXCL13+"
}
adata.obs['celltype_refined'] = adata.obs.leiden.map(celltype_refined)
sc.pl.umap(adata, color = ['celltype_refined', 'status', 'ICB_response', 'patient_id'], frameon=False, ncols = 2, wspace=0.2)

adata.write("D:Scanpy/cd8_elderly_ICB.h5ad")
model.save("D:Scanpy/model_cd8_elderly_ICB.model")


adata.obs.groupby(['patient_id']).count()
num_tot_cells = adata.obs.groupby(['patient_id']).count()
num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.leiden))
num_tot_cells
cell_type_counts = adata.obs.groupby(['patient_id', 'status', 'celltype_refined']).count()
cell_type_counts = cell_type_counts[cell_type_counts.sum(axis = 1) > 0].reset_index()
cell_type_counts = cell_type_counts[cell_type_counts.columns[0:5]]
cell_type_counts
cell_type_counts['total_cells'] = cell_type_counts.patient_id.map(num_tot_cells).astype(int)
cell_type_counts['frequency'] = cell_type_counts.nCount_RNA / cell_type_counts.total_cells
cell_type_counts
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
plt.figure(figsize = (25,25))
sc.pl.umap(adata, color = ['CXCL13', 'ITGAE', 'CXCR6', 'ZNF683', 'ENTPD1', 'CTLA4', 'PDCD1', 'HAVCR2', 'LAG3', 'TIGIT', 'TOX', 'LAYN'],
          frameon=False, cmap='viridis', ncols=3,show=False, size=1, wspace=0.2)

sc.set_figure_params(dpi=300)
plt.figure(figsize = (24,15))
sc.pl.umap(adata, color = ['celltype_refined', 'status', 'ICB_response', 'patient_id'], frameon=False, ncols = 2, wspace=0.3,
          show=False, size=1)


sc.tl.rank_genes_groups(adata, 'celltype_refined', method='wilcoxon')
markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
markers
markers[markers.group=='CD8-IL7R+'].head(30)

markers_gene_dict = {
'CD8-AREG+': ['AREG', 'KLRC3', 'XCL1'],
    'CD8-CXCL13+': ['CXCL13', 'CTLA4', 'PDCD1'],
    'CD8-IL7R+': ['IL7R', 'CD55', 'GPR183'],
    'CD8-KLRB1+': ['KLRB1', 'CCR6', 'NCR3'],
    'CD8-memory': ['ITGAE', 'CXCR6', 'ZNF683'],
    'CD8-effector': ['GZMK', 'CMC1', 'CST7'],
    'CD8-cycling':['MKI67', 'STMN1', 'TIGIT'],
    'DP-CD4+/CD8+':['CD4', 'CD8A', 'BATF'],    
}
sc.set_figure_params(dpi = 200)
sc.pl.dotplot(adata, markers_gene_dict, 'celltype_refined', dendrogram=True)
