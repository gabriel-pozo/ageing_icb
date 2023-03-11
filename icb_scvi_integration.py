#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[2]:


# load sparse matrix:
X = io.mmread("D:Scanpy/new_york/counts.mtx")


# In[3]:


# create anndata object
ny = anndata.AnnData(
 X=X.transpose().tocsr()
)


# In[4]:


# load cell metadata:
cell_meta = pd.read_csv("D:Scanpy/new_york/metadata.csv")


# In[5]:


# load gene names:
with open("D:Scanpy/new_york/genes.csv", 'r') as f:
 gene_names = f.read().splitlines()


# In[6]:


# set anndata observations and index obs by barcodes, var by gene names
ny.obs = cell_meta
ny.obs.index = ny.obs['barcode']
ny.var.index = gene_names


# In[7]:


del(cell_meta, gene_names, X)


# In[8]:


# load sparse matrix:
X = io.mmread("D:Scanpy/md_anderson/counts.mtx")


# In[9]:


# create anndata object
md = anndata.AnnData(
 X=X.transpose().tocsr()
)


# In[10]:


# load cell metadata:
cell_meta = pd.read_csv("D:Scanpy/md_anderson/metadata.csv")


# In[11]:


# load gene names:
with open("D:Scanpy/md_anderson/genes.csv", 'r') as f:
 gene_names = f.read().splitlines()


# In[12]:


# set anndata observations and index obs by barcodes, var by gene names
md.obs = cell_meta
md.obs.index = md.obs['barcode']
md.var.index = gene_names


# In[13]:


del(cell_meta, gene_names, X)


# In[14]:


adata = md.concatenate(ny)


# In[15]:


del(md, ny)


# In[16]:


adata.obs.index=adata.obs['barcode']
adata.obs


# In[17]:


def categorise(row):  
    if row['Age'] > 65:
        return 'Elderly'
    elif row['Age'] < 65:
        return 'Non-elderly'


# In[18]:


adata.obs['status'] = adata.obs.apply(lambda row: categorise(row), axis=1)


# In[19]:


adata.layers['counts'] = adata.X.copy()


# In[20]:


sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata


# In[21]:


condition_key = "sample_id"


# In[22]:


sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=3000, layer='counts', subset=True, batch_key=condition_key)


# In[23]:


scvi.model.SCVI.setup_anndata(adata, layer = "counts", batch_key = condition_key,
 continuous_covariate_keys=['nCount_RNA', 'percent.mt'])


# In[24]:


model = scvi.model.SCVI(adata, n_layers=2, n_latent=30)


# In[25]:


model.train(max_epochs=100)


# In[26]:


adata.obsm['X_scVI'] = model.get_latent_representation()


# In[27]:


sc.pp.neighbors(adata, use_rep = 'X_scVI')


# In[28]:


sc.tl.umap(adata)


# In[29]:


sc.pl.umap(adata, color = 'sample_id', frameon=False)


# In[31]:


sc.pl.umap(adata, color = ['status', 'ICB_response', 'Stage', 'patient_id'], frameon=False, ncols=2)


# In[33]:


sc.pl.umap(adata, color = ['CD4', 'CD8A', 'IL7R', 'CCR7', 'SELL', 'IL2RA', 'FOXP3', 'NKG7'], cmap='viridis', ncols=4, frameon=False)


# In[73]:


sc.tl.leiden(adata, resolution=3)


# In[75]:


sc.pl.umap(adata, color='leiden', frameon=False, legend_loc='on data', legend_fontoutline=1, legend_fontsize = 8)


# In[76]:


sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')


# In[77]:


markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
markers


# In[78]:


markers.to_csv("D:Scanpy/marcadores_elderly_419k.csv")


# In[100]:


#markers[markers.group=='18'].head(10)


# In[99]:


#markers[markers.group=='19'].head(10)


# In[98]:


#markers[markers.group=='25'].head(10)


# In[97]:


#markers[markers.group=='28'].head(10)


# In[96]:


#markers[markers.group=='1'].head(10)


# In[95]:


#markers[markers.group=='22'].head(10)


# In[94]:


#markers[markers.group=='16'].head(10)


# In[93]:


#markers[markers.group=='32'].head(10)


# In[86]:


markers[markers.names=='CD4']


# In[87]:


markers[markers.names=='FOXP3']


# In[88]:


markers[markers.names=='CD8A']


# In[89]:


markers[markers.names=='CD8B']


# In[91]:


sc.pl.umap(adata, color = ['leiden', 'CD4', 'CD8A'], frameon=False, groups=('18', '19', '1', '16', '38'), cmap='viridis', ncols=3)


# In[106]:


sc.pl.umap(adata, color = ['leiden', 'CD4', 'CD8A'], frameon=False, groups=('12'), cmap='viridis', ncols=3)


# In[107]:


cell_type = {"0":"CD4-Treg",
"1":"CD8",
"2":"CD8",
"3":"CD4",
"4":"CD4",
"5":"CD8",
"6":"CD8",
"7":"CD4-Treg",
"8":"CD8",
"9":"CD4",
"10":"CD4",
"11":"CD8",
"12":"CD4",
"13":"CD4",
"14":"CD8",
"15":"CD4",
"16":"CD8",
"17":"CD4",
"18":"CD4",
"19":"CD4",
"20":"CD4",
"21":"CD4",
"22":"CD8",
"23":"CD8",
"24":"CD8",
"25":"CD8",
"26":"CD8",
"27":"CD8",
"28":"MAIT",
"29":"CD8",
"30":"CD8",
"31":"CD8",
"32":"CD8",
"33":"CD4",
"34":"CD8",
"35":"CD4",
"36":"CD4",
"37":"CD8",
"38":"CD8",
"39":"CD4"}


# In[108]:


adata.obs['cell_type'] = adata.obs.leiden.map(cell_type)


# In[109]:


sc.set_figure_params(dpi=200)
sc.pl.umap(adata, color = ['cell_type', 'status', 'Stage', 'ICB_response'], frameon=False, ncols=2)


# In[56]:


adata.write("D:Scanpy/artigo_elderly_ICB.h5ad")
model.save("D:Scanpy/elderly_ICB_model.model")


# In[113]:


adata.obs.groupby(['patient_id']).count()


# In[114]:


num_tot_cells = adata.obs.groupby(['patient_id']).count()
num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.leiden))
num_tot_cells


# In[118]:


cell_type_counts = adata.obs.groupby(['patient_id', 'status', 'cell_type']).count()
cell_type_counts = cell_type_counts[cell_type_counts.sum(axis = 1) > 0].reset_index()
cell_type_counts = cell_type_counts[cell_type_counts.columns[0:5]]
cell_type_counts


# In[119]:


cell_type_counts['total_cells'] = cell_type_counts.patient_id.map(num_tot_cells).astype(int)
cell_type_counts['frequency'] = cell_type_counts.nCount_RNA / cell_type_counts.total_cells
cell_type_counts


# In[124]:


sc.set_figure_params(dpi=100)
import matplotlib.pyplot as plt
plt.figure(figsize = (5,5))
ax = sns.boxplot(data = cell_type_counts, x = 'cell_type', y = 'frequency', hue = 'status')
plt.xticks(rotation = 35, rotation_mode = 'anchor', ha = 'right')
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()


# In[2]:


cell_meta1 = pd.read_csv("D:Scanpy/cell_meta1_elderly.csv", index_col=0)
cell_meta2 = pd.read_csv("D:Scanpy/cell_meta2_elderly.csv", index_col=0)
cell_meta3 = pd.read_csv("D:Scanpy/cell_meta3_elderly.csv", index_col=0)
df = pd.concat([cell_meta1, cell_meta2, cell_meta3])
df


# In[3]:


del(df['barcode.1'])


# In[4]:


df


# In[5]:


adata = sc.read("D:Scanpy/artigo_elderly_ICB.h5ad")
adata.obs


# In[13]:


adata[adata.obs.cell_type=='MAIT']


# In[12]:


df.groupby(['celltype_refined']).count()


# In[14]:


adata.obs['celltype_refined'] = df.celltype_refined


# In[16]:


sc.set_figure_params(dpi=200)
sc.pl.umap(adata, color = ['cell_type', 'celltype_refined'], ncols=2, frameon=False)


# In[17]:


adata.obs.to_csv("D:Scanpy/metadata_elderly_final.csv")


# In[18]:


del(df, cell_meta1, cell_meta2, cell_meta3)


# In[19]:


adata.obs.groupby(['patient_id']).count()


# In[16]:


num_tot_cells = adata.obs.groupby(['patient_id']).count()
num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.leiden))
num_tot_cells


# In[21]:


cell_type_counts = adata.obs.groupby(['patient_id', 'status', 'celltype_refined']).count()
cell_type_counts = cell_type_counts[cell_type_counts.sum(axis = 1) > 0].reset_index()
cell_type_counts = cell_type_counts[cell_type_counts.columns[0:5]]
cell_type_counts


# In[22]:


cell_type_counts['total_cells'] = cell_type_counts.patient_id.map(num_tot_cells).astype(int)
cell_type_counts['frequency'] = cell_type_counts.nCount_RNA / cell_type_counts.total_cells
cell_type_counts


# In[23]:


import matplotlib.pyplot as plt
plt.figure(figsize = (12,5))
ax = sns.boxplot(data = cell_type_counts, x = 'celltype_refined', y = 'frequency', hue = 'status')
plt.xticks(rotation = 45, rotation_mode = 'anchor', ha = 'right')
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()


# In[29]:


cell_type_counts = adata.obs.groupby(['patient_id', 'ICB_response', 'celltype_refined']).count()
cell_type_counts = cell_type_counts[cell_type_counts.sum(axis = 1) > 0].reset_index()
cell_type_counts = cell_type_counts[cell_type_counts.columns[0:5]]
cell_type_counts


# In[30]:


cell_type_counts['total_cells'] = cell_type_counts.patient_id.map(num_tot_cells).astype(int)
cell_type_counts['frequency'] = cell_type_counts.nCount_RNA / cell_type_counts.total_cells
cell_type_counts


# In[31]:


import matplotlib.pyplot as plt
plt.figure(figsize = (12,5))
ax = sns.boxplot(data = cell_type_counts, x = 'celltype_refined', y = 'frequency', hue = 'ICB_response')
plt.xticks(rotation = 45, rotation_mode = 'anchor', ha = 'right')
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()


# In[2]:


adata = sc.read("D:Scanpy/artigo_elderly_ICB_final.h5ad")


# In[3]:


ax = sc.pl.correlation_matrix(adata, 'celltype_refined')


# In[3]:


sc.set_figure_params(dpi=300)
plt.figure(figsize = (24,15))
sc.pl.umap(adata, color = ['cell_type', 'status', 'ICB_response', 'patient_id'], frameon=False, ncols = 2, wspace=0.3,
          show=False, size=1)


# In[4]:


num_total_cells = adata.obs.groupby(['patient_id']).count()
num_total_cells = dict(zip(num_total_cells.index, num_total_cells.leiden))
num_total_cells


# In[10]:


cell_type_counts = adata.obs.groupby(['patient_id', 'status', 'celltype_refined']).count()
cell_type_counts = cell_type_counts[cell_type_counts.sum(axis = 1) > 0].reset_index()
cell_type_counts = cell_type_counts[cell_type_counts.columns[0:5]]
cell_type_counts


# In[11]:


cell_type_counts['total_cells'] = cell_type_counts.patient_id.map(num_total_cells).astype(int)
cell_type_counts['frequency'] = cell_type_counts.nCount_RNA / cell_type_counts.total_cells
cell_type_counts


# In[12]:


sc.set_figure_params(dpi=300)
plt.figure(figsize = (12,5))
ax = sns.boxplot(data = cell_type_counts, x = 'celltype_refined', y = 'frequency', hue = 'status')
plt.xticks(rotation = 45, rotation_mode = 'anchor', ha = 'right')
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()


# In[3]:


df = (adata.obs
      .groupby("status")["cell_type"]
      .value_counts(normalize=True)
      .mul(100)
      .round(2)
      .unstack())
df


# In[6]:


colors = ['#2ca02c', '#1f77b4', '#ff7f0e', '#d62728']
sc.set_figure_params(dpi=300)
fig, ax = plt.subplots(figsize = (6,6))
df.plot.bar(stacked=True, ax=ax, width=0.5, color = colors)
plt.xticks(rotation = 45, rotation_mode = 'anchor', ha = 'right')
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()


# In[22]:


sc.set_figure_params(dpi=300)
plt.figure(figsize = (24,15))
sc.pl.umap(adata, color = ['CD8A', 'CD4', 'FOXP3', 'TRBV18'], frameon=False, ncols = 2,
          show=False, size=1, cmap='viridis')


# In[2]:


adata = sc.read("D:Scanpy/artigo_elderly_ICB_final.h5ad")


# In[3]:


sc.set_figure_params(dpi=200)
plt.figure(figsize = (6,6))
marker_genes_dict = {
    'CD4': ['SELL', 'CCR7', 'KLF2', 'TCF7'],
    'CD4-Treg': ['IL2RA', 'TNFRSF4', 'TNFRSF18', 'BATF'],
    'CD8': ['GZMK', 'GZMA', 'GZMB', 'NKG7'],
    'MAIT': ['TRAV13-1', 'KLRB1', 'CCR6', 'IL4I1'],
}
sc.pl.dotplot(adata, marker_genes_dict, 'cell_type', dendrogram=True)

