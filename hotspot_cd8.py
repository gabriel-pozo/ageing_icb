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


adata = sc.read("D:Scanpy/cd8_elderly_ICB.h5ad")


# In[3]:


adata


# In[4]:


sc.pl.umap(adata, color = 'celltype_refined', size=1)


# In[5]:


hs = hotspot.Hotspot(
    adata,
    layer_key="counts",
    model='danb',
    latent_obsm_key="X_scVI",
    umi_counts_obs_key="nCount_RNA"
)


# In[6]:


hs.create_knn_graph(weighted_graph=False, n_neighbors=30)


# In[7]:


hs_results = hs.compute_autocorrelations()


# In[16]:


hs_results.head(20)


# In[17]:


hs_genes = hs_results.loc[hs_results.FDR < 0.05].sort_values('Z', ascending=False).head(500).index


# In[18]:


local_correlations = hs.compute_local_correlations(hs_genes)


# In[49]:


modules = hs.create_modules(
    min_gene_threshold=15, core_only=True, fdr_threshold=0.05
)


# In[50]:


module_scores = hs.calculate_module_scores()


# In[51]:


hs.plot_local_correlations(vmin=-30, vmax=30)


# In[53]:


modules.value_counts()


# In[76]:


module = 16

results = hs.results.join(hs.modules)
results = results.loc[results.Module == module]

results.sort_values('Z', ascending=False).head(60)


# In[55]:


module_cols = []
for c in module_scores.columns:
    key = f"Module {c}"
    adata.obs[key] = module_scores[c]
    module_cols.append(key)


# In[40]:


sc.pl.umap(adata, color = ['celltype_refined', 'status', 'ICB_response', 'patient_id'], frameon=False, ncols = 2, size=1)


# In[58]:


sc.set_figure_params(dpi=300)
plt.figure(figsize = (25,20))
sc.pl.umap(adata, color = module_cols, frameon=False, cmap = 'viridis', ncols=3, size=1)


# In[77]:


adata.obs


# In[82]:


results = hs.results.join(hs.modules)
results.to_csv('D:Scanpy/cd8_modules_hotspot_2000hvg.csv')


# In[83]:


module_scores = hs.calculate_module_scores()
module_scores.to_csv('D:Scanpy/cd8_modules_scores_hotspot.csv')


# In[84]:


adata.obs.to_csv("D:Scanpy/artigo_ICB/metadata.csv")


# In[2]:


adata = sc.read("D:Scanpy/cd8_elderly_ICB.h5ad")


# In[3]:


adata.obs


# In[6]:


cell_meta = pd.read_csv("D:Scanpy/full_metadata_ICB.csv", index_col=0)
cell_meta


# In[7]:


adata.obs = cell_meta
adata.obs


# In[10]:


module_cols = ['Module.1', 'Module.2', 'Module.3', 'Module.4', 'Module.5', 'Module.6', 'Module.7', 'Module.8', 'Module.9',
              'Module.10', 'Module.11', 'Module.12', 'Module.13', 'Module.14', 'Module.15', 'Module.16']
sc.set_figure_params(dpi=300)
plt.figure(figsize = (25,20))
sc.pl.umap(adata, color = module_cols, frameon=False, cmap = 'viridis', ncols=4, size=1)


# In[15]:


df = (adata.obs
      .groupby("patient_id")["celltype_refined"]
      .value_counts(normalize=True)
      .mul(100)
      .round(2)
      .unstack())
df


# In[16]:


colors = ['#1f77b4',
'#ff7f0e',
'#2ca02c',
'#d62728',
'#9467bd',
'#8c564b',
'#e377c2',
'#7f7f7f',
'#bcbd22',
'#17becf']
sc.set_figure_params(dpi=300)
fig, ax = plt.subplots(figsize = (6,6))
df.plot.bar(stacked=True, ax=ax, width=0.5, color = colors)
plt.xticks(rotation = 45, rotation_mode = 'anchor', ha = 'right')
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()

