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


# In[21]:


metadata = pd.read_csv("D:Scanpy/metadata_elderly.csv", index_col=0)


# In[22]:


metadata


# In[24]:


b = metadata['barcode.1']


# In[25]:


adata[adata.obs.barcode.isin(b)]


# In[26]:


adata = adata[adata.obs.barcode.isin(b)]


# In[27]:


adata.shape


# In[28]:


adata.obs['cell_type'] = metadata['cell_type']


# In[29]:


adata.obs


# In[30]:


adata.write("D:Scanpy/raw_artigo_elderly_ICB.h5ad")


# In[ ]:





# In[2]:


adata = sc.read("D:Scanpy/raw_artigo_elderly_ICB.h5ad")


# In[3]:


adata[adata.obs.cell_type=='CD8']


# In[4]:


adata = adata[adata.obs.cell_type=='CD8']


# In[5]:


adata.layers['counts'] = adata.X.copy()


# In[6]:


sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata


# In[7]:


condition_key = "sample_id"


# In[8]:


sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=2000, layer='counts', subset=True, batch_key=condition_key)


# In[9]:


scvi.model.SCVI.setup_anndata(adata, layer = "counts", batch_key = condition_key,
 continuous_covariate_keys=['nCount_RNA', 'percent.mt'])


# In[10]:


model = scvi.model.SCVI(adata, n_layers=2, n_latent=30)


# In[11]:


model.train(max_epochs=100)


# In[12]:


adata.obsm['X_scVI'] = model.get_latent_representation()


# In[13]:


sc.pp.neighbors(adata, use_rep = 'X_scVI')


# In[14]:


sc.tl.umap(adata)


# In[15]:


sc.pl.umap(adata, color = 'sample_id', frameon=False)


# In[16]:


sc.pl.umap(adata, color = ['status', 'ICB_response', 'Stage', 'patient_id'], frameon=False, ncols=2)


# In[17]:


sc.tl.leiden(adata, resolution = 1)


# In[18]:


sc.pl.umap(adata, color = 'leiden', frameon=False)


# In[19]:


sc.pl.umap(adata, color = ['FGFBP2', 'ITGAE', 'CXCR6', 'KLRB1'], ncols=2, cmap='viridis', frameon=False)


# In[20]:


sc.pl.umap(adata, color = ['PDCD1', 'HAVCR2', 'LAG3', 'CTLA4'], ncols=2, cmap='viridis', frameon=False)


# In[21]:


sc.pl.umap(adata, color = ['GZMK', 'TIGIT', 'ZNF683', 'GNLY'], ncols=2, cmap='viridis', frameon=False)


# In[24]:


sc.pl.umap(adata, color = ['ENTPD1', 'NT5E', 'EOMES', 'TNFRSF9'], ncols=2, cmap='viridis', frameon=False)


# In[32]:


sc.pl.umap(adata, color = ['HLA-DRA', 'CD74', 'CXCL13', 'MKI67'], ncols=2, cmap='viridis', frameon=False)


# In[25]:


sc.pl.umap(adata, color = 'leiden', frameon=False, legend_loc='on data')


# In[59]:


sc.pl.umap(adata, color = 'leiden', frameon=False, groups=('8'),)


# In[41]:


sc.pl.umap(adata, color = ['CD8A', 'CD8B', 'CD4'], frameon=False, cmap='viridis')


# In[43]:


sc.pl.umap(adata, color = ['GZMA', 'GZMB', 'GZMH', 'PRF1'], frameon=False, cmap='viridis', ncols=2)


# In[44]:


sc.pl.umap(adata, color = ['TCF7', 'IL7R', 'CCR7', 'SELL'], frameon=False, cmap='viridis', ncols=2)


# In[58]:


sc.pl.umap(adata, color = ['TOX', 'CD27', 'TIGIT', 'HAVCR2'], frameon=False, cmap='viridis', ncols=2)


# In[26]:


sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')


# In[27]:


markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
markers


# In[64]:


markers.to_csv("D:Scanpy/markers_cd8_elderly_ICB.csv")


# In[83]:


adata.write("D:Scanpy/cd8_elderly_ICB.h5ad")
model.save("D:Scanpy/model_cd8_elderly_ICB.model")


# In[2]:


adata = sc.read("D:Scanpy/cd8_elderly_ICB.h5ad")


# In[3]:


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


# In[4]:


adata.obs['celltype_refined'] = adata.obs.leiden.map(celltype_refined)


# In[5]:


sc.pl.umap(adata, color = ['celltype_refined', 'status', 'ICB_response', 'patient_id'], frameon=False, ncols = 2, wspace=0.2)


# In[6]:


adata.obs.groupby(['patient_id']).count()


# In[7]:


num_tot_cells = adata.obs.groupby(['patient_id']).count()
num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.leiden))
num_tot_cells


# In[8]:


cell_type_counts = adata.obs.groupby(['patient_id', 'status', 'celltype_refined']).count()
cell_type_counts = cell_type_counts[cell_type_counts.sum(axis = 1) > 0].reset_index()
cell_type_counts = cell_type_counts[cell_type_counts.columns[0:5]]
cell_type_counts


# In[9]:


cell_type_counts['total_cells'] = cell_type_counts.patient_id.map(num_tot_cells).astype(int)
cell_type_counts['frequency'] = cell_type_counts.nCount_RNA / cell_type_counts.total_cells
cell_type_counts


# In[11]:


plt.figure(figsize = (10,4))
ax = sns.boxplot(data = cell_type_counts, x = 'celltype_refined', y = 'frequency', hue = 'status')
plt.xticks(rotation = 35, rotation_mode = 'anchor', ha = 'right')
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()


# In[12]:


adata.write("D:Scanpy/cd8_elderly_ICB.h5ad")


# In[14]:


cov_df = pd.DataFrame({"status": ["Elderly", "Elderly", "Elderly", "Elderly", "Elderly", "Elderly", "Elderly", "Elderly",
                                 "Elderly", "Elderly", "Elderly", "Non-elderly", "Non-elderly", "Non-elderly", "Non-elderly",
                                 "Non-elderly"]}, index = ["MD01004", "MD043008", "MD01019", "MD043006", "MD01024", "MD01010",
                                                          "NY016007", "NY016016", "NY016021", "NY016022", "NY016025", "MD01005",
                                                          "MD043011", "MD043003", "NY016014", "NY016015"])


# In[19]:


cov_df


# In[15]:


data_scanpy_1 = dat.from_scanpy(
    adata,
    cell_type_identifier="celltype_refined",
    sample_identifier="patient_id",
    covariate_df=cov_df
)
print(data_scanpy_1)


# In[17]:


# Stacked barplot for each sample
viz.stacked_barplot(data_scanpy_1, feature_name="samples")
plt.xticks(rotation=90, fontsize=10)
plt.show()

# Stacked barplot for the levels of "Condition"
viz.stacked_barplot(data_scanpy_1, feature_name="status")
plt.xticks(rotation=0)
plt.show()


# In[18]:


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


# In[68]:


sc.set_figure_params(dpi=300)
plt.figure(figsize = (25,25))
sc.pl.umap(adata, color = ['CXCL13', 'ITGAE', 'CXCR6', 'ZNF683', 'ENTPD1', 'CTLA4', 'PDCD1', 'HAVCR2', 'LAG3', 'TIGIT', 'TOX', 'LAYN'],
          frameon=False, cmap='viridis', ncols=3,show=False, size=1, wspace=0.2)


# In[66]:


sc.set_figure_params(dpi=300)
plt.figure(figsize = (24,15))
sc.pl.umap(adata, color = ['celltype_refined', 'status', 'ICB_response', 'patient_id'], frameon=False, ncols = 2, wspace=0.3,
          show=False, size=1)


# In[11]:


sc.tl.rank_genes_groups(adata, 'celltype_refined', method='wilcoxon')


# In[12]:


markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
markers


# In[77]:


markers[markers.group=='CD8-IL7R+'].head(30)


# In[86]:


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

