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


adata = sc.read("D:Scanpy/cd8_elderly_ICB.h5ad")

hs = hotspot.Hotspot(
    adata,
    layer_key="counts",
    model='danb',
    latent_obsm_key="X_scVI",
    umi_counts_obs_key="nCount_RNA"
)

hs.create_knn_graph(weighted_graph=False, n_neighbors=30)
hs_results = hs.compute_autocorrelations()
hs_results.head(20)
hs_genes = hs_results.loc[hs_results.FDR < 0.05].sort_values('Z', ascending=False).head(500).index
local_correlations = hs.compute_local_correlations(hs_genes)

modules = hs.create_modules(
    min_gene_threshold=15, core_only=True, fdr_threshold=0.05
)

module_scores = hs.calculate_module_scores()
hs.plot_local_correlations(vmin=-30, vmax=30)
modules.value_counts()

module = 16
results = hs.results.join(hs.modules)
results = results.loc[results.Module == module]
results.sort_values('Z', ascending=False).head(60)

module_cols = []
for c in module_scores.columns:
    key = f"Module {c}"
    adata.obs[key] = module_scores[c]
    module_cols.append(key)

sc.set_figure_params(dpi=300)
plt.figure(figsize = (25,20))
sc.pl.umap(adata, color = module_cols, frameon=False, cmap = 'viridis', ncols=4, size=1)

results = hs.results.join(hs.modules)
results.to_csv('D:Scanpy/cd8_modules_hotspot_2000hvg.csv')

module_scores = hs.calculate_module_scores()
module_scores.to_csv('D:Scanpy/cd8_modules_scores_hotspot.csv')
