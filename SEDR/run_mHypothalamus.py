import scanpy as sc
import pandas as pd
from sklearn import metrics
import torch

import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA  # sklearn PCA is used because PCA in scanpy is not stable. 

import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.append('/home/lytq/Spatial-Transcriptomics-Benchmark/utils')
from evaluate import evaluate_clustering
from load_st_data import load_mHypothalamus
from res_search import search_resolution

import time
import tracemalloc

import SEDR

data_path = '/home/lytq/Spatial-Transcriptomics-Benchmark/data/mHypothalamus'
output_path = '/home/lytq/Spatial-Transcriptomics-Benchmark/Results/MERFISH/mHypothalamus/SEDR'



device = 'cuda:5' if torch.cuda.is_available() else 'cpu'
for section_id in data_names:
    print(f"================ Start Processing {section_id} ======================")
    random_seed = 2025
    SEDR.fix_seed(random_seed)
    
    n_clusters = 8
    
    dir_out = f'{output_path}/{section_id}'
    os.makedirs(dir_out, exist_ok=True)
    
    tracemalloc.start()
    start_time = time.time()
    
    # Load data
    adata = load_mHypothalamus(data_path, section_id)
    df_meta = adata.obs

    sc.pp.filter_genes(adata, min_cells=50)
    sc.pp.filter_genes(adata, min_counts=10)
    sc.pp.normalize_total(adata, target_sum=1, exclude_highly_expressed=True, inplace=False)['X']
    sc.pp.scale(adata)

    n_comps = min(200, len(adata.var.index) - 1)
    adata_X = PCA(n_components=n_comps, random_state=42).fit_transform(adata.X)
    adata.obsm['X_pca'] = adata_X

    graph_dict = SEDR.graph_construction(adata, 12)
    
    sedr_net = SEDR.Sedr(adata.obsm['X_pca'], graph_dict, mode='clustering', device=device)
    using_dec = True
    if using_dec:
        sedr_net.train_with_dec(N=1)
    else:
        sedr_net.train_without_dec(N=1)
    sedr_feat, _, _, _ = sedr_net.process()
    adata.obsm['SEDR'] = sedr_feat
    
    # SEDR.mclust_R(adata, n_clusters, use_rep='SEDR', key_added='SEDR')
    sc.pp.neighbors(adata ,n_neighbors=20)
    eval_res = search_resolution(adata, n_clusters, res_start=0.1, res_end=1, res_step=0.01)
    sc.tl.leiden(adata, key_added='SEDR', resolution=eval_res)
    print('Clustering finished')
    
    # Evaluate clustering
    time_taken = time.time() - start_time
    current, peak = tracemalloc.get_traced_memory()
    memory_used = peak / (1024 ** 2)
    tracemalloc.stop()

    results = evaluate_clustering(adata, df_meta, time_taken, memory_used, dir_out, pred_key='SEDR')
    print(f'ARI = {results["ARI"]:.4f}')

    # Plot clustering
    fig, axes = plt.subplots(1, 1, figsize=(6, 6))
    sc.pl.spatial(adata, color='SEDR', ax=axes, show=False, spot_size=20)
    axes.set_title(f'SEDR (ARI={results["ARI"]:.4f})')
    axes.axis('off')
    plt.tight_layout()
    plt.savefig(os.path.join(dir_out, 'clustering.pdf'), dpi=300, bbox_inches='tight')

    # Plot UMAP
    sc.pp.neighbors(adata, use_rep='SEDR', metric='cosine')
    sc.tl.umap(adata)

    fig, ax = plt.subplots(1, 2, figsize=(8, 3))
    sc.pl.umap(adata, color='layer_guess', ax=ax[0], show=False)
    sc.pl.umap(adata, color='SEDR', ax=ax[1], show=False)
    ax[0].set_title('Manual Annotation')
    ax[1].set_title('SEDR')
    for a in ax:
        a.set_aspect(1)
    plt.tight_layout()
    plt.savefig(os.path.join(dir_out, 'umap.pdf'), format='pdf', dpi=300, bbox_inches='tight')    

    low_dim_data = pd.DataFrame(adata.obsm['SEDR'], index=adata.obs.index)
    low_dim_data.to_csv(f'{dir_out}/low_dim_data.csv')
    adata.obs.to_csv(f'{dir_out}/cell_metadata.csv')
    umap_coords = adata.obsm["X_umap"]
    spot_ids = adata.obs_names
    umap_df = pd.DataFrame(umap_coords, columns=["UMAP1", "UMAP2"])
    umap_df["spot_id"] = spot_ids
    umap_df = umap_df[["spot_id", "UMAP1", "UMAP2"]]
    umap_df.to_csv(f'{dir_out}/spatial_umap_coords.csv')    
    
    print(f"================ Finished Processing {section_id} ======================")