import os
import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt
import anndata

import sys
sys.path.append('/home/lytq/Spatial-Transcriptomics-Benchmark/utils')
from evaluate import evaluate_clustering
from load_st_data import load_mHypothalamus

import time
import psutil
import tracemalloc

sys.path.append('/home/lytq/Spatial-Transcriptomics-Benchmark/SpaceFlow/SpaceFlow-main')
from SpaceFlow import SpaceFlow

def search_resolution(adata, n_clusters):
    for res in sorted(list(np.arange(0.1, 1, 0.01)), reverse=True):
        sc.tl.leiden(adata, random_state=0, resolution=res)
        if len(pd.DataFrame(adata.obs['leiden']).leiden.unique()) == n_clusters:
            print(f"Resolution: {res}")
            break
    return res

BASE_PATH = Path('/home/lytq/Spatial-Transcriptomics-Benchmark/data/mHypothalamus')
output_path = Path('/home/lytq/Spatial-Transcriptomics-Benchmark/Results/MERFISH/mHypothalamus/SpaceFlow')

sample_list = ['-0.04', '-0.09', '-0.14', '-0.19', '-0.24']
        

for sample in sample_list:
    print(f"================ Start Processing {sample} ======================")
    dir_input = Path(f'{BASE_PATH}/{sample}/')
    dir_output = Path(f'{output_path}/{sample}/')
    dir_output.mkdir(parents=True, exist_ok=True)
    
    n_clusters = 8
    
    start_time = time.time()
    tracemalloc.start()
    
    adata = load_mHypothalamus(BASE_PATH, sample)
    gt_df = adata.obs

    sc.pp.filter_genes(adata, min_cells=3)
    sf = SpaceFlow.SpaceFlow(adata=adata)
    sf.preprocessing_data(n_top_genes=3000)
    
    sf.train(
        spatial_regularization_strength=0.1,
        z_dim=50,
        lr=1e-3,
        epochs=1000,
        max_patience=50,
        min_stop=100,
        random_seed=42,
        gpu=4,
        regularization_acceleration=True,
        edge_subset_sz=1000000,
        embedding_save_filepath=os.path.join(dir_output, "low_dim_data.csv"),    
    )
    
    sc.pp.neighbors(adata, n_neighbors=50)
    sc.tl.umap(adata)
    
    embedding = anndata.AnnData(sf.embedding)
    sc.pp.neighbors(embedding, n_neighbors=50, use_rep='X')
    res = search_resolution(embedding, n_clusters)
        
    sf.segmentation(domain_label_save_filepath=os.path.join(dir_output, "domain_labels.csv"),                   
                    n_neighbors=50,
                    resolution=res)

    pred = pd.read_csv(os.path.join(dir_output, "domain_labels.csv"), header=None)
    adata.obs['pred'] = pred.values
    adata.obs["pred"] = adata.obs["pred"].astype("category") 
    adata.obsm['SpaceFlow'] = sf.embedding
    
    end_time = time.time()
    time_taken = end_time - start_time
    current, peak = tracemalloc.get_traced_memory()
    memory_used = peak / (1024 ** 2) # MB
    tracemalloc.stop()

    metrics = evaluate_clustering(adata, gt_df, time_taken, memory_used, dir_output, pred_key='pred')
    
    # Plot spatial clusters
    # fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    # sc.pl.spatial(adata, color='layer_guess', ax=axes[0], show=False)
    fig, axes = plt.subplots(1, 1, figsize=(6, 6))
    sc.pl.spatial(adata, color='pred', ax=axes, show=False, spot_size=20)
    # axes[0].set_title('Manual Annotation')
    axes.set_title(f'SpaceFlow (ARI: {metrics["ARI"]:.4f})')
    handles, labels = axes.get_legend_handles_labels()
    new_labels = [str(int(label) + 1) if label.isdigit() else label for label in labels]
    axes.legend(handles, new_labels, loc='center left', frameon=False, bbox_to_anchor=(1, 0.5))
    axes.axis('off')
    plt.tight_layout()
    plt.savefig(dir_output / 'clustering.pdf', format='pdf', dpi=300, bbox_inches='tight')
    
    # Plot UMAP
    fig, ax = plt.subplots(1, 2, figsize=(8, 3))
    sc.pl.umap(adata, color='ground_truth', ax=ax[0], show=False)
    sc.pl.umap(adata, color='pred', ax=ax[1], show=False)
    ax[0].set_title('Manual Annotation')
    ax[1].set_title('SpaceFlow')
    handles, labels = ax[1].get_legend_handles_labels()
    new_labels = [str(int(label) + 1) if label.isdigit() else label for label in labels]
    ax[1].legend(handles, new_labels, loc='center left', frameon=False, bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig(dir_output / 'umap.pdf', format='pdf', dpi=300, bbox_inches='tight')    
    
    low_dim_data = pd.DataFrame(adata.obsm["SpaceFlow"], index=adata.obs.index)
    low_dim_data.to_csv(dir_output / 'low_dim_data.csv', index=False) 
    adata.obs.to_csv(dir_output / 'cell_metadata.csv', index=False)
    umap_coords = adata.obsm["X_umap"]
    spot_ids = adata.obs_names
    umap_df = pd.DataFrame(umap_coords, columns=["UMAP1", "UMAP2"])
    umap_df["spot_id"] = spot_ids
    umap_df = umap_df[["spot_id", "UMAP1", "UMAP2"]]
    umap_df.to_csv(dir_output / "spatial_umap_coords.csv", index=False)
    
    print(f"================ End Processing {sample} ======================")
    