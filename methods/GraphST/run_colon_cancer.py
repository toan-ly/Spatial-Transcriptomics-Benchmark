import warnings
warnings.filterwarnings("ignore")

import json
import os
import random
import time
import tracemalloc
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import torch

from GraphST import GraphST
from GraphST.utils import clustering

import sys
sys.path.append('/home/lytq/Spatial-Transcriptomics-Benchmark/utils')
from evaluate import evaluate_clustering
from load_st_data import load_colon_visium_hd

# SEEDS = [42, 123, 456, 789, 2024]
SEEDS = [2024]
SAMPLE_NAME = 'visium_hd_cancer_colon_square_016um'
N_CLUSTERS = 6
REFINE_RADIUS = 6   # label refinement neighbors (not microns); match STAGATE k≈6


def set_seed(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)


def main():
    data_path = Path('/home/lytq/Spatial-Transcriptomics-Benchmark/data/colon_cancer/')
    output_path = Path('/home/lytq/Spatial-Transcriptomics-Benchmark/results/colon_cancer/')
    dir_input = data_path / SAMPLE_NAME
    spatial_dir = dir_input / 'spatial'
    device = torch.device('cuda:4' if torch.cuda.is_available() else 'cpu')

    print(f'Loading {SAMPLE_NAME}...')
    adata_raw, Ann_df = load_colon_visium_hd(dir_input, spatial_dir)

    for seed in SEEDS:
        print(f"\n==============================\nRUNNING SEED: {seed}\n==============================")
        set_seed(seed)

        dir_out = output_path / f'{seed}/GraphST'
        dir_out.mkdir(parents=True, exist_ok=True)

        time_start = time.time()
        tracemalloc.start()

        adata = adata_raw.copy()
        adata.layers['count'] = adata.X.copy()
        sc.pp.highly_variable_genes(adata, flavor='seurat_v3', layer='count', n_top_genes=3000)

        model = GraphST.GraphST(
            adata,
            device=device,
            random_seed=seed,
            datatype='10X',   
        )
        adata = model.train()

        clustering(
            adata, n_clusters=N_CLUSTERS, radius=REFINE_RADIUS,
            method='mclust', refinement=True,
        )

        time_taken = time.time() - time_start
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        results = evaluate_clustering(
            adata, Ann_df, time_taken, peak / (1024 ** 2), dir_out,
            pred_key='domain', gt_df_key='fine_annot_type',
        )
        ARI = results['ARI']
        print(f'ARI = {ARI:.4f}')

        # Plot spatial clustering
        fig = plt.figure(figsize=(12, 4))
        ax1 = plt.subplot(1, 2, 1)
        sc.pl.spatial(adata, 
                img_key="hires", 
                color="ground_truth",
                title="Ground truth",
                show=False,
                ax=ax1)
        
        ax2 = plt.subplot(1, 2, 2)
        sc.pl.spatial(adata, 
                img_key="hires", 
                color="domain",
                title=f"GraphST (ARI={ARI:.4f})",
                show=False,
                ax=ax2)

        plt.tight_layout()
        plt.savefig(dir_out / 'clustering.pdf', bbox_inches='tight', dpi=300)

        # UMAP visualization
        sc.pp.neighbors(adata, use_rep='emb_pca', n_neighbors=10)
        sc.tl.umap(adata)
        fig, axes = plt.subplots(1, 2, figsize=(10, 3))
        sc.pl.umap(adata, color='Ground Truth', ax=axes[0], show=False)
        sc.pl.umap(adata, color='domain', ax=axes[1], show=False)
        axes[0].set_title('Manual Annotation')
        axes[1].set_title('GraphST')
        plt.tight_layout()
        plt.savefig(dir_out / 'umap.pdf', format='pdf', dpi=300, bbox_inches='tight')

        pd.DataFrame(adata.obsm['emb'], index=adata.obs.index).to_csv(dir_out / 'low_dim_data.csv')
        adata.obs.to_csv(dir_out / 'cell_metadata.csv')
        umap_df = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
        umap_df['spot_id'] = adata.obs_names
        umap_df[['spot_id', 'UMAP1', 'UMAP2']].to_csv(dir_out / 'spatial_umap_coords.csv', index=False)


if __name__ == '__main__':
    main()