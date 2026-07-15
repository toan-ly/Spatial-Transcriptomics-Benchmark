import warnings
warnings.filterwarnings('ignore')

import os
import random
import time
import tracemalloc
from pathlib import Path

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import torch

from SpaceFlow.SpaceFlow import SpaceFlow

import sys
sys.path.append('/home/lytq/Spatial-Transcriptomics-Benchmark/utils')
from evaluate import evaluate_clustering
from load_st_data import load_colon_visium_hd

SEEDS = [42, 123, 456, 789, 2024]
SAMPLE_NAME = 'visium_hd_cancer_colon_square_016um'
N_CLUSTERS = 6
GPU = 4


def set_seed(seed):
    np.random.seed(seed)
    random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def search_resolution(adata, n_clusters):
    for res in sorted(list(np.arange(0.01, 2.0, 0.02)), reverse=True):
        sc.tl.leiden(adata, random_state=0, resolution=res)
        if len(adata.obs['leiden'].unique()) == n_clusters:
            print(f'Resolution: {res}')
            return res
    print('Warning: exact cluster count not found; using last resolution')
    return res


def main():
    data_path = Path('/home/lytq/Spatial-Transcriptomics-Benchmark/data/colon_cancer/')
    output_path = Path('/home/lytq/Spatial-Transcriptomics-Benchmark/results/colon_cancer/')
    dir_input = data_path / SAMPLE_NAME
    spatial_dir = dir_input / 'spatial'

    print(f'Loading {SAMPLE_NAME}...')
    adata_raw, Ann_df = load_colon_visium_hd(dir_input, spatial_dir)

    for seed in SEEDS:
        print(f"\n==============================\nRUNNING SEED: {seed}\n==============================")
        set_seed(seed)

        dir_out = output_path / f'{seed}/SpaceFlow'
        dir_out.mkdir(parents=True, exist_ok=True)

        time_start = time.time()
        tracemalloc.start()

        adata = adata_raw.copy()
        sc.pp.filter_genes(adata, min_cells=3)

        sf = SpaceFlow.SpaceFlow(adata=adata)
        sf.preprocessing_data(n_top_genes=3000, n_neighbors=10)

        sf.train(
            spatial_regularization_strength=0.1,
            z_dim=50,
            lr=1e-3,
            epochs=1000,
            max_patience=50,
            min_stop=100,
            random_seed=seed,
            gpu=GPU,
            regularization_acceleration=True,
            edge_subset_sz=1000000,
            embedding_save_filepath=str(dir_out / 'low_dim_data.tsv'),
        )

        embedding = anndata.AnnData(sf.embedding)
        sc.pp.neighbors(embedding, n_neighbors=50, use_rep='X')
        res = search_resolution(embedding, N_CLUSTERS)

        sf.segmentation(
            domain_label_save_filepath=str(dir_out / 'domain_labels.csv'),
            n_neighbors=50,
            resolution=res,
        )

        pred = pd.read_csv(dir_out / 'domain_labels.csv', header=None)
        adata.obs['pred'] = pred.iloc[:, 0].values
        adata.obs['pred'] = adata.obs['pred'].astype('category')
        adata.obsm['SpaceFlow'] = sf.embedding

        sc.pp.neighbors(adata, n_neighbors=50, use_rep='SpaceFlow')
        sc.tl.umap(adata)

        time_taken = time.time() - time_start
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        memory_used = peak / (1024 ** 2)

        results = evaluate_clustering(
            adata, Ann_df, time_taken, memory_used, str(dir_out),
            pred_key='pred', gt_df_key='fine_annot_type',
        )
        ARI = results['ARI']
        print(f'ARI = {ARI:.4f}')

        adata.obs['pred_shift'] = (adata.obs['pred'].astype(int) + 1).astype(str)

        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        sc.pl.spatial(
            adata, color='pred_shift', library_id=SAMPLE_NAME,
            ax=ax, show=False, spot_size=2,
        )
        ax.set_title(f'SpaceFlow (ARI={ARI:.4f})')
        ax.axis('off')
        plt.tight_layout()
        plt.savefig(dir_out / 'clustering.pdf', dpi=300, bbox_inches='tight')

        fig, axes = plt.subplots(1, 2, figsize=(10, 3))
        sc.pl.umap(adata, color='Ground Truth', ax=axes[0], show=False)
        sc.pl.umap(adata, color='pred_shift', ax=axes[1], show=False)
        axes[0].set_title('Manual Annotation')
        axes[1].set_title('SpaceFlow')
        plt.tight_layout()
        plt.savefig(dir_out / 'umap.pdf', dpi=300, bbox_inches='tight')

        pd.DataFrame(adata.obsm['SpaceFlow'], index=adata.obs.index).to_csv(dir_out / 'low_dim_data.csv')
        adata.obs.to_csv(dir_out / 'cell_metadata.csv')
        umap_df = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
        umap_df['spot_id'] = adata.obs_names
        umap_df[['spot_id', 'UMAP1', 'UMAP2']].to_csv(dir_out / 'spatial_umap_coords.csv', index=False)

        del sf, adata, embedding
        torch.cuda.empty_cache()

        print(f'Finished seed {seed}')


if __name__ == '__main__':
    main()