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

SEEDS = [42, 123, 456, 789, 2024]
SAMPLE_NAME = 'visium_hd_cancer_colon_square_016um'
N_CLUSTERS = 6
MAX_SPOTS = None
REFINE_RADIUS = 6   # label refinement neighbors (not microns); match STAGATE k≈6


def set_seed(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)


def _resolve_data_dir(dir_input: Path) -> Path:
    dir_input = Path(dir_input)
    qc_dir = dir_input / 'qc'
    return qc_dir if (qc_dir / 'counts.mtx').exists() else dir_input


def load_colon_visium_hd(dir_input, spatial_dir, max_spots=None, seed=42):
    dir_input = Path(dir_input)
    spatial_dir = Path(spatial_dir)
    data_dir = _resolve_data_dir(dir_input)
    obs = pd.read_csv(data_dir / 'observations.tsv', sep='\t', index_col=0)
    var = pd.read_csv(data_dir / 'features.tsv', sep='\t', index_col=0)
    coords = pd.read_csv(data_dir / 'coordinates.tsv', sep='\t', index_col=0)
    labels = pd.read_csv(dir_input / 'labels.tsv', sep='\t', index_col=0)
    print('Loading count matrix (may take several minutes)...')
    X = scipy.io.mmread(data_dir / 'counts.mtx').tocsr()
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs_names = obs.index.astype(str)
    adata.var_names = var.index.astype(str)
    adata.var_names_make_unique()
    selected = adata.obs['selected'].astype(str).str.lower() == 'true'
    adata = adata[selected].copy()
    adata.obs['Ground Truth'] = labels.loc[adata.obs_names, 'label'].values
    adata = adata[adata.obs['Ground Truth'] != 'Outside'].copy()
    adata.obsm['spatial'] = coords.loc[adata.obs_names, ['x', 'y']].to_numpy()
    sf = json.loads((spatial_dir / 'scalefactors_json.json').read_text())
    library_id = dir_input.name
    adata.uns['spatial'] = {
        library_id: {
            'images': {
                'hires': plt.imread(spatial_dir / 'tissue_hires_image.png'),
                'lowres': plt.imread(spatial_dir / 'tissue_lowres_image.png'),
            },
            'scalefactors': sf,
        }
    }
    if max_spots is not None and adata.n_obs > max_spots:
        from sklearn.model_selection import train_test_split
        idx = np.arange(adata.n_obs)
        chosen, _ = train_test_split(
            idx, train_size=max_spots,
            stratify=adata.obs['Ground Truth'].values, random_state=seed,
        )
        adata = adata[chosen].copy()
    Ann_df = pd.DataFrame({'fine_annot_type': adata.obs['Ground Truth'].values})
    print(f'Loaded {adata.n_obs} spots x {adata.n_vars} genes')
    return adata, Ann_df


def main():
    data_path = Path('/home/lytq/Spatial-Transcriptomics-Benchmark/data/colon_cancer/visium_hd_cancer_colon')
    output_path = Path('/home/lytq/Spatial-Transcriptomics-Benchmark/results')
    dir_input = data_path / SAMPLE_NAME
    spatial_dir = dir_input / 'spatial'
    device = torch.device('cuda:4' if torch.cuda.is_available() else 'cpu')

    print(f'Loading {SAMPLE_NAME}...')
    adata_raw, Ann_df = load_colon_visium_hd(dir_input, spatial_dir, max_spots=MAX_SPOTS)

    for seed in SEEDS:
        print(f"\n==============================\nRUNNING SEED: {seed}\n==============================")
        set_seed(seed)

        dir_out = output_path / 'colon_cancer' / str(seed) / 'GraphST'
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

        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        sc.pl.spatial(adata, color='domain', library_id=SAMPLE_NAME,
                      ax=ax, show=False, spot_size=2)
        ax.set_title(f'GraphST (ARI={ARI:.4f})')
        ax.axis('off')
        plt.tight_layout()
        plt.savefig(dir_out / 'umap.pdf' if False else dir_out / 'clustering.pdf', dpi=300)

        sc.pp.neighbors(adata, use_rep='emb_pca', n_neighbors=10)
        sc.tl.umap(adata)
        fig, axes = plt.subplots(1, 2, figsize=(10, 3))
        sc.pl.umap(adata, color='Ground Truth', ax=axes[0], show=False)
        sc.pl.umap(adata, color='domain', ax=axes[1], show=False)
        axes[0].set_title('Manual Annotation')
        axes[1].set_title('GraphST')
        plt.tight_layout()
        plt.savefig(dir_out / 'umap.pdf', dpi=300, bbox_inches='tight')

        pd.DataFrame(adata.obsm['emb'], index=adata.obs.index).to_csv(dir_out / 'low_dim_data.csv')
        adata.obs.to_csv(dir_out / 'cell_metadata.csv')
        umap_df = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
        umap_df['spot_id'] = adata.obs_names
        umap_df[['spot_id', 'UMAP1', 'UMAP2']].to_csv(dir_out / 'spatial_umap_coords.csv', index=False)


if __name__ == '__main__':
    main()