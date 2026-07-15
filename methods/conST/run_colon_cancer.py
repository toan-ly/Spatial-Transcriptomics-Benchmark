import warnings
warnings.filterwarnings('ignore')

import argparse
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

from conST.src.graph_func import graph_construction
from conST.src.utils_func import adata_preprocess, mk_dir, res_search_fixed_clus
from conST.src.training import conST_training

import sys
sys.path.append('/home/lytq/Spatial-Transcriptomics-Benchmark/utils')
from evaluate import evaluate_clustering
from load_st_data import load_colon_visium_hd

SEEDS = [42, 123, 456, 789, 2024]
SAMPLE_NAME = 'visium_hd_cancer_colon_square_016um'
N_CLUSTERS = 6


def set_seed(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True


def refine(sample_id, pred, dis, shape='hexagon'):
    refined_pred = []
    pred = pd.DataFrame({'pred': pred}, index=sample_id)
    dis_df = pd.DataFrame(dis, index=sample_id, columns=sample_id)
    if shape == 'hexagon':
        num_nbs = 6
    elif shape == 'square':
        num_nbs = 4
    else:
        print("Shape not recongized, shape='hexagon' for Visium data, 'square' for ST data.")
    for i in range(len(sample_id)):
        index = sample_id[i]
        dis_tmp = dis_df.loc[index, :].sort_values(ascending=False)
        nbs = dis_tmp[0:num_nbs + 1]
        nbs_pred = pred.loc[nbs.index, 'pred']
        self_pred = pred.loc[index, 'pred']
        v_c = nbs_pred.value_counts()
        if (v_c.loc[self_pred] < num_nbs / 2) and (np.max(v_c) > num_nbs / 2):
            refined_pred.append(v_c.idxmax())
        else:
            refined_pred.append(self_pred)
    return refined_pred


def build_params(device):
    parser = argparse.ArgumentParser()
    parser.add_argument('--k', type=int, default=10)
    parser.add_argument('--knn_distanceType', type=str, default='euclidean')
    parser.add_argument('--epochs', type=int, default=200)
    parser.add_argument('--cell_feat_dim', type=int, default=300)
    parser.add_argument('--feat_hidden1', type=int, default=100)
    parser.add_argument('--feat_hidden2', type=int, default=20)
    parser.add_argument('--gcn_hidden1', type=int, default=32)
    parser.add_argument('--gcn_hidden2', type=int, default=8)
    parser.add_argument('--p_drop', type=float, default=0.2)
    parser.add_argument('--use_img', type=bool, default=False)
    parser.add_argument('--img_w', type=float, default=0.1)
    parser.add_argument('--use_pretrained', type=bool, default=True)
    parser.add_argument('--using_mask', type=bool, default=False)
    parser.add_argument('--feat_w', type=float, default=10)
    parser.add_argument('--gcn_w', type=float, default=0.1)
    parser.add_argument('--dec_kl_w', type=float, default=10)
    parser.add_argument('--gcn_lr', type=float, default=0.01)
    parser.add_argument('--gcn_decay', type=float, default=0.01)
    parser.add_argument('--dec_cluster_n', type=int, default=10)
    parser.add_argument('--dec_interval', type=int, default=20)
    parser.add_argument('--dec_tol', type=float, default=0.0)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--beta', type=float, default=100)
    parser.add_argument('--cont_l2l', type=float, default=0.3)
    parser.add_argument('--cont_l2c', type=float, default=0.1)
    parser.add_argument('--cont_l2g', type=float, default=0.1)
    parser.add_argument('--edge_drop_p1', type=float, default=0.1)
    parser.add_argument('--edge_drop_p2', type=float, default=0.1)
    parser.add_argument('--node_drop_p1', type=float, default=0.2)
    parser.add_argument('--node_drop_p2', type=float, default=0.3)
    parser.add_argument('--eval_resolution', type=int, default=1)
    parser.add_argument('--eval_graph_n', type=int, default=20)

    params = parser.parse_args(args=['--k', '20', '--knn_distanceType', 'euclidean', '--epochs', '200'])
    params.device = device
    return params


def main():
    data_path = Path('/home/lytq/Spatial-Transcriptomics-Benchmark/data/colon_cancer/')
    output_path = Path('/home/lytq/Spatial-Transcriptomics-Benchmark/results/colon_cancer/')
    dir_input = data_path / SAMPLE_NAME
    spatial_dir = dir_input / 'spatial'
    device = 'cuda:4' if torch.cuda.is_available() else 'cpu'
    print(f'Using device: {device}')

    print(f'Loading {SAMPLE_NAME}...')
    adata_raw, Ann_df = load_colon_visium_hd(dir_input, spatial_dir)

    params = build_params(device)

    for seed in SEEDS:
        print(f"\n==============================\nRUNNING SEED: {seed}\n==============================")
        set_seed(seed)

        dir_out = output_path / f'{seed}/conST'
        params.save_path = mk_dir(str(dir_out))

        time_start = time.time()
        tracemalloc.start()

        adata_h5 = adata_raw.copy()
        adata_X = adata_preprocess(adata_h5, min_cells=5, pca_n_comps=params.cell_feat_dim)
        graph_dict = graph_construction(adata_h5.obsm['spatial'], adata_h5.shape[0], params)
        params.cell_num = adata_h5.shape[0]

        conST_net = conST_training(adata_X, graph_dict, params, N_CLUSTERS)
        conST_net.pretraining()
        conST_net.major_training()
        conST_embedding = conST_net.get_embedding()

        adata_conST = anndata.AnnData(conST_embedding, obs=adata_h5.obs)
        adata_conST.uns['spatial'] = adata_h5.uns['spatial']
        adata_conST.obsm['spatial'] = adata_h5.obsm['spatial']

        sc.pp.neighbors(adata_conST, n_neighbors=params.eval_graph_n)
        eval_resolution = res_search_fixed_clus(adata_conST, N_CLUSTERS)
        print(f'Resolution: {eval_resolution}')
        cluster_key = 'conST_leiden'
        sc.tl.leiden(adata_conST, key_added=cluster_key, resolution=eval_resolution)

        index = [str(x) for x in range(adata_X.shape[0])]
        dis_matrix = graph_dict['adj_norm'].to_dense().numpy() + np.eye(graph_dict['adj_norm'].shape[0])
        adata_conST.obs['refine'] = refine(
            sample_id=index,
            pred=adata_conST.obs[cluster_key].tolist(),
            dis=dis_matrix,
            shape='square',
        )

        time_taken = time.time() - time_start
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        memory_used = peak / (1024 ** 2)

        results = evaluate_clustering(
            adata_conST, Ann_df, time_taken, memory_used, str(dir_out),
            pred_key='refine', gt_df_key='fine_annot_type',
        )
        ARI = results['ARI']
        print(f'ARI = {ARI:.4f}')

        adata_conST.obs['refine_shift'] = (adata_conST.obs['refine'].astype(int) + 1).astype(str)

        ################# Plotting #################
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        sc.pl.spatial(
            adata_conST, color='refine_shift', library_id=SAMPLE_NAME,
            ax=ax, show=False, spot_size=2,
        )
        ax.set_title(f'conST (ARI={ARI:.4f})')
        ax.axis('off')
        plt.tight_layout()
        plt.savefig(dir_out / 'clustering.pdf', format='pdf', dpi=300, bbox_inches='tight')

        fig, ax = plt.subplots(1, 2, figsize=(8, 3))
        sc.pl.umap(adata_conST, color='Ground Truth', ax=ax[0], show=False)
        sc.pl.umap(adata_conST, color='refine_shift', ax=ax[1], show=False)
        ax[0].set_title('Manual Annotation')
        ax[1].set_title('conST')
        for a in ax:
            a.set_aspect(1)
        plt.tight_layout()
        plt.savefig(dir_out / 'umap.pdf', format='pdf', dpi=300, bbox_inches='tight')

        sc.tl.umap(adata_conST)
        fig, axes = plt.subplots(1, 2, figsize=(10, 3))
        sc.pl.umap(adata_conST, color='Ground Truth', ax=axes[0], show=False)
        sc.pl.umap(adata_conST, color='refine_shift', ax=axes[1], show=False)
        axes[0].set_title('Manual Annotation')
        axes[1].set_title('conST')
        plt.tight_layout()
        plt.savefig(dir_out / 'umap.pdf', dpi=300, bbox_inches='tight')

        adata_conST.obsm['conST'] = conST_embedding
        pd.DataFrame(adata_conST.obsm['conST'], index=adata_conST.obs.index).to_csv(dir_out / 'low_dim_data.csv')
        adata_conST.obs.to_csv(dir_out / 'cell_metadata.csv')
        umap_df = pd.DataFrame(adata_conST.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
        umap_df['spot_id'] = adata_conST.obs_names
        umap_df[['spot_id', 'UMAP1', 'UMAP2']].to_csv(dir_out / 'spatial_umap_coords.csv')

        print(f'Finished running {SAMPLE_NAME} with seed {seed}')

if __name__ == '__main__':
    main()
