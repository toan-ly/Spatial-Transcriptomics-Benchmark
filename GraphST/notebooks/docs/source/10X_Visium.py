import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
from pathlib import Path
import matplotlib.pyplot as plt

from GraphST import GraphST
from GraphST.utils import clustering

import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.append('/home/lytq/Spatial-Transcriptomics-Benchmark/utils')
from sdmbench import compute_ARI, compute_NMI, compute_CHAOS, compute_PAS, compute_ASW, compute_HOM, compute_COM

import time
import psutil
import tracemalloc


def load_data(file_fold: str, is_breast_cancer: bool) -> sc.AnnData:
    """Load spatial data and metadata"""
    adata = sc.read_visium(file_fold, count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.var_names_make_unique()
    
    metadata_file = Path(file_fold) / 'metadata.tsv'
    df_meta = pd.read_csv(metadata_file, sep='\t')
    if not is_breast_cancer:
        adata.obs['layer_guess'] = df_meta['layer_guess']
    else:
        adata.obs['fine_annot_type'] = df_meta['fine_annot_type']
    
    return adata, df_meta

def train_model(adata: sc.AnnData, device: torch.device) -> sc.AnnData:
    """Train GraphST model"""
    model = GraphST.GraphST(adata, device=device)
    adata = model.train()
    return adata

def perform_clustering(adata: sc.AnnData, n_clusters: int, tool: str, radius: int = 50):
    """Perform spatial clustering using specified tool"""
    if tool == 'mclust':
        clustering(adata, n_clusters, radius=radius, method=tool, refinement=True) # For DLPFC dataset, we use optional refinement step.
    elif tool in ['leiden', 'louvain']:
        clustering(adata, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False)
    else:
        raise ValueError(f"Unsupported clustering tool: {tool}")

def save_results(adata: sc.AnnData, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    low_dim_data = pd.DataFrame(adata.obsm['feat'], index=adata.obs.index)
    # expression_data = pd.DataFrame(adata.layers['count'], index=adata.obs.index, columns=adata.var.index)
    cell_metadata = adata.obs

    low_dim_data.to_csv(f"{output_dir}/low_dim_data.csv")
    # expression_data.T.to_csv(f"{output_dir}/expression_matrix.csv")
    cell_metadata.to_csv(f"{output_dir}/cell_metadata.csv")

def evaluate_clustering(adata: sc.AnnData, df_meta, time_taken: float, memory_used: float, output_dir: str, dataset: str) -> dict:
    """Evaluate clustering using sdmbench"""
    gt_key = 'ground_truth'
    pred_key = 'domain'
    if dataset == 'DLPFC':
        adata.obs['ground_truth'] = df_meta['layer_guess'].values
    elif dataset == 'BRCA':
        adata.obs['ground_truth'] = df_meta['fine_annot_type'].values
    adata = adata[~pd.isnull(adata.obs['ground_truth'])]
    
    results = {
        "ARI": compute_ARI(adata, gt_key, pred_key),
        "AMI": compute_NMI(adata, gt_key, pred_key),
        "Homogeneity": compute_HOM(adata, gt_key, pred_key),
        "Completeness": compute_COM(adata, gt_key, pred_key),
        "ASW": compute_ASW(adata, pred_key),
        "CHAOS": compute_CHAOS(adata, pred_key),
        "PAS": compute_PAS(adata, pred_key),
        "Time": time_taken,
        "Memory": memory_used
    }
    
    df_results = pd.DataFrame([results])
    df_results.to_csv(os.path.join(output_dir, "metrics.csv"), index=False)
    return results

def visulize_results(adata: sc.AnnData, ARI: float, output_dir: str, dataset: str = None):
    os.makedirs(output_dir, exist_ok=True)

    # Spatial clustering visualization
    spatial_file = os.path.join(output_dir, "clustering.pdf")
    fig, axes = plt.subplots(1, 1, figsize=(8, 4))
    # fig, axes = plt.subplots(figsize=(6, 6), dpi=300)
    if dataset == 'DLPFC':
        sc.pl.spatial(adata, 
                img_key="hires", 
                color=["ground_truth", "domain"], 
                title=[f"Ground truth ({dataset})", "GraphST (ARI=%.4f)"%ARI],
                show=True)
    elif dataset == 'BRCA':
        # sc.pl.spatial(adata, 
        #         img_key="hires", 
        #         color=["ground_truth", "domain"], 
        #         # palette="tab20",
        #         title=[f"Ground truth", "GraphST (ARI=%.4f)"%ARI],
        #         show=True)
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
    plt.savefig(spatial_file, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # UMAP visualization
    sc.pp.neighbors(adata, use_rep='emb_pca', n_neighbors=10)
    sc.tl.umap(adata) 
    if dataset == 'DLPFC':
        fig, axes = plt.subplots(1, 2, figsize=(8, 4))
        sc.pl.umap(adata, color='layer_guess', ax=axes[0], show=False)
        sc.pl.umap(adata, color='domain', ax=axes[1], show=False)

        axes[0].set_title('Manual Annotation')
        axes[1].set_title('GraphST')
        for ax in axes:
            ax.set_aspect(1)
        plt.tight_layout()
    else:
        fig, axes = plt.subplots(1, 1, figsize=(6, 6))
        sc.pl.umap(adata, color='domain', ax=axes, show=False)
        axes.set_title('GraphST')
        plt.tight_layout()     
        
    umap_file = os.path.join(output_dir, "umap.pdf")
    plt.savefig(umap_file, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
        
    # Save UMAP coordinates
    umap_coords = adata.obsm["X_umap"]
    spot_ids = adata.obs_names
    umap_df = pd.DataFrame(umap_coords, columns=["UMAP1", "UMAP2"])
    umap_df["spot_id"] = spot_ids
    umap_df = umap_df[["spot_id", "UMAP1", "UMAP2"]]
    umap_df.to_csv(os.path.join(output_dir, "spatial_umap_coords.csv"), index=False)

    print(f"Visualizations saved to {output_dir}")


def process_dataset(dataset: str, n_clusters: int, radius: int, tool: str, base_dir: str, output_dir: str, device: torch.device):
    start_time = time.time()
    # start_mem = psutil.Process().memory_info().rss  # Initial memory usage
    tracemalloc.start()
    
    file_fold = os.path.join(base_dir, dataset) if dataset != 'BRCA' else base_dir
    output_dir = os.path.join(output_dir, dataset)
    is_breast_cancer = 'BRCA' in base_dir
    
    # Load data
    adata, df_meta = load_data(file_fold, is_breast_cancer)
    
    # Train model
    adata = train_model(adata, device)
    
    # Perform clustering
    perform_clustering(adata, n_clusters, tool, radius)

    # Save results
    save_results(adata, output_dir)
    
    end_time = time.time()
    end_mem = psutil.Process().memory_info().rss  # Final memory usage
    time_taken = end_time - start_time
    # memory_used = (end_mem - start_mem) / (1024 ** 2)  # Convert to MB    
    size, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    memory_used = peak / (1024 ** 2)  # Convert to MB
    
    # Evaluate clustering
    # if not is_breast_cancer:
    results = evaluate_clustering(adata, df_meta, time_taken, memory_used, output_dir, dataset)
    ARI = results['ARI']
    print(f'Dataset {dataset} ARI: {ARI}')
    visulize_results(adata, ARI, output_dir, dataset)
    # else:
    #     results = evaluate_clustering(adata, df_meta, time_taken, memory_used, output_dir, dataset)
    #     visulize_results(adata, 0, output_dir)
    
def main():
    # Device setup
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

    # base_dir = '/home/lytq/Spatial-Transcriptomics-Benchmark/data/DLPFC'
    # datasets = os.listdir(base_dir)
    # datasets = [d for d in datasets if d.isdigit()]
    # # print(datasets)
    # print(len(datasets))
    
    base_dir = '/home/lytq/Spatial-Transcriptomics-Benchmark/data/BRCA1'
    datasets = ['BRCA']
    
    radius = 50
    tool = 'mclust'
    data_name = base_dir.split('/')[-1]    
    output_dir = "/home/lytq/Spatial-Transcriptomics-Benchmark/RESULTS/" + data_name + "/GraphST"

    for dataset in datasets:
        # n_clusters = 5 if dataset in ['151669', '151670', '151671', '151672'] else 7       
        n_clusters = 20
        print(f'Processing dataset {dataset}...')
        # print(f'Number of clusters: {n_clusters}')
        process_dataset(dataset, n_clusters, radius, tool, base_dir, output_dir, device)
        

if __name__ == '__main__':
    main()