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
from load_st_data import load_mHypothalamus
from evaluate import evaluate_clustering

import time
import psutil
import tracemalloc


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

def visulize_results(adata: sc.AnnData, ARI: float, output_dir: str, dataset: str = None):
    os.makedirs(output_dir, exist_ok=True)

    # Spatial clustering visualization
    spatial_file = os.path.join(output_dir, "clustering.pdf")
    fig, axes = plt.subplots(1, 1, figsize=(6, 6))
    sc.pl.spatial(adata, 
            color='domain',
            title=f"GraphST (ARI=%.4f)"%ARI,
            show=False,
            spot_size=20,)
    plt.tight_layout()
    plt.savefig(spatial_file, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # UMAP visualization
    sc.pp.neighbors(adata, use_rep='emb_pca', n_neighbors=10)
    sc.tl.umap(adata) 

    fig, axes = plt.subplots(1, 2, figsize=(8, 3))
    sc.pl.umap(adata, color='layer_guess', ax=axes[0], show=False)
    sc.pl.umap(adata, color='domain', ax=axes[1], show=False)

    axes[0].set_title('Manual Annotation')
    axes[1].set_title('GraphST')
    for ax in axes:
        ax.set_aspect(1)
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
    
    output_dir = os.path.join(output_dir, dataset)
    
    # Load data
    adata = load_mHypothalamus(base_dir, dataset)
    df_meta = adata.obs
    
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
    results = evaluate_clustering(adata, df_meta, time_taken, memory_used, output_dir, pred_key='domain')
    ARI = results['ARI']
    print(f'Dataset {dataset} ARI: {ARI}')
    visulize_results(adata, ARI, output_dir, dataset)
    # else:
    #     results = evaluate_clustering(adata, df_meta, time_taken, memory_used, output_dir, dataset)
    #     visulize_results(adata, 0, output_dir)
    
def main():
    # Device setup
    device = torch.device('cuda:5' if torch.cuda.is_available() else 'cpu')

    base_dir = '/home/lytq/Spatial-Transcriptomics-Benchmark/data/mHypothalamus'
    datasets = ['-0.04', '-0.09', '-0.14', '-0.19', '-0.24']          

    radius = 50
    tool = 'mclust'
    output_dir = "/home/lytq/Spatial-Transcriptomics-Benchmark/Results/MERFISH/mHypothalamus/GraphST/"

    for dataset in datasets:
        n_clusters = 8
        print(f'Processing dataset {dataset}...')
        process_dataset(dataset, n_clusters, radius, tool, base_dir, output_dir, device)
        

if __name__ == '__main__':
    main()