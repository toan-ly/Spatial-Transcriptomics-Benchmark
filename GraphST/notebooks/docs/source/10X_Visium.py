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


def load_data(file_fold: str, is_breast_cancer: bool) -> sc.AnnData:
    """Load spatial data and metadata"""
    adata = sc.read_visium(file_fold, count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.var_names_make_unique()
    
    metadata_file = Path(file_fold) / 'metadata.tsv'
    df_meta = pd.read_csv(metadata_file, sep='\t')
    if not is_breast_cancer:
        adata.obs['layer_guess'] = df_meta['layer_guess']
    
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
    expression_data = pd.DataFrame(adata.layers['count'], index=adata.obs.index, columns=adata.var.index)
    cell_metadata = adata.obs

    low_dim_data.to_csv(f"{output_dir}/low_dim_data.csv")
    expression_data.T.to_csv(f"{output_dir}/expression_matrix.csv")
    cell_metadata.to_csv(f"{output_dir}/cell_metadata.csv")

def evaluate_clustering(adata: sc.AnnData, df_meta) -> float:
    """Evaluate clustering using ARI"""
    df_meta_layer = df_meta['layer_guess']
    adata.obs['ground_truth'] = df_meta_layer.values
    adata = adata[~pd.isnull(adata.obs['ground_truth'])]

    ARI = metrics.adjusted_rand_score(adata.obs['domain'], adata.obs['ground_truth'])
    adata.uns['ARI'] = ARI
    return ARI 

def visulize_results(adata: sc.AnnData, ARI: float, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    # Spatial clustering visualization
    spatial_file = os.path.join(output_dir, "spatial_clustering.png")
    fig, axes = plt.subplots(1, 1, figsize=(8, 4))
    if ARI:
        sc.pl.spatial(adata, 
                img_key="hires", 
                color=["ground_truth", "domain"], 
                title=["Ground truth", "ARI=%.4f"%ARI],
                show=True)
    else:
        sc.pl.spatial(adata, 
                img_key="hires", 
                color="domain", 
                palette="tab20",
                title="Clustering",
                show=True)
    plt.tight_layout()
    plt.savefig(spatial_file)
    plt.close()
    
    # UMAP visualization
    sc.pp.neighbors(adata, use_rep='emb_pca', n_neighbors=10)
    sc.tl.umap(adata) 
    if ARI:
        fig, axes = plt.subplots(1, 2, figsize=(8, 3))
        sc.pl.umap(adata, color='layer_guess', ax=axes[0], show=False)
        sc.pl.umap(adata, color='domain', ax=axes[1], show=False)

        axes[0].set_title('Manual Annotation')
        axes[1].set_title('Clustering')
        for ax in axes:
            ax.set_aspect(1)
        plt.tight_layout()
    else:
        fig, axes = plt.subplots(1, 1, figsize=(4, 3))
        sc.pl.umap(adata, color='domain', ax=axes, show=False)
        axes.set_title('Clustering')
        plt.tight_layout()     
    umap_file = os.path.join(output_dir, "umap.png")
    plt.savefig(umap_file)
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
    file_fold = os.path.join(base_dir, dataset)
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
    
    # Evaluate clustering
    if not is_breast_cancer:
        ARI = evaluate_clustering(adata, df_meta)
        print(f'Dataset {dataset} ARI: {ARI}')
        visulize_results(adata, ARI, output_dir)
    else:
        visulize_results(adata, 0, output_dir)
    
def main():
    # Device setup
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

    # base_dir = '/home/lytq/GraphST/data/DLPFC'
    # datasets = os.listdir(base_dir)
    # datasets = [d for d in datasets if d.isdigit()]
    # # print(datasets)
    # print(len(datasets))
    
    base_dir = '/home/lytq/GraphST/data/BRCA1'
    datasets = ['V1_Human_Breast_Cancer_Block_A_Section_1']
    
    radius = 50
    tool = 'mclust'
    data_name = base_dir.split('/')[-1]    
    output_dir = "/home/lytq/GraphST/results/" + data_name

    for dataset in datasets:
        # n_clusters = 5 if dataset in ['151669', '151670', '151671', '151672'] else 7       
        n_clusters = 20
        print(f'Processing dataset {dataset}...')
        # print(f'Number of clusters: {n_clusters}')
        process_dataset(dataset, n_clusters, radius, tool, base_dir, output_dir, device)
        

if __name__ == '__main__':
    main()