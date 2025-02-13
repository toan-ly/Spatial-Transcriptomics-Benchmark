import os
import torch
import pandas as pd
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt
import psutil
import warnings

from GraphST import GraphST
from GraphST.utils import clustering
from sdmbench import compute_ARI, compute_NMI, compute_CHAOS, compute_PAS, compute_ASW, compute_HOM, compute_COM

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
    clustering(adata, n_clusters, radius=radius, method=tool, refinement=True)

def evaluate_clustering(adata: sc.AnnData, df_meta) -> dict:
    """Evaluate clustering using sdmbench"""
    gt_key = 'layer_guess'
    pred_key = 'domain'
    adata.obs['ground_truth'] = df_meta[gt_key].values
    adata = adata[~pd.isnull(adata.obs['ground_truth'])]
    
    results = {
        "ARI": compute_ARI(adata, gt_key, pred_key),
        "AMI": compute_NMI(adata, gt_key, pred_key),
        "Homogeneity": compute_HOM(adata, gt_key, pred_key),
        "Completeness": compute_COM(adata, gt_key, pred_key),
        "ASW": compute_ASW(adata, pred_key),
        "CHAOS": compute_CHAOS(adata, pred_key),
        "PAS": compute_PAS(adata, pred_key)
    }
    return results

def process_dataset(dataset: str, n_clusters: int, radius: int, tool: str, base_dir: str, output_dir: str, device: torch.device):
    file_fold = os.path.join(base_dir, dataset)
    output_dir = os.path.join(output_dir, dataset)
    is_breast_cancer = 'BRCA' in base_dir
    
    adata, df_meta = load_data(file_fold, is_breast_cancer)
    adata = train_model(adata, device)
    perform_clustering(adata, n_clusters, tool, radius)
    
    if not is_breast_cancer:
        results = evaluate_clustering(adata, df_meta)
        print(f'Dataset {dataset} Clustering Metrics: {results}')
    
    print(f'Processing completed for dataset {dataset}')

def main():
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    base_dir = '/home/lytq/GraphST/data/BRCA1'
    datasets = ['V1_Human_Breast_Cancer_Block_A_Section_1']
    output_dir = "/home/lytq/GraphST/results/BRCA1"
    
    for dataset in datasets:
        n_clusters = 20
        process_dataset(dataset, n_clusters, radius=50, tool='mclust', base_dir=base_dir, output_dir=output_dir, device=device)

if __name__ == '__main__':
    main()