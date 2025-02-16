from sdmbench import compute_ARI, compute_NMI, compute_CHAOS, compute_PAS, compute_ASW, compute_HOM, compute_COM
import scanpy as sc
import pandas as pd
import os

def evaluate_clustering(
    adata: sc.AnnData, 
    df_meta, 
    time_taken: float, 
    memory_used: float, 
    output_dir: str = None,
    gt_key: str ='ground_truth', 
    pred_key: str ='pred', 
    gt_df_key: str ='layer_guess',
) -> dict:
    
    """Evaluate clustering using sdmbench"""
    adata.obs[gt_key] = df_meta[gt_df_key].values
    adata = adata[~pd.isnull(adata.obs[gt_key])]
    
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
    if output_dir:
        df_results.to_csv(os.path.join(output_dir, "metrics.csv"), index=False)
    return df_results