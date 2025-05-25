import os,csv,re,sys
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
import random, torch
from sklearn import metrics
import cv2
import matplotlib.pyplot as plt
from pathlib import Path
import sys
sys.path.append('/home/lytq/Spatial-Transcriptomics-Benchmark/utils')
from evaluate import evaluate_clustering
from load_st_data import load_mHypothalamus

import time
import psutil
import tracemalloc


BASE_PATH = Path('/home/lytq/Spatial-Transcriptomics-Benchmark/data/mHypothalamus')
output_path = Path('/home/lytq/Spatial-Transcriptomics-Benchmark/Results/MERFISH/mHypothalamus/SpaGCN')

sample_list = ['-0.04', '-0.09', '-0.14', '-0.19', '-0.24']
ARI_list = []
for sample_name in sample_list:
    print(f"================ Start Processing {sample_name} ======================")

    dir_input = Path(f'{BASE_PATH}/{sample_name}/')
    dir_output = Path(f'{output_path}/{sample_name}/')
    dir_output.mkdir(parents=True, exist_ok=True)

    n_clusters = 8
        
    time_start = time.time()
    tracemalloc.start()
    
    ##### read data
    adata = load_mHypothalamus(BASE_PATH, sample_name)

    #Set coordinates
    x_array=adata.obs["x"].tolist()
    y_array=adata.obs["y"].tolist()
    x_pixel=x_array
    y_pixel=y_array

    #Calculate adjacent matrix
    b=49
    a=1
    adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, image=None, beta=b, alpha=a, histology=False)
    # np.savetxt(f'{dir_output}/adj.csv', adj, delimiter=',')


    ##### Spatial domain detection using SpaGCN
    spg.prefilter_genes(adata, min_cells=3) # avoiding all genes are zeros
    spg.prefilter_specialgenes(adata)
    #Normalize and take log for UMI
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)

    ### 4.2 Set hyper-parameters
    p=0.5 
    spg.test_l(adj,[1, 10, 100, 500, 1000])
    l=spg.find_l(p=p,adj=adj,start=0.1, end=1000,sep=0.01, tol=0.01)
    n_clusters=n_clusters
    r_seed=t_seed=n_seed=100
    res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, 
                        t_seed=t_seed, n_seed=n_seed)

    ### 4.3 Run SpaGCN
    clf=spg.SpaGCN()
    clf.set_l(l)
    #Set seed
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    #Run
    clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
    y_pred, prob=clf.predict()
    adata.obs["pred"]= y_pred
    adata.obs["pred"]=adata.obs["pred"].astype('category')
    #Do cluster refinement(optional)
    adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
    refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
    adata.obs["refined_pred"]=refined_pred
    adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
    
    df_meta = adata.obs

    sc.pp.neighbors(adata, n_neighbors=10)
    sc.tl.umap(adata)
    
    time_taken = time.time() - time_start
    size, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    memory_used = peak / (1024 ** 2) # MB
    
    clustering_results = evaluate_clustering(adata, df_meta, time_taken, memory_used, dir_output, pred_key='refined_pred')
    
    fig, axes = plt.subplots(1, 1, figsize=(6, 6))
    sc.pl.spatial(adata, color='refined_pred', ax=axes, spot_size=20, show=False)
    axes.set_title(f'SpaGCN (ARI = {clustering_results["ARI"]:.4f})')
    axes.axis('off')
    plt.tight_layout()
    plt.savefig(f'{dir_output}/clustering.pdf', dpi=300, bbox_inches='tight')
    
    
    fig, axes = plt.subplots(1,2,figsize=(4*2, 3))
    sc.pl.umap(adata, color='layer_guess', ax=axes[0], show=False)
    sc.pl.umap(adata, color='refined_pred', ax=axes[1], show=False)
    axes[0].set_title('Manual Annotation')
    axes[1].set_title('SpaGCN')
    for ax in axes:
        ax.set_aspect(1)
    plt.tight_layout()
    plt.savefig(f'{dir_output}/umap.pdf', dpi=300, bbox_inches='tight')
    
    low_dim_data = pd.DataFrame(adata.obsm['X_pca'], index=adata.obs.index)
    cell_metadata = adata.obs
    low_dim_data.to_csv(f"{dir_output}/low_dim_data.csv")
    cell_metadata.to_csv(f"{dir_output}/cell_metadata.csv")
    
    umap_coords = adata.obsm["X_umap"]
    spot_ids = adata.obs_names
    umap_df = pd.DataFrame(umap_coords, columns=["UMAP1", "UMAP2"])
    umap_df["spot_id"] = spot_ids
    umap_df = umap_df[["spot_id", "UMAP1", "UMAP2"]]
    umap_df.to_csv(os.path.join(dir_output, "spatial_umap_coords.csv"))
    
    # df_meta = df_meta[~pd.isnull(df_meta['layer_guess'])]
    ARI = clustering_results['ARI']
    print('===== Project: {} ARI score: {:.3f}'.format(sample_name, ARI))
