import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os

import torch

import STAGATE_pyG as STAGATE

import time
import tracemalloc
import sys
sys.path.append('/home/lytq/Spatial-Transcriptomics-Benchmark/utils')
from evaluate import evaluate_clustering
from load_st_data import load_mHypothalamus
    
def main():
    # the location of R (used for the mclust clustering)
    os.environ['R_HOME'] = '/home/lytq/.conda/envs/stagate/lib/R'
    os.environ['R_USER'] = '/.conda/envs/stagate/lib/python3.10/site-packages/rpy2'

    data_path = '/home/lytq/Spatial-Transcriptomics-Benchmark/data/mHypothalamus'
    output_path = '/home/lytq/Spatial-Transcriptomics-Benchmark/Results/MERFISH/mHypothalamus/STAGATE'

    data_names = ['-0.04', '-0.09', '-0.14', '-0.19', '-0.24']    
    device = torch.device('cuda:5' if torch.cuda.is_available() else 'cpu')
    for section_id in data_names:
        print(f'Processing {section_id}...')
        n_clusters = 8

        dir_out = f'{output_path}/{section_id}'
        os.makedirs(dir_out, exist_ok=True)

        time_start = time.time()
        tracemalloc.start()

        # Load data
        adata = load_mHypothalamus(data_path, section_id)
        df_meta = adata.obs

        #Normalization
        sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        ## Constructing the spatial network
        STAGATE.Cal_Spatial_Net(adata, rad_cutoff=150)
        STAGATE.Stats_Spatial_Net(adata)

        ## Runing STAGATE
        adata = STAGATE.train_STAGATE(adata, device=device)
        print(f'Finish training STAGATE for {section_id}')
        
        ## Perform clustering
        sc.pp.neighbors(adata, use_rep='STAGATE')
        sc.tl.umap(adata)
        adata = STAGATE.mclust_R(adata, used_obsm='STAGATE', num_cluster=n_clusters)

        ## Calculate clustering metrics
        obs_df = adata.obs.dropna()
        
        time_taken = time.time() - time_start
        current, peak = tracemalloc.get_traced_memory()
        memory_used = peak / (1024 ** 2)
        tracemalloc.stop()
        
        clustering_results = evaluate_clustering(adata, df_meta, time_taken, memory_used, dir_out, pred_key='mclust')
        ARI = clustering_results["ARI"]
        print('Adjusted rand index = %.2f' %ARI)

        ## Save UMAP
        # plt.rcParams["figure.figsize"] = (4, 3)
        # sc.pl.umap(adata, color=["Ground Truth", "mclust"], title=["Ground Truth", 'STAGATE (ARI=%.2f)'%ARI], frameon=False)

        fig, axes = plt.subplots(1, 2, figsize=(8, 3))
        sc.pl.umap(adata, color='layer_guess', ax=axes[0], show=False)
        sc.pl.umap(adata, color='mclust', ax=axes[1], show=False)
        axes[0].set_title('Manual Annotation')
        axes[1].set_title(f'STAGATE')
        for ax in axes:
            ax.set_aspect(1)
        plt.tight_layout()
        plt.savefig(os.path.join(dir_out, 'umap.pdf'), dpi=300, bbox_inches='tight')
        
        # Plot spatial clustering
        plt.rcParams["figure.figsize"] = (6, 6)
        sc.pl.spatial(adata, 
                      color=["mclust"], 
                      title=['STAGATE (ARI=%.4f)'%ARI],
                      spot_size=20,)
        plt.axis('off')
        plt.savefig(os.path.join(dir_out, 'clustering.pdf'), bbox_inches='tight', dpi=300)

        ## Save results
        adata.obs['STAGATE'] = adata.obs['mclust']
        # adata.write(f'{dir_out}/result.h5ad')
        # adata.obs.to_csv(f'{dir_out}/metadata.tsv', sep='\t')

        # df = pd.DataFrame(data=adata.obsm['STAGATE'], index=adata.obs.index)
        # df.to_csv(f'{dir_out}/PCs.tsv', sep='\t')
        
        ## Save the trajectory
        # used_adata = adata[adata.obs['Ground Truth']!='nan',]
        # used_adata = used_adata[~used_adata.obs['Ground Truth'].isna()]
        # sc.tl.paga(used_adata, groups='Ground Truth')
        # plt.rcParams["figure.figsize"] = (4, 3)
        # sc.pl.paga_compare(used_adata, legend_fontsize=10, frameon=False, size=20,
        #                 title=section_id+'_STAGATE', legend_fontoutline=2, show=False)
        # plt.savefig(os.path.join(dir_out, f'{section_id}_trajectory.png'), bbox_inches='tight', dpi=300)        

        low_dim_data = pd.DataFrame(adata.obsm['STAGATE'], index=adata.obs.index)
        # expression_data = pd.DataFrame(adata.layers['count'], index=adata.obs.index, columns=adata.var.index)
        cell_metadata = adata.obs

        low_dim_data.to_csv(f"{dir_out}/low_dim_data.csv")
        # expression_data.T.to_csv(f"{dir_out}/expression_matrix.csv")
        cell_metadata.to_csv(f"{dir_out}/cell_metadata.csv")
        
        umap_coords = adata.obsm["X_umap"]
        spot_ids = adata.obs_names
        umap_df = pd.DataFrame(umap_coords, columns=["UMAP1", "UMAP2"])
        umap_df["spot_id"] = spot_ids
        umap_df = umap_df[["spot_id", "UMAP1", "UMAP2"]]
        umap_df.to_csv(os.path.join(dir_out, "spatial_umap_coords.csv"))



if __name__ == '__main__':
    main()