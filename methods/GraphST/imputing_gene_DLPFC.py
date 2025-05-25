import os
import torch
import pandas as pd
import scanpy as sc
import numpy as np
from sklearn import metrics
import multiprocessing as mp
import warnings
warnings.filterwarnings("ignore")

from scipy.stats import pearsonr
from skimage.metrics import structural_similarity as ssim
from sklearn.metrics import mean_squared_error as mse
from scipy.sparse import csr_matrix

from libpysal.weights import KNN
from esda.moran import Moran

import glob 

from GraphST import GraphST

def compute_morans_i(exp_matrix, w):
    # morans_i_list = []
    # for i in range(exp_matrix.shape[1]): # Iterate over genes
    #     moran = Moran(exp_matrix[:, i], w)
    #     morans_i_list.append(moran.I)
    # return np.mean(morans_i_list) # Return the average Moran's I value

    return np.mean([Moran(exp_matrix[:, i], w).I for i in range(exp_matrix.shape[1])])


device = torch.device("cuda:5" if torch.cuda.is_available() else "cpu")

data_folder = '/home/lytq/Spatial-Transcriptomics-Benchmark/data/DLPFC/'
output_path = '/home/lytq/Spatial-Transcriptomics-Benchmark/Results/Imputation/GraphST/metrics.csv'


files = glob.glob(data_folder + '/*')
idx = 1
for file in files:
    print('=='*20)
    print(f"[{idx}] Processing file: {file}")
    idx += 1
    
    # Load data
    section_id = file.split('/')[-1]
    adata = sc.read_visium(file, count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.var_names_make_unique()
    
    adata_before = adata
    
    # Model training
    model = GraphST.GraphST(adata, device=device)
    adata = model.train()
    
    adata_after = adata.obsm['emb']
    orig_genes = adata_before.var_names.copy()
    highly_variable_genes = adata.var_names[adata.var['highly_variable']]
    adata_before = adata_before[:, highly_variable_genes].copy()


    gene_df2 = pd.DataFrame(adata_after,
                            index=adata.obs.index,
                            columns=highly_variable_genes)

    gene_df1 = pd.DataFrame(adata_before.X.toarray(),
                            index=adata.obs.index,
                            columns=adata_before.var.index)
                        
    
    # Compute Pearson Correlation Coefficient
    assert gene_df1.shape == gene_df2.shape
    vec1 = gene_df1.values.flatten()
    vec2 = gene_df2.values.flatten()

    pcc, _ = pearsonr(vec1, vec2)
    print(f"Pearson Correlation Coefficient: {pcc:.4f}")

    # Compute Structural Similarity Index (SSIM) 
    arr1 = gene_df1.values
    arr2 = gene_df2.values

    ssim_value = ssim(arr1, arr2, data_range=arr2.max() - arr2.min())
    print(f"Structural Similarity Index: {ssim_value:.4f}")
    
    # Compute Root Mean Squared Error (RMSE)
    rmse_value = np.sqrt(mse(arr1, arr2))
    print(f"Root Mean Squared Error: {rmse_value:.4f}")
    
    # Compute Moran's I
    # Create a KNN spatial weights object
    spatial_coords = adata.obsm['spatial']
    w = KNN.from_array(spatial_coords, k=5) 

    morans_i_before = compute_morans_i(arr1, w)
    morans_i_after = compute_morans_i(arr2, w)
    print(f"Moran's I before: {morans_i_before:.4f}")
    print(f"Moran's I after: {morans_i_after:.4f}")
    
    # Save the results
    results = {
        'Model': 'GraphST',
        'Dataset': section_id,
        "Moran's I before": morans_i_before,
        "Moran's I after": morans_i_after,
        'PCC': pcc,
        'SSIM': ssim_value,
        'RMSE': rmse_value
    }
    
    # Write results to CSV
    results_df = pd.DataFrame([results])
    results_df.to_csv(output_path, mode='a', header=False, index=False)