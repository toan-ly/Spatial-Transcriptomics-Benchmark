import os
import pandas as pd
import anndata
import numpy as np

def load_mHypothalamus(root_dir='/home/lytq/Spatial-Transcriptomics-Benchmark/data/mHypothalamus', section_id='0.26'):
    cnts_file = os.path.join(root_dir, 'MERFISH_Animal1_cnts.xlsx')
    info_file = os.path.join(root_dir, 'MERFISH_Animal1_info.xlsx')

    cnts_xls = pd.ExcelFile(cnts_file)
    df_cnts = pd.read_excel(cnts_xls, section_id)

    info_xls = pd.ExcelFile(info_file)
    df_info = pd.read_excel(info_xls, section_id)
    
    obs_ = df_info
    if len(df_info.columns) == 5:
        obs_.columns = ['psuedo_barcodes', 'x', 'y', 'cell_class', 'Neuron_cluster_ID']
    elif len(df_info.columns) == 6:
        obs_.columns = ['psuedo_barcodes', 'x', 'y', 'cell_class', 'Neuron_cluster_ID', 'layer_guess']
    
    obs_.index = obs_['psuedo_barcodes'].tolist()
    var_ = df_cnts.iloc[:, 0]
    var_ = pd.DataFrame(var_)
    # print(var_)
    
    ad = anndata.AnnData(X=df_cnts.iloc[:,1:].T, obs=obs_, var=var_)
    spatial = np.vstack((ad.obs['x'].to_numpy(), ad.obs['y'].to_numpy()))
    ad.obsm['spatial'] = spatial.T
    return ad
    
    