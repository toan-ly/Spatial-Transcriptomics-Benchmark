import os
import pandas as pd
import anndata as ad
import numpy as np
import scipy
import json
from pathlib import Path
import matplotlib.pyplot as plt


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
    
    adata = ad.AnnData(X=df_cnts.iloc[:,1:].T, obs=obs_, var=var_)
    spatial = np.vstack((adata.obs['x'].to_numpy(), adata.obs['y'].to_numpy()))
    adata.obsm['spatial'] = spatial.T
    return adata
    
    
def _resolve_data_dir(dir_input: Path) -> Path:
    dir_input = Path(dir_input)
    qc_dir = dir_input / 'qc'
    return qc_dir if (qc_dir / 'counts.mtx').exists() else dir_input

def load_colon_visium_hd(dir_input, spatial_dir):
    dir_input = Path(dir_input)
    spatial_dir = Path(spatial_dir)
    data_dir = _resolve_data_dir(dir_input)

    obs = pd.read_csv(data_dir / 'observations.tsv', sep='\t', index_col=0)
    var = pd.read_csv(data_dir / 'features.tsv', sep='\t', index_col=0)
    coords = pd.read_csv(data_dir / 'coordinates.tsv', sep='\t', index_col=0)
    labels = pd.read_csv(dir_input / 'labels.tsv', sep='\t', index_col=0)

    print('Loading count matrix ...')
    X = scipy.io.mmread(data_dir / 'counts.mtx').tocsr()
    assert X.shape[0] == len(obs) == len(coords), (
        f'shape mismatch: mtx {X.shape[0]}, obs {len(obs)}, coords {len(coords)}'
    )
  
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs_names = obs.index.astype(str)
    adata.var_names = var.index.astype(str)
    adata.var_names_make_unique()
    
    selected = adata.obs['selected'].astype(str).str.lower() == 'true'
    adata = adata[selected].copy()
    adata.obs['Ground Truth'] = labels.loc[adata.obs_names, 'label'].values
    adata = adata[adata.obs['Ground Truth'] != 'Outside'].copy()

    adata.obsm['spatial'] = coords.loc[adata.obs_names, ['0', '1']].to_numpy()
    sf = json.loads((spatial_dir / 'scalefactors_json.json').read_text())
    library_id = dir_input.name
    adata.uns['spatial'] = {
        library_id: {
            'images': {
                'hires': plt.imread(spatial_dir / 'tissue_hires_image.png'),
                'lowres': plt.imread(spatial_dir / 'tissue_lowres_image.png'),
            },
            'scalefactors': sf,
        }
    }
   
    Ann_df = pd.DataFrame(
        {'fine_annot_type': adata.obs['Ground Truth'].values},
        index=adata.obs_names,
    )
    print(f'Loaded {adata.n_obs} spots x {adata.n_vars} genes')
    print(adata.obs['Ground Truth'].value_counts())
    return adata, Ann_df