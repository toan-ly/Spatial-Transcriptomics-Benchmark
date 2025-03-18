import torch
import argparse
import random
import numpy as np
import pandas as pd
import sys

sys.path.append('/home/lytq/Spatial-Transcriptomics-Benchmark/conST/conST-main')
from src.graph_func import graph_construction
from src.utils_func import mk_dir, adata_preprocess, load_ST_file, res_search_fixed_clus, plot_clustering
from src.training import conST_training

import anndata
from sklearn import metrics
import matplotlib.pyplot as plt
import scanpy as sc
import os
import warnings
warnings.filterwarnings('ignore')

sys.path.append('/home/lytq/Spatial-Transcriptomics-Benchmark/utils')
from evaluate import evaluate_clustering
from load_st_data import load_mHypothalamus

import time
import tracemalloc


parser = argparse.ArgumentParser()
parser.add_argument('--k', type=int, default=10, help='parameter k in spatial graph')
parser.add_argument('--knn_distanceType', type=str, default='euclidean',
                    help='graph distance type: euclidean/cosine/correlation')
parser.add_argument('--epochs', type=int, default=200, help='Number of epochs to train.')
parser.add_argument('--cell_feat_dim', type=int, default=300, help='Dim of PCA')
parser.add_argument('--feat_hidden1', type=int, default=100, help='Dim of DNN hidden 1-layer.')
parser.add_argument('--feat_hidden2', type=int, default=20, help='Dim of DNN hidden 2-layer.')
parser.add_argument('--gcn_hidden1', type=int, default=32, help='Dim of GCN hidden 1-layer.')
parser.add_argument('--gcn_hidden2', type=int, default=8, help='Dim of GCN hidden 2-layer.')
parser.add_argument('--p_drop', type=float, default=0.2, help='Dropout rate.')
parser.add_argument('--use_img', type=bool, default=False, help='Use histology images.')
parser.add_argument('--img_w', type=float, default=0.1, help='Weight of image features.')
parser.add_argument('--use_pretrained', type=bool, default=True, help='Use pretrained weights.')
parser.add_argument('--using_mask', type=bool, default=False, help='Using mask for multi-dataset.')
parser.add_argument('--feat_w', type=float, default=10, help='Weight of DNN loss.')
parser.add_argument('--gcn_w', type=float, default=0.1, help='Weight of GCN loss.')
parser.add_argument('--dec_kl_w', type=float, default=10, help='Weight of DEC loss.')
parser.add_argument('--gcn_lr', type=float, default=0.01, help='Initial GNN learning rate.')
parser.add_argument('--gcn_decay', type=float, default=0.01, help='Initial decay rate.')
parser.add_argument('--dec_cluster_n', type=int, default=10, help='DEC cluster number.')
parser.add_argument('--dec_interval', type=int, default=20, help='DEC interval nnumber.')
parser.add_argument('--dec_tol', type=float, default=0.00, help='DEC tol.')

parser.add_argument('--seed', type=int, default=0, help='random seed')
parser.add_argument('--beta', type=float, default=100, help='beta value for l2c')
parser.add_argument('--cont_l2l', type=float, default=0.3, help='Weight of local contrastive learning loss.')
parser.add_argument('--cont_l2c', type=float, default= 0.1, help='Weight of context contrastive learning loss.')
parser.add_argument('--cont_l2g', type=float, default= 0.1, help='Weight of global contrastive learning loss.')

parser.add_argument('--edge_drop_p1', type=float, default=0.1, help='drop rate of adjacent matrix of the first view')
parser.add_argument('--edge_drop_p2', type=float, default=0.1, help='drop rate of adjacent matrix of the second view')
parser.add_argument('--node_drop_p1', type=float, default=0.2, help='drop rate of node features of the first view')
parser.add_argument('--node_drop_p2', type=float, default=0.3, help='drop rate of node features of the second view')

# ______________ Eval clustering Setting ______________
parser.add_argument('--eval_resolution', type=int, default=1, help='Eval cluster number.')
parser.add_argument('--eval_graph_n', type=int, default=20, help='Eval graph kN tol.') 

params =  parser.parse_args(args=['--k', '20', '--knn_distanceType', 'euclidean', '--epochs', '200'])

np.random.seed(params.seed)
torch.manual_seed(params.seed)
torch.cuda.manual_seed(params.seed)
device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
print('Using device: ' + device)
params.device = device

def refine(sample_id, pred, dis, shape="hexagon"):
    refined_pred=[]
    pred=pd.DataFrame({"pred": pred}, index=sample_id)
    dis_df=pd.DataFrame(dis, index=sample_id, columns=sample_id)
    if shape=="hexagon":
        num_nbs=6 
    elif shape=="square":
        num_nbs=4
    else:
        print("Shape not recongized, shape='hexagon' for Visium data, 'square' for ST data.")
    for i in range(len(sample_id)):
        index=sample_id[i]
        dis_tmp=dis_df.loc[index, :].sort_values(ascending=False)
        nbs=dis_tmp[0:num_nbs+1]
        nbs_pred=pred.loc[nbs.index, "pred"]
        self_pred=pred.loc[index, "pred"]
        v_c=nbs_pred.value_counts()
        if (v_c.loc[self_pred]<num_nbs/2) and (np.max(v_c)>num_nbs/2):
            refined_pred.append(v_c.idxmax())
        else:           
            refined_pred.append(self_pred)
    return refined_pred

# set seed before every run
def seed_torch(seed):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
seed_torch(params.seed)

sample_list = ['-0.04', '-0.09', '-0.14', '-0.19', '-0.24']
# sample_list = ['-0.19']

data_root = '/home/lytq/Spatial-Transcriptomics-Benchmark/data/mHypothalamus'
save_root = '/home/lytq/Spatial-Transcriptomics-Benchmark/Results/MERFISH/mHypothalamus/conST'

for data_name in sample_list:

    print(f"================ Start Processing {data_name} ======================")
    path = f'{data_root}/{data_name}'    
    print(path)
    save_path = os.path.join(save_root, data_name)
    params.save_path = mk_dir(save_path)    
    
    start_time = time.time()
    tracemalloc.start()
    

    adata_h5 = load_mHypothalamus(data_root, data_name) 
    
    if params.cell_feat_dim > len(adata_h5.var.index):
        params.cell_feat_dim = len(adata_h5.var.index) - 1
    adata_X = adata_preprocess(adata_h5, min_cells=5, pca_n_comps=params.cell_feat_dim)
    graph_dict = graph_construction(adata_h5.obsm['spatial'], adata_h5.shape[0], params)

    params.cell_num = adata_h5.shape[0]

    n_clusters = 8
    
    conST_net = conST_training(adata_X, graph_dict, params, n_clusters)
    conST_net.pretraining()
    conST_net.major_training()

    conST_embedding = conST_net.get_embedding()

    # clustering
    adata_conST = anndata.AnnData(conST_embedding, obs=adata_h5.obs)
    # adata_conST.uns['spatial'] = adata_h5.uns['spatial']
    adata_conST.obs['layer_guess'] = adata_h5.obs['layer_guess']
    adata_conST.obsm['spatial'] = adata_h5.obsm['spatial']
    df_meta = adata_conST.obs

    sc.pp.neighbors(adata_conST, n_neighbors=params.eval_graph_n)

    print('Finding resolution...')
    eval_resolution = res_search_fixed_clus(adata_conST, n_clusters, increment=0.01, start=0.05, end=1.5)
    print('Resolution:', eval_resolution)
    cluster_key = "conST_leiden"
    sc.tl.leiden(adata_conST, key_added=cluster_key, resolution=eval_resolution)

    print(adata_conST)
    # plot_clustering(adata_conST, cluster_key)
    
    
    # Refine clusters
    # index = np.arange(start=0, stop=adata_X.shape[0]).tolist()
    # index = [str(x) for x in index]
    
    # dis_matrix = graph_dict['adj_norm'].to_dense().numpy() + np.eye(graph_dict['adj_norm'].shape[0])
    # refine = refine(sample_id = index, pred = adata_conST.obs['leiden'].tolist(), dis=dis_matrix)
    # adata_conST.obs['refine'] = refine
    
    # cluster_key = 'refine'
    # # plot_clustering(adata_conST, cluster_key)
    # df_meta['conST_refine'] = adata_conST.obs['refine'].tolist()
    # # df_meta.to_csv(f'{params.save_path}/metadata.tsv', sep='\t', index=False)
    # ARI = metrics.adjusted_rand_score(df_meta['layer_guess'], df_meta['conST_refine'])
    # print('===== Project: {} refined ARI score: {:.3f}'.format(data_name, ARI))
    
    current, peak = tracemalloc.get_traced_memory()
    memory_used = peak / (1024 ** 2) # MB
    time_taken = time.time() - start_time
    tracemalloc.stop()
    
    results = evaluate_clustering(adata_conST, df_meta, time_taken, memory_used, save_path, pred_key=cluster_key)
    
    ################# Plotting #################
    sc.tl.umap(adata_conST)
    # adata_conST.uns['refine_colors'] = adata_conST.uns['layer_guess_colors']
    # adata_conST.obs['refine_shift'] = (adata_conST.obs[cluster_key].astype(int) + 1).astype('category')    
    
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    sc.pl.spatial(adata_conST, color=cluster_key, ax=ax, show=False, spot_size=20)
    ax.set_title(f'conST (ARI = {results["ARI"]:.4f})')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(f'{params.save_path}/clustering.pdf', format='pdf', dpi=300, bbox_inches='tight')
    
    fig, ax = plt.subplots(1, 2, figsize=(8, 3))
    sc.pl.umap(adata_conST, color='layer_guess', ax=ax[0], show=False)
    sc.pl.umap(adata_conST, color=cluster_key, ax=ax[1], show=False)
    ax[0].set_title('Manual Annotation')
    ax[1].set_title('conST')
    for a in ax:
        a.set_aspect(1)
    plt.tight_layout()
    plt.savefig(f'{params.save_path}/umap.pdf', format='pdf', dpi=300, bbox_inches='tight')
    
    adata_conST.obsm['conST'] = conST_embedding
    low_dim_data = pd.DataFrame(adata_conST.obsm['conST'], index=adata_conST.obs.index)
    low_dim_data.to_csv(f'{params.save_path}/low_dim_data.csv')
    adata_conST.obs.to_csv(f'{params.save_path}/cell_metadata.csv')
    umap_coords = adata_conST.obsm["X_umap"]
    spot_ids = adata_conST.obs_names
    umap_df = pd.DataFrame(umap_coords, columns=["UMAP1", "UMAP2"])
    umap_df["spot_id"] = spot_ids
    umap_df = umap_df[["spot_id", "UMAP1", "UMAP2"]]
    umap_df.to_csv(f'{params.save_path}/spatial_umap_coords.csv')
    
    
    print(f"================ Finish Processing {data_name} ======================")
    print(f"Time taken: {time_taken}")
    print(f"Memory used: {memory_used}")
    print("=====================================================================")
