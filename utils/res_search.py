import numpy as np
import pandas as pd
import scanpy as sc

def search_resolution(adata, n_clusters, res_start=0.1, res_end=0.5, res_step=0.01):
    for res in sorted(list(np.arange(res_start, res_end, res_step)), reverse=True):
        sc.tl.leiden(adata, random_state=0, resolution=res)
        if len(pd.DataFrame(adata.obs['leiden']).leiden.unique()) == n_clusters:
            print(f"Resolution: {res}")
            break
    return res