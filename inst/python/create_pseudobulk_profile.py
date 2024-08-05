# check dependencies
REQUIRED_PACKAGES = ['pandas', 'numpy', 'scanpy','anndata','decoupler','scipy']
for package in REQUIRED_PACKAGES:
	try:
		importlib.import_module(package)
		print(f'{package} is installed')
	except ImportError:
		print(f'{package} not installed. Installing now...')
		subprocess.check_call([sys.executable, "-m", "pip", "install", package])
		print(f'{package} installed successfully')
		#importlib.import_module(package)

# import dependencies
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
import sys
import os
import decoupler as dc
from scipy.io import mmread

# get args
sparse_matrix_infile = sys.argv[1]
outdir = sys.argv[2]
meta_infile = sys.argv[3]
cluster_col_str = sys.argv[4]
gene_names_infile = sys.argv[5]
condition_col_str = sys.argv[6]
pid_col_str = sys.argv[7]
test_method = sys.argv[8]

# read temp files
rna_ct = mmread(sparse_matrix_infile)
rna_genes = pd.read_csv(gene_names_infile)['v1']
rna_obs = pd.read_csv(meta_infile)

# remove temp files
try:
    os.remove(sparse_matrix_infile)
except OSError:
    pass
try:
    os.remove(meta_infile)
except OSError:
    pass
try:
    os.remove(gene_names_infile)
except OSError:
    pass

# assemble the anndata object
adata_ct = rna_ct.T.tocsr()
adata = anndata.AnnData(X = adata_ct, obs = rna_obs)
adata.var_names = rna_genes.astype(str)
adata.obs_names = adata.obs.barcode
adata.var_names_make_unique()
adata.obs_names_make_unique()

# form id column
adata.obs['joined_pid'] = adata.obs[pid_col_str]
if test_method == 'cluster_identity':
    adata.obs['joined_pid'] = adata.obs['joined_pid'].astype(str)
if test_method == 'cluster_by_category':
    adata.obs['joined_pid'] = adata.obs['joined_pid'].astype(str) + "_" + adata.obs[condition_col_str].astype(str)
if test_method == 'cluster_by_condition':
    adata.obs['joined_pid'] = adata.obs['joined_pid'].astype(str) + "_" + adata.obs[condition_col_str].astype(str)

# create pseudobulk
adata.obs = adata.obs.astype(str)
adata.layers['counts'] = adata.X.copy()
pdata_sum = dc.get_pseudobulk(
    adata,
    sample_col='joined_pid',
    groups_col=cluster_col_str,
    layer='counts',
    mode='sum',
    min_cells=10,
    min_counts=1000,
    skip_checks=True)

# get elements
pdata_ct_sum = pdata_sum.X
pdata_obs = pdata_sum.obs
pdata_var = pdata_sum.var

# write elements
np.savetxt(os.path.join(outdir,'__pseudobulk_sum_counts__.csv'), pdata_ct_sum, delimiter=",")
pdata_obs.to_csv(os.path.join(outdir,'__pseudobulk_obs__.csv'))
pdata_var.to_csv(os.path.join(outdir,'__pseudobulk_var__.csv'))