## nohup.py
import scanpy as sc
import os
from spatial_analysis_functions import *
from collections import Counter
from multiprocessing import Pool

## load expression matrix and clinical information
expPath = "/mnt/data/lyx/IMC/analysis/clustering/annotate_allcells.csv"
clinicalPath = "/mnt/data/lyx/IMC/clinical.csv"
qcPath = "/mnt/data/lyx/IMC/IMC_CRC_QC.csv"

adata = MergeClinical(expPath, clinicalPath, qcPath, 35)
adata

savePath = "/mnt/data/lyx/IMC/analysis/spatial/"
if not os.path.exists(savePath):
    os.mkdir(savePath)

adata.write(os.path.join(savePath, "adata.h5ad"))

## IM ROIs spatial analysis
adata = sc.read_h5ad(os.path.join(savePath,"adata.h5ad"))
adata = adata[adata.obs['Tissue'] == "IM",]

## Permutation test
IDs = [x for x in Counter(adata.obs["ID"]).keys()]

multiprocess = True
if multiprocess:
    multi_threads = 65
    with Pool(multi_threads) as p:
        p.starmap(PermutationTestMain, [(adata, ID, 22, os.path.join("/mnt/data/lyx/IMC/analysis/spatial/permutation",ID)) for ID in IDs])