## nohup.py
import scanpy as sc
import os
import scimap as sm
from sklearn.cluster import MiniBatchKMeans
from spatial_analysis_functions import *
from multiprocessing import Pool

## load expression matrix and clinical information
expPath = "/mnt/data/lyx/IMC/analysis/clustering/annotate_allcells.csv"
clinicalPath = "/mnt/data/lyx/IMC/clinical.csv"
qcPath = "/mnt/data/lyx/IMC/IMC_CRC_QC.csv"

adata = MergeClinical(expPath, clinicalPath, qcPath, 35)

adata = sc.read_h5ad("/mnt/data/lyx/IMC/analysis/spatial/adata.h5ad")

for tissue in ["All","IM","CT","TAT"]:
    if tissue == "All":
        adata_ = adata
    else:
        adata_ = adata[adata.obs['Tissue'].isin([tissue]),]

    n_neighborsVec = [5,10,15,20,25,30]
    for n_neighbors in n_neighborsVec:
        adata_ = sm.tl.spatial_count(adata_, x_coordinate='x', y_coordinate='y',
                                    phenotype='SubType', method='knn',knn=n_neighbors,
                                    imageid='ID', subset=None, label=('spatial_count_knn_'+str(n_neighbors)))



    for n_neighbors in n_neighborsVec:
        kmeans = MiniBatchKMeans(n_clusters=10, random_state=0).fit(adata_.uns[('spatial_count_knn_'+str(n_neighbors))])

        cluster_labels = list(map(str,kmeans.labels_))
        cluster_labels = list(map(lambda orig_string: 'CNP' + '_' + orig_string, cluster_labels))
        adata_.obs[('CNP'+str(n_neighbors))] = cluster_labels

    #### save the clustering result
    savePath = "/mnt/data/lyx/IMC/analysis/spatial/"
    if not os.path.exists(savePath):
        os.mkdir(savePath)
    adata_.obs.to_csv(savePath+"cellular_exp_neighbor_"+tissue+".csv")

    ## Permutation test
    IDs = [x for x in Counter(adata_.obs["ID"]).keys()]

    multiprocess = True
    if multiprocess:
       multi_threads = 24
       with Pool(multi_threads) as p:
           p.starmap(PermutationTestMain, [(adata_, ID, 22,os.path.join(("/mnt/data/lyx/IMC/analysis/spatial/permutation_"+tissue),ID)) for ID in IDs])