# Spatial Analysis

## create scanpy object
import scanpy as sc
import os
import scimap as sm
import matplotlib.pyplot as plt
from sklearn.cluster import MiniBatchKMeans
from spatial_analysis_functions import *
from collections import Counter
from multiprocessing import Pool

## load expression matrix and clinical information
expPath = "/mnt/data/lyx/IMC/analysis/clustering/annotate_allcells.csv"
clinicalPath = "/mnt/data/lyx/IMC/clinical.csv"
qcPath = "/mnt/data/lyx/IMC/IMC_CRC_QC.csv"

adata = MergeClinical(expPath, clinicalPath, qcPath, 35)
adata

## IM/CT ROIs to scimap spatial analysis
adata = sc.read_h5ad("/mnt/data/lyx/IMC/analysis/spatial/adata.h5ad")
for tissue in ["CT","TAT","IM"]:
    adata_ = adata[adata.obs['Tissue'].isin([tissue]),]
#adata = adata[adata.obs['Tissue'].isin(["IM","CT"]),]

## scanpy pre-processing
#adata.var_names = adata.var["Marker"]


## scimap analysis pipeline
### 1. spatial count: compute a neighbourhood matrix using radius method or KNN method
#n_radiusVec = [10,20,30]
#for radius in n_radiusVec:
#    adata = sm.tl.spatial_count(adata, x_coordinate='x', y_coordinate='y',
#                                phenotype='SubType', method='radius', radius=radius,
#                                imageid='ID', subset=None, label=('spatial_count_radius_'+str(radius)))

n_neighborsVec = [10,20,30]
for n_neighbors in n_neighborsVec:
    adata = sm.tl.spatial_count(adata, x_coordinate='x', y_coordinate='y',
                                phenotype='SubType', method='knn',knn=n_neighbors,
                                imageid='ID', subset=None, label=('spatial_count_knn_'+str(n_neighbors)))

### 2. spatial expression: compute a neighbourhood weighted matrix based on the expression values.
#adata = sm.tl.spatial_expression(adata, x_coordinate='x', y_coordinate='y',
#                                 method='radius', radius=30, imageid='ID',
#                                 use_raw=False, log=False, subset=None, label='spatial_expression_radius')

#adata = sm.tl.spatial_expression(adata, x_coordinate='x', y_coordinate='y',
#                                 method='knn', knn=n_neighbors, imageid='ID',
#                                 use_raw=False, subset=None, label='spatial_expression_knn_'+str(n_neighbors))

### 3. spatial aggregates: find regions of aggregration of similar cells.
# Use the purity parameter to fine-tune percent of similar cells within a given radius
#adata = sm.tl.spatial_aggregate(adata, x_coordinate='x', y_coordinate='y',
#                                phenotype='SubType', method='radius', radius=30, purity=95,
#                                imageid='ID', subset=None, label='spatial_aggregate_95')

#adata = sm.tl.spatial_aggregate(adata, x_coordinate='x', y_coordinate='y',
#                                phenotype='SubType', method='knn', knn=10, purity=60,
#                                imageid='ID', subset=None, label='spatial_aggregate_knn')

### 4. Clustering the spatial count and expression
#for radius in n_radiusVec:
#    kmeans = MiniBatchKMeans(n_clusters=10, random_state=0).fit(adata.uns[('spatial_count_radius_'+str(radius))])

#    cluster_labels = list(map(str,kmeans.labels_))
#    cluster_labels = list(map(lambda orig_string: 'kmeans' + '_' + orig_string, cluster_labels))
#    adata.obs[('kmeans_radius_'+str(radius))] = cluster_labels

for n_neighbors in n_neighborsVec:
    kmeans = MiniBatchKMeans(n_clusters=10, random_state=0).fit(adata.uns[('spatial_count_knn_'+str(n_neighbors))])

    cluster_labels = list(map(str,kmeans.labels_))
    cluster_labels = list(map(lambda orig_string: 'kmeans' + '_' + orig_string, cluster_labels))
    adata.obs[('kmeans_knn_'+str(n_neighbors))] = cluster_labels

#### save the clustering result
savePath = "/mnt/data/lyx/IMC/analysis/spatial/"
if not os.path.exists(savePath):
    os.mkdir(savePath)
adata.obs.to_csv(savePath+"cellular_neighbor_"+tissue+".csv")

### 5. Motif: spatial_lda
adata = sm.tl.spatial_lda(adata, x_coordinate='x', y_coordinate='y',
                          phenotype='SubType', imageid='ID', num_motifs=10, radius=50)

lda_res = adata.uns['spatial_lda']
motif_res = [lda_res.loc[ind].sort_values(
    ascending=False).index[0] for ind in (lda_res.index)]
adata.obs['motif'] = motif_res

sm.pl.stacked_barplot(adata, x_axis='motif', y_axis='SubType',
                      method='percent', figsize=(15, 12))
plt.savefig(os.path.join(savePath,"Celltype_Motif.png"),bbox_inches='tight')

# Official Pipeline
### 1. Spatial distance analysis between cell type
### Calculate distances between cell types
adata = sm.tl.spatial_distance(adata,
                               x_coordinate='x', y_coordinate='y',
                               z_coordinate=None,
                               phenotype='SubType',
                               subset=None,
                               imageid='ID',
                               label='spatial_distance')

plt.rcParams['figure.figsize'] = [3, 1]
sm.pl.spatial_distance(adata, phenotype='SubType')
plt.savefig(os.path.join(savePath,"Celltype_distance.png"),bbox_inches='tight')

# from certain celltype
sm.pl.spatial_distance(adata, method='numeric', imageid='RFS_status',
                       phenotype='SubType', distance_from='CD366+ Tumor cell', log=True)
plt.savefig(os.path.join(savePath,"Celltype_distance_CD366Tumorcell.png"),bbox_inches='tight')

### 2. Spatial co-occurance analysis
# Using the radius method to identify local neighbours compute P-values
adata = sm.tl.spatial_interaction(
    adata, x_coordinate='x', y_coordinate='y',
    phenotype='SubType', method='radius', imageid='ID',  # method can be choose as knn
    radius=30, label='spatial_interaction_radius')

# view results
# spatial_interaction heatmap
sm.pl.spatial_interaction(adata,
                          summarize_plot=True,
                          spatial_interaction='spatial_interaction_radius',
                          row_cluster=True, linewidths=0.75, linecolor='black')

plt.rcParams['figure.figsize'] = [3, 1]
plt.savefig(os.path.join(savePath,"Celltype_spatial_interaction.png"),bbox_inches='tight')

### 3. Quantifying the proximity score
# Calculate the score for proximity between two cell types
adata = sm.tl.spatial_pscore(adata, proximity=['MC4','SC1'],
                             score_by='ID',
                             x_coordinate='x', y_coordinate='y',
                             imageid='ID',
                             phenotype='SubType',
                             method='radius',
                             radius=50,
                             subset=None,
                             label='spatial_pscore')

# Plot only `Proximity Volume` scores
plt.figure(figsize=(10, 5))
sm.pl.spatial_pscore(adata, color='Black', plot_score='Proximity Volume')
plt.savefig(os.path.join(savePath,"Celltype_spatial_ProximityVolume.png"),bbox_inches='tight')

# Plot only `Proximity Density` scores
plt.figure(figsize=(10, 5))
sm.pl.spatial_pscore(adata, color='Black', plot_score='Proximity Density')
plt.savefig(os.path.join(savePath,"Celltype_spatial_ProximityDensity.png"),bbox_inches='tight')

# voronoi plot
plt.rcParams['figure.figsize'] = [15, 10]
sm.pl.voronoi(adata, color_by='spatial_pscore',
              x_coordinate='x', y_coordinate='y',
              voronoi_edge_color='black',
              voronoi_line_width=0.3,
              voronoi_alpha=0.8,
              size_max=5000,
              overlay_points=None,
              plot_legend=True,
              legend_size=6)

plt.savefig("/mnt/data/lyx/IMC/spatial/Celltype_spatial_Proximity_Voronoi.png")

adata.write(os.path.join(savePath, "adata_spatial.h5ad"))

adata = sc.read_h5ad(os.path.join(savePath, "adata_IM_spatial.h5ad"))

## Permutation test
IDs = [x for x in Counter(adata.obs["ID"]).keys()]

multiprocess = True
if multiprocess:
    multi_threads = 70
    with Pool(multi_threads) as p:
        p.starmap(PermutationTestMain, [(adata, ID, 22,os.path.join("/mnt/data/lyx/IMC/analysis/spatial/permutation",ID)) for ID in IDs])

## Merge permutation results
ResultPath = "/mnt/data/lyx/IMC/analysis/spatial/permutation/"
types = [str(x) for x in Counter(adata.obs["SubType"]).keys()]

## all 
adata = sc.read_h5ad("/mnt/data/lyx/IMC/analysis/spatial/adata.h5ad")
for tissue in ["CT","TAT","IM"]:
    adata_ = adata[adata.obs['Tissue'].isin([tissue]),]

    n_neighborsVec = [10,20,30]
    for n_neighbors in n_neighborsVec:
        adata_ = sm.tl.spatial_count(adata_, x_coordinate='x', y_coordinate='y',
                                    phenotype='SubType', method='knn',knn=n_neighbors,
                                    imageid='ID', subset=None, label=('spatial_count_knn_'+str(n_neighbors)))

    for n_neighbors in n_neighborsVec:
        kmeans = MiniBatchKMeans(n_clusters=10, random_state=0).fit(adata_.uns[('spatial_count_knn_'+str(n_neighbors))])

        cluster_labels = list(map(str,kmeans.labels_))
        cluster_labels = list(map(lambda orig_string: 'CNP' + '_' + orig_string, cluster_labels))
        adata_.obs[('kmeans_knn_'+str(n_neighbors))] = cluster_labels

    #### save the clustering result
    savePath = "/mnt/data/lyx/IMC/analysis/spatial/"
    if not os.path.exists(savePath):
        os.mkdir(savePath)
    adata_.obs.to_csv(savePath+"cellular_neighbor_"+tissue+".csv")

    ## Permutation test
    IDs = [x for x in Counter(adata_.obs["ID"]).keys()]

    multiprocess = True
    if multiprocess:
        multi_threads = 99
        with Pool(multi_threads) as p:
            p.starmap(PermutationTestMain, [(adata_, ID, 22,os.path.join(("/mnt/data/lyx/IMC/analysis/spatial/permutation_"+tissue),ID)) for ID in IDs])
