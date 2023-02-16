# Unsupervised clustering and gating annotation
## scanpy - leiden pipeline

import scanpy as sc
import pandas as pd
import numpy as np
import scimap as sm
import matplotlib.pyplot as plt
import os
from random import shuffle

from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import confusion_matrix
import itertools
import joblib

## load scanpy anndata, to obtain this object, run 3.clustering.r and 5. spatial analysis
adata = sc.read_h5ad("/mnt/data/lyx/IMC/spatial/adata.h5ad")

## random select cells
randomIndex = adata.obs.index.tolist()
shuffle(randomIndex)
randomIndex = randomIndex[0:30000]
adata_ = adata[randomIndex]

del randomIndex

## scanpy pipeline
sc.tl.pca(adata_, svd_solver='arpack')
sc.pp.neighbors(adata_, n_neighbors=30, n_pcs=10) # Computing the neighborhood graph
sc.tl.umap(adata_) # Build a UMAP to visualize the neighbourhood graph
sc.tl.leiden(adata_)

plt.rcParams['figure.figsize'] = [8, 6]
sc.pl.umap(adata_, color=['leiden'], use_raw=False, s=30,show=False) # Plot the UMAP
plt.savefig("/mnt/data/lyx/IMC/clustering/leiden_clustering.png",bbox_inches='tight')

sc.tl.rank_genes_groups(adata_, 'leiden', method='t-test')
plt.rcParams['figure.figsize'] = [8, 6]
sc.pl.rank_genes_groups(adata_, n_genes=10, sharey=False, fontsize=16)
plt.savefig("/mnt/data/lyx/IMC/clustering/leiden_markers.png",bbox_inches='tight')

plt.rcParams['figure.figsize'] = [16, 12]
sc.pl.umap(adata_, color=['CD45', 'CD3', 'CD163','CD68', 'CollagenI', 'FAP', 'SMA', 'VEGF','EpCAM'], use_raw=False, s=30,show=False) # Plot the UMAP
plt.savefig("/mnt/data/lyx/IMC/clustering/leiden_markers.png",bbox_inches='tight')

## Major annotation
new_cluster_names = ['Myeloid cell', 'Tumor cell','Myeloid cell', 'Lymphocyte',
    'Stromal cell', 'Tumor cell','Unknown cell', 'Tumor cell',
    'Tumor cell','Myeloid cell','Myeloid cell','Lymphocyte']
adata_.obs['MajorType']=[new_cluster_names[int(x)] for x in adata_.obs["leiden"].tolist()]
adata_.obs['MajorType']=adata_.obs['MajorType'].astype('category')

plt.rcParams['figure.figsize'] = [8, 6]
sc.pl.umap(adata_, color=['MajorType'], use_raw=False, s=30,show=False) # Plot the UMAP
plt.savefig("/mnt/data/lyx/IMC/clustering/leiden_clustering_majortype.png",bbox_inches='tight')

## Minor annotation
adata_.obs['SubType']='0'

### Myeloid
adata_tem=adata_[adata_.obs["MajorType"]=='Myeloid cell']
sc.pp.neighbors(adata_tem, n_neighbors=30, n_pcs=10) # Computing the neighborhood graph
sc.tl.umap(adata_tem) # Build a UMAP to visualize the neighbourhood graph
sc.tl.leiden(adata_tem)
sc.tl.rank_genes_groups(adata_tem, 'leiden', method='t-test')
plt.rcParams['figure.figsize'] = [8, 6]
sc.pl.rank_genes_groups(adata_tem, n_genes=10, sharey=False, fontsize=16)
plt.savefig("/mnt/data/lyx/IMC/clustering/leiden_Myeloid_markers.png",bbox_inches='tight')

new_cluster_names = ['UMC1', 'UMC2','HLADR+ Myeloid cell', 'cDC',
    'CD11b+ Macrophage', 'CD163+ Macrophage','CD14_16+ Monocyte', 'UMC3','B3GAT1+ Myeloid cell','UMC4']
adata_tem.obs['MinorType']=[new_cluster_names[int(x)] for x in adata_tem.obs["leiden"].tolist()]

adata_.obs.loc[adata_tem.obs.index,"SubType"]=adata_tem.obs['MinorType']

### Tumor cell
adata_tem=adata_[adata_.obs["MajorType"]=='Tumor cell']
sc.pp.neighbors(adata_tem, n_neighbors=30, n_pcs=10) # Computing the neighborhood graph
sc.tl.umap(adata_tem) # Build a UMAP to visualize the neighbourhood graph
sc.tl.leiden(adata_tem)
sc.tl.rank_genes_groups(adata_tem, 'leiden', method='t-test')
plt.rcParams['figure.figsize'] = [8, 6]
sc.pl.rank_genes_groups(adata_tem, n_genes=10, sharey=False, fontsize=16)
plt.savefig("/mnt/data/lyx/IMC/clustering/leiden_Tumor_markers.png",bbox_inches='tight')

new_cluster_names = ['Ki67+ Tumor cell', 'UTC1','GLUT1+ Tumor cell', 'CAIX+ Tumor cell',
    'GLUT1+ Tumor cell', 'IL7R+ Tumor cell','PRPS1 Tumor cell', 'UTC2','Ki67_CAIX+ Tumor cell','PDL1+ Tumor cell']
adata_tem.obs['MinorType']=[new_cluster_names[int(x)] for x in adata_tem.obs["leiden"].tolist()]

adata_.obs.loc[adata_tem.obs.index,"SubType"]=adata_tem.obs['MinorType']

### Lymphocyte
adata_tem=adata_[adata_.obs["MajorType"]=='Lymphocyte']
sc.pp.neighbors(adata_tem, n_neighbors=30, n_pcs=10) # Computing the neighborhood graph
sc.tl.umap(adata_tem) # Build a UMAP to visualize the neighbourhood graph
sc.tl.leiden(adata_tem)
sc.tl.rank_genes_groups(adata_tem, 'leiden', method='t-test')
plt.rcParams['figure.figsize'] = [8, 6]
sc.pl.rank_genes_groups(adata_tem, n_genes=10, sharey=False, fontsize=16)
plt.savefig("/mnt/data/lyx/IMC/clustering/leiden_Lymphocyte_markers.png",bbox_inches='tight')

new_cluster_names = ['CD8+ T cell', 'B cell','ULC1', 'TNFRSF7+ T cell',
    'CD4+ T cell', 'CD4+ T cell','ULC2', 'CD8+ T cell','NKT cell']
adata_tem.obs['MinorType']=[new_cluster_names[int(x)] for x in adata_tem.obs["leiden"].tolist()]

adata_.obs.loc[adata_tem.obs.index,"SubType"]=adata_tem.obs['MinorType']

### Stromal cell
adata_tem=adata_[adata_.obs["MajorType"]=='Stromal cell']
sc.pp.neighbors(adata_tem, n_neighbors=30, n_pcs=10) # Computing the neighborhood graph
sc.tl.umap(adata_tem) # Build a UMAP to visualize the neighbourhood graph
sc.tl.leiden(adata_tem)
sc.tl.rank_genes_groups(adata_tem, 'leiden', method='t-test')
plt.rcParams['figure.figsize'] = [8, 6]
sc.pl.rank_genes_groups(adata_tem, n_genes=10, sharey=False, fontsize=16)
plt.savefig("/mnt/data/lyx/IMC/clustering/leiden_Stromal_markers.png",bbox_inches='tight')

new_cluster_names = ['CollagenI+ Stromal cell', 'USC1','Vimentin_SMA+ Stromal cell', 'SMA+ Stromal cell',
    'USC2', 'GLUT1+ Stromal cell','USC3', 'FAP+ Stromal cell']
adata_tem.obs['MinorType']=[new_cluster_names[int(x)] for x in adata_tem.obs["leiden"].tolist()]

adata_.obs.loc[adata_tem.obs.index,"SubType"]=adata_tem.obs['MinorType']

### Unkown cell
adata_tem=adata_[adata_.obs["MajorType"]=='Unknown cell']
sc.pp.neighbors(adata_tem, n_neighbors=30, n_pcs=10) # Computing the neighborhood graph
sc.tl.umap(adata_tem) # Build a UMAP to visualize the neighbourhood graph
sc.tl.leiden(adata_tem)
sc.tl.rank_genes_groups(adata_tem, 'leiden', method='t-test')
plt.rcParams['figure.figsize'] = [8, 6]
sc.pl.rank_genes_groups(adata_tem, n_genes=10, sharey=False, fontsize=16)
plt.savefig("/mnt/data/lyx/IMC/clustering/leiden_Unknown_markers.png",bbox_inches='tight')

new_cluster_names = ['CD11c+ cell', 'PRPS1+ cell','UUC1', 'UUC2',
    'PD1+ cell', 'UUC3','UUC4', 'CD80+ cell']
adata_tem.obs['MinorType']=[new_cluster_names[int(x)] for x in adata_tem.obs["leiden"].tolist()]

adata_.obs.loc[adata_tem.obs.index,"SubType"]=adata_tem.obs['MinorType']

adata_.obs["SubType"]=adata_.obs["SubType"].astype('category')

adata_.write("/mnt/data/lyx/IMC/clustering/adata.h5ad")

plt.rcParams['figure.figsize'] = [12, 9]
sc.pl.umap(adata_, color=['SubType'], use_raw=False, s=30,show=False) # Plot the UMAP
plt.savefig("/mnt/data/lyx/IMC/clustering/leiden_clustering_alltype.png",bbox_inches='tight')

## train random-forest classifier
adata=sc.read_h5ad("/mnt/data/lyx/IMC/clustering/leiden/adata.h5ad")
true_Major = adata.obs["MajorType"].tolist()
true_Minor = adata.obs["SubType"].tolist()

x = adata.X
featuresName=adata.var_names

### label encoder
le = preprocessing.LabelEncoder()
le.fit(np.unique(true_Major))
true_Minor_le=le.transform(true_Major)

### splite matrix
train_X, test_X, train_y, test_y = train_test_split(x, true_Minor_le, test_size=0.15, random_state=619)

### train model
rf = RandomForestRegressor(n_estimators=1000)
rf.fit(train_X, train_y)
pred_y = rf.predict(test_X)
pred_y = np.round(pred_y)

### save model
joblib.dump(rf, '/mnt/data/lyx/IMC/clustering/leiden/rfMajormodel.pkl')  
rf = joblib.load('/mnt/data/lyx/IMC/clustering/leiden/rfMajormodel.pkl')

### confusion matrix
def plot_confuse_data(true_y, pred_y):
    confusion = confusion_matrix(y_true=true_y, y_pred=pred_y)

    plt.rcParams['figure.figsize'] = [12, 9]

    #颜色风格为绿。。。。
    plt.imshow(confusion, cmap=plt.cm.Greens)

    confusion = confusion.astype('float') / confusion.sum(axis=1)[:, np.newaxis]# 第一个是迭代对象，表示坐标的显示顺序，第二个参数是坐标轴显示列表

    confusion=confusion[0:5,0:5]

    plt.xlabel('Predicted label')
    plt.ylabel('True label')
    plt.title('Confusion matrix')

# 显示数据
    for i, j in itertools.product(range(confusion.shape[0]), range(confusion.shape[1])):
            plt.text(j, i, '{:.2f}'.format(confusion[i, j]), horizontalalignment="center",
                    color="white" if confusion[i, j] > 0.7 else "black")

    plt.savefig("test.png")

plot_confuse_data(test_y, pred_y)

confusion = confusion_matrix(y_true=test_y, y_pred=pred_y)
confusion = confusion.astype('float') / confusion.sum(axis=1)[:, np.newaxis]
print(np.diagonal(confusion)) 
