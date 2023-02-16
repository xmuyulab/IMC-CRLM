# functions for spatial analysis
import os
import numpy as np
import pandas as pd
from anndata import AnnData
from sklearn.metrics.pairwise import euclidean_distances
from collections import Counter
from random import shuffle

def MergeClinical(expPath, clinicalPath, qcPath, metacolStart):
    # load data
    exp = pd.read_csv(expPath)
    clinical = pd.read_csv(clinicalPath)
    qc = pd.read_csv(qcPath)

    PID = [x.split("_")[0] for x in exp["ID"].tolist()]
    exp["PID"] = PID

    clinical.index = clinical["PID"]

    imageid=[]
    tistype=[]
    for i in range(qc.shape[0]):
        imageid.append(qc["sample"][i]+'_'+qc["ROI"][i])
        if qc.iloc[i,2] == 1:
            tistype.append("TAT")
        elif qc.iloc[i,3] == 1:
            tistype.append("CT")
        elif qc.iloc[i,4] == 1:
            tistype.append("IM")
        elif qc.iloc[i,5] == 1:
            tistype.append("TLS")
        else:
            tistype.append("NA")
    tistype[128] = "IM"

    qc.index = imageid
    qc["Tissue"]=tistype

    # merge recurrence information
    Rec = [str(clinical.loc[str(x), "RFS_status"]) for x in PID]
    exp["RFS_status"] = Rec

    # Tissue
    Tis = [str(qc.loc[str(x), "Tissue"]) for x in exp["ID"]]
    exp["Tissue"] = Tis

    # create anndata
    meta = exp.iloc[:, metacolStart:]
    exp = exp.iloc[:, :metacolStart]

    var = pd.DataFrame()
    var["Marker"] = exp.columns

    X = np.array(exp)

    adata = AnnData(X=X, obs=meta, var=var)

    # spatial information
    x = [str.split(a,',')[0] for a in meta["Position"]]
    y = [str.split(a,',')[1] for a in meta["Position"]]

    x_ = [str.split(a,'(')[1] for a in x]
    y_ = [str.split(a,')')[0] for a in y]

    spatial = pd.DataFrame()
    spatial["x"] = [float(a) for a in x_]
    spatial["y"] = [float(a) for a in y_]

    adata.uns["spatial"] = spatial

    del adata.obs["Position"]
    adata.obs["x"] = [float(a) for a in x_]
    adata.obs["y"] = [float(a) for a in y_]

    return adata

## Permutation test
### calculate the cell distance
def GetDistance(adata):
    x = adata.obs["x"]
    y = adata.obs["y"]

    coordinates = np.array([x,y]).T

    disMat = euclidean_distances(coordinates,coordinates)

    return disMat

### calculate neighbors, distance less than threshold
def GetNeighbors(disMat,dis_threshold):
    disMat2 = np.where(disMat<=dis_threshold,1,0)
    disMat2 = np.triu(disMat2) - np.eye(disMat2.shape[0])
    disMat2 = disMat2.astype(np.int0)

    cell1 = []
    cell2 = []

    nrow = disMat2.shape[0]
    ncol = disMat2.shape[1]

    for i in range(nrow):
        a = disMat2[i,:]

        cell1Temp = [i for _ in range(a.sum())]
        cell1.extend(cell1Temp)

        for j in range(i+1,ncol):
            if a[j] == 1:
                cell2.append(j)
    
    neiMat = pd.DataFrame()
    neiMat["cell1"] = cell1
    neiMat["cell2"] = cell2

    return neiMat

def Permutation(neiMat,celllabel,perm_n=1001):
    ### permutation
    alllabel = [x for x in Counter(celllabel).keys()]
    numlabel = len(alllabel)
    numPair = neiMat.shape[0]

    permuLabel = celllabel.copy()

    PermutationResut = dict()

    for i in range(perm_n):
        if i == 0:
            permuLabel = celllabel.copy()
        else:
            shuffle(permuLabel)

        print("Perform permutation test: "+str(i))
        interactNumMat = np.zeros(shape=(numlabel,numlabel)).astype(np.int0)
        interactNumMat = pd.DataFrame(interactNumMat)
        interactNumMat.columns = alllabel
        interactNumMat.index = alllabel

        cell1Temp = [permuLabel[x] for x in neiMat["cell1"]]
        cell2Temp = [permuLabel[x] for x in neiMat["cell2"]]
        
        pairTemp = [cell1Temp[a]+'__'+cell2Temp[a] for a in range(numPair)]
        pairCount = Counter(pairTemp)

        for j in pairCount.keys():
            c1,c2 = j.split("__")
            interactNumMat.loc[c1,c2] += pairCount[j]

        diag_ = np.tril(interactNumMat)-np.tril(interactNumMat,-1)
        interactNumMat_ = interactNumMat + np.transpose(interactNumMat) - diag_

        PermutationResut["permutation"+str(i)] = interactNumMat_

    return PermutationResut

def PermutationTest(PermutationResut,Path):
    numlabel = PermutationResut["permutation0"].shape[0]

    result = np.zeros(shape=(numlabel,numlabel))
    result = pd.DataFrame(result)
    result.columns = PermutationResut["permutation0"].columns
    result.index = PermutationResut["permutation0"].index

    clonsness = result.copy()    
    avoidance = result.copy()

    trueInteraction = PermutationResut["permutation0"]

    permuNum = len(PermutationResut)-1

    for i in range(numlabel):
        for j in range(numlabel):

            intergdTmep = trueInteraction.iloc[i,j]
            if intergdTmep==0:
                clonsness.iloc[i,j]=1
                avoidance.iloc[i,j]=1
                continue;

            permuList = [PermutationResut["permutation"+str(z)].iloc[i,j] for z in range(1,len(PermutationResut))]
            cloTemp = [1 if intergdTmep<permuList[z] else 0 for z in range(len(permuList))]
            avoiTemp = [1 if intergdTmep>permuList[z] else 0 for z in range(len(permuList))]

            clonsness.iloc[i,j] = (sum(cloTemp)+1)/(permuNum+1)
            avoidance.iloc[i,j] = (sum(avoiTemp)+1)/(permuNum+1)
    
    clonsness.to_csv(os.path.join(Path,"ClosePvalue","ClosePvalue.csv"))
    avoidance.to_csv(os.path.join(Path,"AvoidPvalue","AvoidPvalue.csv"))
    trueInteraction.to_csv(os.path.join(Path,"TureInteraction","TureInteraction.csv"))

    return

def PermutationTestMain(adata,ID,dis_threshold,Path):
    adata = adata[adata.obs["ID"]==ID]
    print("Permutation test on sample "+ ID)

    if not os.path.exists(Path):
        os.mkdir(Path)
        os.mkdir(os.path.join(Path,"TureInteraction"))
        os.mkdir(os.path.join(Path,"ClosePvalue"))
        os.mkdir(os.path.join(Path,"AvoidPvalue"))

    disMat = GetDistance(adata)
    neiMat = GetNeighbors(disMat,dis_threshold)
    PermutationResut = Permutation(neiMat,celllabel = adata.obs["SubType"])
    PermutationTest(PermutationResut,Path)

    return

