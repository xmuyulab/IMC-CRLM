# TC_CA9 and KRAS mutation

library(RColorBrewer)
library(ggplot2)
library(SingleCellExperiment)

source("/home/lyx/project/IMC/clustering_functions.r")
source("/home/lyx/project/IMC/FlowSOM_metaClustering.r")
source("/home/lyx/project/IMC/abundance_functions.r")
source("/home/lyx/project/IMC/spatial_analysis_functions.r")
source("/home/lyx/project/IMC/structural_analysis_functions.r")

savePath <- "/mnt/data/lyx/IMC/analysis/KRAS/"
if (!dir.exists(savePath)) {
    dir.create(savePath)
}

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
sce <- sce[, (sce$Tissue == "IM" | sce$Tissue == "CT")]
sce

clinical <- read.csv("/mnt/data/lyx/IMC/clinical.csv")

## spatial interation
GroupInfo <- GetGroupInfo(sce, clinical)
celltypes <- names(table(sce$SubType))

for (tissue in c("IM", "CT")) {
    ResultPath <- paste0("/mnt/data/lyx/IMC/analysis/spatial/permutation_", tissue, "/")

    sceTemp <- sce[, (sce$Tissue == tissue)]

    ## Get ROIs
    ROI_KRASMut <- names(table(sceTemp[, sceTemp$KRAS_mutation == 1]$ID))
    ROI_KRASWT <- names(table(sceTemp[, sceTemp$KRAS_mutation == 0]$ID))

    ### KRAS mutation
    list_ <- getResult(ResultPath, ROIs = ROI_KRASMut, celltypes)
    MergeDF1 <- list_[[1]]
    labelDF1 <- list_[[2]]
    rm(list_)

    ### WT
    list_ <- getResult(ResultPath, ROIs = ROI_KRASWT, celltypes)
    MergeDF2 <- list_[[1]]
    labelDF2 <- list_[[2]]
    rm(list_)

    DoubleHeat(MergeDF1, labelDF1, group1 = "KRAS Mut", MergeDF2, labelDF2, group2 = "WT", plot = "groupheatmap", savePath = paste0(savePath, "KRAS_spatial_", tissue, ".pdf"))

    ### Different
    array1 <- LoadAnalysisResult(ROI_KRASMut, celltypes, sceTemp)
    array2 <- LoadAnalysisResult(ROI_KRASWT, celltypes, sceTemp)

    TestDiff(array1, array2, celltypes, savepath = paste0(savePath, "KRAS_Spatial_Diff_", tissue, ".pdf"))
}

## Cellular neighborhood pattern
sce <- sce[, sce$Tissue == "IM"]

scimapResult <- read.csv("/mnt/data/lyx/IMC/analysis/spatial/cellular_neighbor_IM.csv")
scimapResult <- scimapResult[, -1]

colName <- c("kmeans_knn_10", "kmeans_knn_20", "kmeans_knn_30")
sce <- BindResult(sce, scimapResult, colName)

## Cellular pattern difference in Relaps and Non-Relaps
CompareCellularPattern(sce, sep = "KRAS_mutation", countcol = "kmeans_knn_20", n_cluster = 10, savePath = paste0(savePath, "knn20_celluarPat/"))

## Cellular pattenr fraction
CNP_countsDF <- GetAbundance(sce, countcol = "kmeans_knn_20", is.fraction = F, is.reuturnMeans = T)
CNPFraction(CNP_countsDF, groupBy = "KRAS_mutation", xCol = c(1, 10), savePath = paste0(savePath, "knn20_celluarPat/"))

## seperate ROIs into TC_CA9 high and low groups
clinicalFeatures <- c("Tissue", "RFS_status", "RFS_time", "KRAS_mutation")
CA9_countsDF <- GetAbundance(sce, countcol = "SubType", is.fraction = T, is.reuturnMeans = F)
CA9_countsDF <- CA9_countsDF[,c("TC_CAIX","PID","Tissue","RFS_status","RFS_time","KRAS_mutation")]

### grouping 
mean_ <- mean(CA9_countsDF[,1])
CA9_countsDF$CA9Label <- ifelse(CA9_countsDF[, 1] > mean_, "High", "Low")

tissue <- "IM"
celltypes <- names(table(sce$SubType))
ResultPath <- paste0("/mnt/data/lyx/IMC/analysis/spatial/permutation_", tissue, "/")

for (group in c("High", "Low")) {

    CA9_countsDFTemp <- CA9_countsDF[which(CA9_countsDF[,"CA9Label"]==group),]

    ## Get ROIs
    ROI_KRASMut <- rownames(CA9_countsDFTemp[which(CA9_countsDFTemp$KRAS_mutation==1),])
    ROI_KRASWT <- rownames(CA9_countsDFTemp[which(CA9_countsDFTemp$KRAS_mutation==0),])

    ### KRAS mutation
    list_ <- getResult(ResultPath, ROIs = ROI_KRASMut, celltypes)
    MergeDF1 <- list_[[1]]
    labelDF1 <- list_[[2]]
    rm(list_)

    ### WT
    list_ <- getResult(ResultPath, ROIs = ROI_KRASWT, celltypes)
    MergeDF2 <- list_[[1]]
    labelDF2 <- list_[[2]]
    rm(list_)

    DoubleHeat(MergeDF1, labelDF1, group1 = "KRAS Mut", MergeDF2, labelDF2, group2 = "WT", plot = "circle", savePath = paste0(savePath, "KRAS_spatial_", tissue, " of CA9-",group,".pdf"))
}
