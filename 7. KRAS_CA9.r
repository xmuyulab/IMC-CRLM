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

for (tissue in c("IM", "CT")) {
    ## Cellular neighborhood pattern
    sce_ <- sce[, sce$Tissue == tissue]

    scimapResult <- read.csv(paste0("/mnt/data/lyx/IMC/analysis/spatial/cellular_neighbor_", tissue, ".csv"))
    scimapResult <- scimapResult[, -1]

    colName <- colnames(scimapResult)[14:20]
    sce_ <- BindResult(sce_, scimapResult, colName)

    structure <- "CNP5"
    if (!file.exists(paste0(savePath, structure, "/", tissue, "/"))) {
        dir.create(paste0(savePath, structure, "/", tissue, "/"))
    }

    ## Cellular pattern difference in Relaps and Non-Relaps
    CompareCellularPattern(sce_, sep = "KRAS_mutation", countcol = structure, n_cluster = 10, savePath = paste0(savePath, structure, "/", tissue, "/"))

    ## Cellular pattenr fraction
    CNP_countsDF <- GetAbundance(sce_, countcol = structure, is.fraction = F, is.reuturnMeans = T)
    CNPFraction(CNP_countsDF, groupBy = "KRAS_mutation", xCol = c(1, 10), savePath = paste0(savePath, structure, "/", tissue, "/"))

    ## seperate ROIs into TC_CA9 high and low groups
    clinicalFeatures <- c("Tissue", "RFS_status", "RFS_time", "KRAS_mutation")
    CA9_countsDF <- GetAbundance(sce_, countcol = "SubType", is.fraction = T, is.reuturnMeans = F)
    CA9_countsDF <- CA9_countsDF[, c("TC_CAIX", "PID", "Tissue", "RFS_status", "RFS_time", "KRAS_mutation")]

    ### grouping
    cutoff <- mean(CA9_countsDF[, 1])
    CA9_countsDF$CA9Label <- ifelse(CA9_countsDF[, 1] > cutoff, "High", "Low")

    celltypes <- names(table(sce_$SubType))
    ResultPath <- paste0("/mnt/data/lyx/IMC/analysis/spatial/permutation_", tissue, "/")

    for (group in c("High", "Low")) {
        CA9_countsDFTemp <- CA9_countsDF[which(CA9_countsDF[, "CA9Label"] == group), ]

        ## Get ROIs
        ROI_KRASMut <- rownames(CA9_countsDFTemp[which(CA9_countsDFTemp$KRAS_mutation == 1), ])
        ROI_KRASWT <- rownames(CA9_countsDFTemp[which(CA9_countsDFTemp$KRAS_mutation == 0), ])

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

        DoubleHeat(MergeDF1, labelDF1, group1 = "KRAS Mut", MergeDF2, labelDF2, group2 = "WT", plot = "circle", savePath = paste0(savePath, "KRAS_spatial_", tissue, " of CA9-", group, ".pdf"))
    }
}
