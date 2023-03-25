library(SingleCellExperiment)
library(pheatmap)
library(survival)
library(survminer)
source("./spatial_analysis_functions.r")
source("./structural_analysis_functions.r")

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
savePath <- "/mnt/data/lyx/IMC/analysis/spatial/"

## remain Tissue associated Tumor
sce <- sce[, sce$Tissue == "IM"]
# sce <- sce[, sce$Tissue == "IM" | sce$Tissue == "CT"]
scimapResult <- read.csv("/mnt/data/lyx/IMC/analysis/spatial/cellular_neighbor_IM.csv")
scimapResult <- scimapResult[, -1]

colnames(scimapResult)

colName <- c("kmeans_knn_10", "kmeans_knn_20", "kmeans_knn_30")
sce <- BindResult(sce, scimapResult, colName)

reclusMarkers <- c(
    # "CD3", "CD20", "CD4", "CD8a", "FoxP3", "CD57", "CLEC9A", ## Cellidentity
    "CLEC9A",
    "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
    "VEGF", "CAIX", ## Hypoxia
    "CD279", "CD274", "CD366", "TIGIT", "CD127", ## Immune-checkpoint
    "CD27", "CD80" ## Immune-activation
)

ReMajorType <- c("Myeloid")
ReclusterName <- "Myeloid"
ReSubType <- NULL

sce_ <- Reclustering(sce, reclusMarkers, ReMajorType, ReclusterName, ReSubType = ReSubType, ncluster = 15, savePath = paste0(savePath, "knn20_celluarPat/reclustering/"))

CLEC9A_Pos <- c("2")

subtype <- "B"
structure <- "kmeans_3"
cutoffValues <- c(5:25)
cat("The search type is: ", subtype, " structure is: ", structure, " search from ", cutoffValues[1], " to ", cutoffValues[2], "\n")

for (cutoffValue in c(5:25)) {
    TempList <- AssignNewlabel(sce_, allROIs = names(table(sce$ID)), phenoLabel = "CLEC9A_Type", ReclusterName = ReclusterName, interstType = CLEC9A_Pos, cutoffType = "manual", cutoffValue = cutoffValue, numGroup = 3)
    sce_Temp <- TempList[[1]]
    CLEC9Ahigh_ROI <- TempList[[2]]
    CLEC9Alow_ROI <- TempList[[3]]
    CLEC9Anone_ROI <- TempList[[4]]

    ### plot the cell subtype number different
    sce_Temp1 <- sce[, sce$ID %in% CLEC9Ahigh_ROI]
    CLEC9AHighCountMat1 <- GetAbundance(sce_Temp1, countcol = "SubType", is.fraction = T)
    CLEC9AHighCountMat2 <- GetAbundance(sce_Temp1, countcol = "kmeans_knn_20", is.fraction = T)

    sce_Temp2 <- sce[, sce$ID %in% CLEC9Alow_ROI]
    CLEC9AlowCountMat1 <- GetAbundance(sce_Temp2, countcol = "SubType", is.fraction = T)
    CLEC9AlowCountMat2 <- GetAbundance(sce_Temp2, countcol = "kmeans_knn_20", is.fraction = T)

    sce_Temp3 <- sce[, sce$ID %in% CLEC9Anone_ROI]
    CLEC9AnoneCountMat1 <- GetAbundance(sce_Temp3, countcol = "SubType", is.fraction = T)
    CLEC9AnoneCountMat2 <- GetAbundance(sce_Temp3, countcol = "kmeans_knn_20", is.fraction = T)

    BMat1 <- as.numeric(CLEC9AHighCountMat1[, subtype])
    KMat1 <- as.numeric(CLEC9AHighCountMat2[, structure])

    BMat2 <- as.numeric(CLEC9AlowCountMat1[, subtype])
    KMat2 <- as.numeric(CLEC9AlowCountMat2[, structure])

    BMat3 <- as.numeric(CLEC9AnoneCountMat1[, subtype])
    KMat3 <- as.numeric(CLEC9AnoneCountMat2[, structure])

    a1 <- t.test(BMat1, BMat2)$p.value
    b1 <- t.test(BMat2, BMat3)$p.value
    c1 <- t.test(BMat1, BMat3)$p.value

    a2 <- t.test(KMat1, KMat2)$p.value
    b2 <- t.test(KMat2, KMat3)$p.value
    c2 <- t.test(KMat1, KMat3)$p.value
    if ((a1 <= 0.05 | b1 <= 0.05) & (a2 <= 0.05 | b2 <= 0.05)) {
        cat("Cutoff is ", cutoffValue, "\n")
        cat("T-test of ",subtype," cells: ", "High-Low: ", a1, " Low-none: ", b1, " High-none: ", c1, "\n")
        cat("T-test of ",structure," cells: ", "High-Low: ", a2, " Low-none: ", b2, " High-none: ", c2, "\n")
    }
}
