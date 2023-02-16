# continue analysis cDC celltypes

library(RColorBrewer)
library(pheatmap)
library(ComplexHeatmap)
library(factoextra)
library(cluster)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)
library(ggpubr)
library(stringr)
library(reshape2)
library(SingleCellExperiment)

source("/home/lyx/project/IMC/clustering_functions.r")
source("/home/lyx/project/IMC/FlowSOM_metaClustering.r")
source("/home/lyx/project/IMC/abundance_functions.r")
source("/home/lyx/project/IMC/spatial_analysis_functions.r")
source("/home/lyx/project/IMC/structural_analysis_functions.r")

savePath <- "/mnt/data/lyx/IMC/analysis/reclustering/"

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
sce$cellid <- sapply(1:ncol(sce), function(x) {
    return(paste0("cell", as.character(x)))
})
sce_ <- sce[, sce$SubType == "cDC"]
sce_

exp <- t(assays(sce_)[[1]])
dim(exp)
exp[1:5, 1:5]

## calculate gap statistic based on number of clusters
gap_stat <- clusGap(exp, FUN = kmeans, nstart = 25, K.max = 10, B = 50)

## plot number of clusters vs. gap statistic
fviz_gap_stat(gap_stat)

## k-means
fit <- kmeans(exp, centers = 8, nstart = 25)
table(fit$cluster)

sce_$cDCSubtype <- fit$cluster

## plot marker heatmaps
SubtypeHeatmap(sce_, "cDCSubtype", paste0(savePath, "cDC_subtypes_Marker_heatmap.pdf"))

sce_$SubType <- paste0("cDC_", sce_$cDCSubtype)

colData(sce)$SubType[match(colData(sce_)$cellid, colData(sce)$cellid)] <- sce_$SubType
table(sce$SubType)

## abundance analysis
CellCountMat <- Transform_CellCountMat(sce_)

pthreshold <- 0.05
fcthreshold <- 1.5
maxfc <- 5

mat1 <- subset(CellCountMat, Tissue == "IM")
mat1_foldchangeMat <- FCandPvalueCal(mat1, x = c(1, 8))
mat1_foldchangeMat$Foldchange <- ifelse(mat1_foldchangeMat$Foldchange > maxfc, maxfc, mat1_foldchangeMat$Foldchange)
VolcanoPlot(mat1_foldchangeMat, pthreshold = pthreshold, fcthreshold = fcthreshold, filename = paste0(savePath, "IM_cDC_CellFraction_Volcano.pdf"))

mat2 <- subset(CellCountMat, Tissue == "CT")
mat2_foldchangeMat <- FCandPvalueCal(mat2, x = c(1, 8))
mat2_foldchangeMat$Foldchange <- ifelse(mat2_foldchangeMat$Foldchange > maxfc, maxfc, mat2_foldchangeMat$Foldchange)
VolcanoPlot(mat2_foldchangeMat, pthreshold = pthreshold, fcthreshold = fcthreshold, filename = paste0(savePath, "CT_cDC_CellFraction_Volcano.pdf"))

rm(mat1, mat2, mat1_foldchangeMat, mat2_foldchangeMat)

## meta analysis
clinical <- read.csv("/mnt/data/lyx/IMC/clinical.csv")

for (celltype in celltypes2Plot) {
    KMVisualize(IM_plotdf, celltype, cutoff = "mean", savePath = paste0(savePath, "IM_", celltype, "_KMCurve.pdf"))
}

## spatial analysis
sce <- sce[, sce$Tissue == "IM"]
GroupInfo <- load_clinical(sce = sce, clinicalFilePath = "/mnt/data/lyx/IMC/clinical.csv")

## merge celltype abundance into a dataframe
AbundanceDF <- MergeAbundanceResult(sce)

selectCelltypes <- c("B", "CD8T", "NK", "MC5", "cDC_5", "cDC_6", "cDC_1") ## results from abundance analysis volcano plot
celltypes <- names(table(sce$SubType))
celltypes <- celltypes[!celltypes %in% selectCelltypes]

k <- 4
distMat <- dist(t(AbundanceDF[selectCelltypes, ]), method = "euclidean")
colclust <- hclust(distMat, method = "ward.D")
TMEClusterID <- cutree(colclust, k = k)
table(TMEClusterID)

plotDF <- AbundanceDF[c(selectCelltypes, celltypes), ]

annotationCol <- matrix(data = NA, nrow = ncol(plotDF), ncol = 4)
annotationCol <- as.data.frame(annotationCol)

rownames(annotationCol) <- colnames(plotDF)
colnames(annotationCol) <- c("TME Archetypes", "RFSS", "KRAS Mutation", "CRC Site")

annotationCol$`TME Archetypes` <- as.factor(as.numeric(TMEClusterID))
annotationCol$RFSS <- ifelse(GroupInfo$RFS_status == 1, "Relaps", "Non-Relaps")
annotationCol$`KRAS Mutation` <- ifelse(GroupInfo$KRAS_mutation == 1, "Mutate", "WT")
annotationCol$`CRC Site` <- ifelse(GroupInfo$CRC_site == 1, "Right Colon", ifelse(GroupInfo$CRC_site == 2, "Left Colon", "Rectum"))

p <- pheatmap(plotDF,
    scale = "column", gaps_row = length(selectCelltypes), cutree_cols = k,
    annotation_col = annotationCol, annotation_legend = TRUE,
    cluster_rows = FALSE, cluster_cols = colclust, clustering_distance_cols = "euclidean", clustering_method = "complete",
    show_rownames = TRUE, show_colnames = TRUE, fontsize_col = 6, fontsize_row = 6,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    angle_col = "270", cellwidth = 6, cellheight = 6
)
pdf(paste0(savePath, "TME Archetypes.pdf"), height = 5, width = 16)
print(p)
dev.off()
