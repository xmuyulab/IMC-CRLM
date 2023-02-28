# Identify tumor microenviroment archetypes
library(SingleCellExperiment)
library(pheatmap)
library(survival)
library(survminer)
source("./structural_analysis_functions.r")
source("./spatial_analysis_functions.r")

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
names(table(sce$Tissue))

## omit Tissue associated Tumor
sce <- sce[, sce$Tissue == "CT" | sce$Tissue == "IM"]
sce <- sce[, sce$Tissue == "IM"]
sce

## load clinical information
GroupInfo <- load_clinical(sce = sce, clinicalFilePath = "/mnt/data/lyx/IMC/clinical.csv")

## merge celltype abundance into a dataframe
AbundanceDF <- MergeAbundanceResult(sce, return.fraction = T)

selectCelltypes <- c("DPT", "Mono_CLEC9A", "Macro_Multi", "Macro_CD11b", "SC_Vimentin", "TC_Ki67") ## results from abundance analysis volcano plot

celltypes <- names(table(sce$SubType))
celltypes <- celltypes[!celltypes %in% selectCelltypes]

## cluster patients
distMat <- dist(t(AbundanceDF[selectCelltypes, ]), method = "euclidean")
colclust <- hclust(distMat, method = "complete")

k <- 15
TMEClusterID <- cutree(colclust, k = k)
table(TMEClusterID)

## column annotation bar
plotDF <- AbundanceDF[c(selectCelltypes, celltypes), ]

annotationCol <- matrix(data = NA, nrow = ncol(plotDF), ncol = 3)
annotationCol <- as.data.frame(annotationCol)

rownames(annotationCol) <- colnames(plotDF)
colnames(annotationCol) <- c("TME Archetypes", "RFSS", "KRAS Mutation")

annotationCol$`TME Archetypes` <- as.factor(as.numeric(TMEClusterID))
annotationCol$RFSS <- ifelse(GroupInfo$RFS_status == 1, "Relaps", "Non-Relaps")
annotationCol$`KRAS Mutation` <- ifelse(GroupInfo$KRAS_mutation == 1, "Mutate", "WT")
# annotationCol$`CRC Site` <- ifelse(GroupInfo$CRC_site == 1, "Right Colon", ifelse(GroupInfo$CRC_site == 2, "Left Colon", "Rectum"))

p <- pheatmap(plotDF,
    scale = "column", gaps_row = length(selectCelltypes), cutree_cols = k,
    annotation_col = annotationCol, annotation_legend = TRUE,
    cluster_rows = FALSE, cluster_cols = colclust, clustering_distance_cols = "euclidean", clustering_method = "complete",
    show_rownames = TRUE, show_colnames = TRUE, fontsize_col = 6, fontsize_row = 6,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    angle_col = "270", cellwidth = 6, cellheight = 6
)
pdf("/mnt/data/lyx/IMC/analysis/spatial/TME Archetypes.pdf", height = 5, width = 16)
print(p)
dev.off()

## Assign TME Archetypes label (Voting)
PlaeblDF <- as.data.frame(TMEClusterID)
PlaeblDF$PID <- sapply(rownames(PlaeblDF), function(x) {
    strsplit(x, "_")[[1]][1]
})
PIDs <- names(table(PlaeblDF$PID))

label <- c()
for (i in PIDs) {
    PlaeblDFTemp <- subset(PlaeblDF, PID == i)
    label <- c(label, names(table(PlaeblDFTemp$TMEClusterID))[1])
}
names(label) <- PIDs

## survival analysis
label <- as.data.frame(label)
label$RFS_time <- GroupInfo$RFS_time[match(rownames(label), GroupInfo$PID)]
label$RFS_status <- GroupInfo$RFS_status[match(rownames(label), GroupInfo$PID)]
label <- label[-match("W21", rownames(label)), ]

km <- survfit(Surv(RFS_time, RFS_status) ~ label, data = label)
p <- ggsurvplot(km,
    data = label,
    linetype = c("solid", "solid"),
    surv.median.line = "hv", surv.scale = "percent",
    pval = T, risk.table = T,
    conf.int = T, conf.int.alpha = 0.1, conf.int.style = "ribbon",
    risk.table.y.text = T,
    palette = c("#3300CC", "#CC3300"),
    xlab = "Recurrence time"
)

pdf("/mnt/data/lyx/IMC/analysis/spatial/TME archetypes survival analysis.pdf", width = 8, height = 6)
print(p)
dev.off()

## visualize TME archetypes
savePath <- "/mnt/data/lyx/IMC/analysis/spatial/archetypes/"
if (!dir.exists(savePath)) {
    dir.create(savePath)
}
ROIs <- names(table(sce$ID))
for (ROI in ROIs) {
    PlotCelltypes(sce, ROI, selectCelltypes, SavePath = paste0(savePath, ROI, " TME archetypes.pdf"))
}





# Cellular neighbors analysis
sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
savePath <- "/mnt/data/lyx/IMC/analysis/spatial/"

## remain Tissue associated Tumor
sce <- sce[, sce$Tissue == "IM"]
# sce <- sce[, sce$Tissue == "IM" | sce$Tissue == "CT"]
scimapResult <- read.csv("/mnt/data/lyx/IMC/analysis/spatial/cellular_neighbor.csv")
scimapResult <- scimapResult[, -1]

colnames(scimapResult)

colName <- c("kmeans_knn_10", "kmeans_knn_20", "kmeans_knn_30")
sce <- BindResult(sce, scimapResult, colName)

colnames(colData(sce))
table(colData(sce)$MajorType)

## Cell subtype fraction in cellular neighbors pattern
HeatmapForCelltypeInNeighbor(sce, "SubType", "kmeans_knn_20", savePath)

## Cellular pattern difference in Relaps and Non-Relaps
CompareCellularPattern(sce, sep = "RFS_status", countcol = "kmeans_knn_20", n_cluster = 10, savePath = savePath)

## Recluster
### clutering via metabolize molecular
rownames(sce)
metaMarkers <- c("Ki67", "VEGF", "CAIX", "HK2", "FASN", "CD80", "CD274", "PRPS1", "CD279", "GLUT1", "CD27")
ReMajorType <- c("Monocyte")
ReclusterName <- "Monocyte"

sce_ <- Reclustering(sce, metaMarkers, ReMajorType, ReclusterName, ncluster = 10, savePath)
# sce_ <- readRDS("/home/lyx/project/IMC/test_sce.rds")

## Certain reclustering types in cellular pattern
interstType <- c("9")
PlotCertainTypeinPattern(sce_, Col1 = ReclusterName, types1 = interstType, Col2 = "kmeans_knn_20", groupCol = "RFS_status", savePath)

## assign the new label
sce_$phenoLabel <- "Pheno_minus"
sce_$phenoLabel[which(colData(sce_)[, ReclusterName] %in% interstType)] <- "Pheno_plus"
table(colData(sce_)$phenoLabel)

phenoLabelCountMat <- GetAbundance(sce_, countcol = "phenoLabel", is.fraction = FALSE, is.reuturnMeans = FALSE)
phenoLabelCountMat <- TransformIntoPlotMat(phenoLabelCountMat, valueCol = c(1:2))
head(phenoLabelCountMat)
BoxPlotForPhenoAssCell(phenoLabelCountMat, savePath)
SurvivalForPhenoAssCell(phenoLabelCountMat, savePath)


## plot the differential expression genes
sce_ <- sce_[, sce_$MajorType == ReMajorType]
mat <- as.data.frame(t(assay(sce_)[metaMarkers, ]))
mat$phenoLabel <- ifelse(sce_$phenoLabel == "Pheno_minus", 0, 1)
FCDF <- FCandPvalueCal(mat, xCol = c(1, 11), yCol = 12, need.sample = TRUE)
VolcanoPlot(FCDF, pthreshold = 0.01, fcthreshold = 3, feature = "Phenotype-Associated cells", filename = paste0(savePath, "Phenotype-associated differential markers in ", ReMajorType, ".pdf"))
