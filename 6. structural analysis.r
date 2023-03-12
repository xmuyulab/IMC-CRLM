# Identify tumor microenviroment pattern
library(SingleCellExperiment)
library(pheatmap)
library(survival)
library(survminer)
source("./spatial_analysis_functions.r")
source("./structural_analysis_functions.r")

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
names(table(sce$Tissue))

## omit Tissue associated Tumor
sce <- sce[, sce$Tissue == "CT" | sce$Tissue == "IM"]
sce <- sce[, sce$Tissue == "IM"]
sce

## load clinical information
GroupInfo <- load_clinical(sce = sce, clinicalFilePath = "/mnt/data/lyx/IMC/clinical.csv")

## TME archetype analysis
if (F) {
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
}

## Cellular neighbors analysis
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

colnames(colData(sce))
table(colData(sce)$MajorType)

## Cell subtype fraction in cellular neighbors pattern
HeatmapForCelltypeInNeighbor(sce, "SubType", "kmeans_knn_20", savePath)

## Cellular pattern difference in Relaps and Non-Relaps
CompareCellularPattern(sce, sep = "RFS_status", countcol = "kmeans_knn_20", n_cluster = 10, savePath = paste0(savePath, "knn20_celluarPat/"))

## Celllular neighborhood pattern survival analysis
CNP_countsDF <- GetAbundance(sce, countcol = "kmeans_knn_20", is.fraction = TRUE, is.reuturnMeans = T)
CNPs <- names(table(colData(sce)$kmeans_knn_20))
for (CNP in CNPs) {
    plotdf <- CNP_countsDF[, c(CNP, "RFS_status", "RFS_time")]
    KMForCNP(plotdf, CNP, savePath = paste0(savePath, "knn20_celluarPat/KM/", "Cellular Neighborhood pattern Suvival analysis of ", CNP, ".pdf"))
}

## Cellular pattenr fraction
CNP_countsDF <- GetAbundance(sce, countcol = "kmeans_knn_20", is.fraction = F, is.reuturnMeans = T)
CNPFraction(CNP_countsDF, groupBy = "RFS_status", xCol = c(1, 10), savePath = paste0(savePath, "knn20_celluarPat/"))

## CNP in image
SavePath1 <- paste0(savePath, "knn20_celluarPat/", "CNP_oncells/")
if (!dir.exists(SavePath1)) {
    dir.create(SavePath1, recursive = T)
}
colData(sce)[, "kmeans_knn_20"] <- as.factor(colData(sce)[, "kmeans_knn_20"])

ROIs <- names(table(colData(sce)$ID))
for (ROI in ROIs) {
    PlotCelltypes(sce, ROI, TypeCol = "kmeans_knn_20", SavePath = paste0(SavePath1, ROI, "_"))
}

## Recluster
### clutering via metabolize molecular
rownames(sce)
metaMarkers <- c("Ki67", "VEGF", "CAIX", "HK2", "FASN", "CD80", "CD274", "PRPS1", "CD279", "GLUT1", "CD27")
ReMajorType <- c("Monocyte")
ReclusterName <- "Monocyte"

sce_ <- Reclustering(sce, metaMarkers, ReMajorType, ReclusterName, ncluster = 7, savePath = paste0(savePath, "knn20_celluarPat/reclustering/"))
# sce_ <- readRDS("/home/lyx/project/IMC/test_sce.rds")

## Certain reclustering types in cellular pattern
interstType <- c("2")
PlotCertainTypeinPattern(sce_, Col1 = ReclusterName, types1 = interstType, Col2 = "kmeans_knn_20", groupCol = "RFS_status", savePath = paste0(savePath, "knn20_celluarPat/reclustering/"))

## assign the new label
if (F) {
    sce_$phenoLabel <- "Pheno_minus"
    sce_$phenoLabel[which(colData(sce_)[, ReclusterName] %in% interstType)] <- "Pheno_plus"
    table(colData(sce_)$phenoLabel)

    phenoLabelCountMat <- GetAbundance(sce_, countcol = "phenoLabel", is.fraction = FALSE, is.reuturnMeans = FALSE)
    phenoLabelCountMat <- TransformIntoPlotMat(phenoLabelCountMat, valueCol = c(1:2))
    head(phenoLabelCountMat)
    BoxPlotForPhenoAssCell(phenoLabelCountMat, savePath)
    SurvivalForPhenoAssCell(phenoLabelCountMat, savePath)
}

## plot the differential expression genes
sce_ <- sce_[, sce_$MajorType == ReMajorType]
mat <- as.data.frame(t(assay(sce_)[metaMarkers, ]))
mat$phenoLabel <- ifelse(sce_$phenoLabel == "Pheno_minus", 0, 1)
FCDF <- FCandPvalueCal(mat, xCol = c(1, 11), yCol = 12, need.sample = TRUE)
VolcanoPlot(FCDF, pthreshold = 0.01, fcthreshold = 3, feature = "Phenotype-Associated cells", filename = paste0(savePath, "Phenotype-associated differential markers in ", ReMajorType, ".pdf"))


## Assign the new label
### CD279 (PD-1)
CD279_Mono <- c("2", "4")
TempList <- AssignNewlabel(sce_, phenoLabel = "CD279_Mono", ReclusterName = ReclusterName, interstType = CD279_Mono, cutoffType = "mean")
sce_ <- TempList[[1]]
CD279high_ROI <- TempList[[2]]
CD279low_ROI <- TempList[[3]]

ResultPath <- paste0("/mnt/data/lyx/IMC/analysis/spatial/permutation_IM/")
celltypes <- names(table(sce$SubType))
celltypes

### CD279 high ROIs
list_ <- getResult(ResultPath, ROIs = CD279high_ROI, celltypes)
MergeDF1 <- list_[[1]]
labelDF1 <- list_[[2]]
rm(list_)

### CD279 low ROIs
list_ <- getResult(ResultPath, ROIs = CD279low_ROI, celltypes)
MergeDF2 <- list_[[1]]
labelDF2 <- list_[[2]]
rm(list_)

DoubleHeat(MergeDF1, labelDF1, group1 = "Mono_CD279 high", MergeDF2, labelDF2, group2 = "Mono_CD279 low", plot = "groupheatmap", savePath = paste0(savePath, "Mono_CD279 Interaction Heatmap.pdf"))

### plot the cell subtype number different
sce_Temp <- sce[, sce$ID %in% CD279high_ROI]
CD279HighCountMat1 <- GetAbundance(sce_Temp, countcol = "SubType", is.fraction = T)
CD279HighCountMat2 <- GetAbundance(sce_Temp, countcol = "kmeans_knn_20", is.fraction = T)

sce_Temp <- sce[, sce$ID %in% CD279low_ROI]
CD279lowCountMat1 <- GetAbundance(sce_Temp, countcol = "SubType", is.fraction = T)
CD279lowCountMat2 <- GetAbundance(sce_Temp, countcol = "kmeans_knn_20", is.fraction = T)

#### celltype
AbundanceSwarmPlot(CD279HighCountMat1, CD279lowCountMat1, groupsName = c("high", "low"), celltype = "Mono_CD57", marker = "CD279", savePath = paste0(savePath, "knn20_celluarPat/reclustering/"))
AbundanceSwarmPlot(CD279HighCountMat1, CD279lowCountMat1, groupsName = c("high", "low"), celltype = "Mono_CLEC9A", marker = "CD279", savePath = paste0(savePath, "knn20_celluarPat/reclustering/"))

#### k-means
AbundanceSwarmPlot(CD279HighCountMat2, CD279lowCountMat2, groupsName = c("high", "low"), celltype = "kmeans_7", marker = "CD279", savePath = paste0(savePath, "knn20_celluarPat/reclustering/"))

#### survival
SurvivalForPhenoAssoLabel(CD279HighCountMat2, CD279lowCountMat2, time = "RFS_time", status = "RFS_status", marker = "CD279", savePath = paste0(savePath, "knn20_celluarPat/reclustering/"))

### CD274 (PD-L1)
CD274_Mono <- c("2")
TempList <- AssignNewlabel(sce_, phenoLabel = "CD274_Mono", ReclusterName = ReclusterName, interstType = CD274_Mono, cutoffType = "mean")
sce_ <- TempList[[1]]
CD274high_ROI <- TempList[[2]]
CD274low_ROI <- TempList[[3]]

### CD274 high ROIs
list_ <- getResult(ResultPath, ROIs = CD274high_ROI, celltypes)
MergeDF1 <- list_[[1]]
labelDF1 <- list_[[2]]
rm(list_)

### CD274 low ROIs
list_ <- getResult(ResultPath, ROIs = CD274low_ROI, celltypes)
MergeDF2 <- list_[[1]]
labelDF2 <- list_[[2]]
rm(list_)

DoubleHeat(MergeDF1, labelDF1, group1 = "Mono_CD274 high", MergeDF2, labelDF2, group2 = "Mono_CD274 low", plot = "groupheatmap", savePath = paste0(savePath, "Mono_CD274 Interaction Heatmap.pdf"))

### plot the cell subtype number different
sce_Temp <- sce[, sce$ID %in% CD274high_ROI]
CD274HighCountMat1 <- GetAbundance(sce_Temp, countcol = "SubType", is.fraction = T)
CD274HighCountMat2 <- GetAbundance(sce_Temp, countcol = "kmeans_knn_20", is.fraction = T)

sce_Temp <- sce[, sce$ID %in% CD274low_ROI]
CD274lowCountMat1 <- GetAbundance(sce_Temp, countcol = "SubType", is.fraction = T)
CD274lowCountMat2 <- GetAbundance(sce_Temp, countcol = "kmeans_knn_20", is.fraction = T)

#### celltype
AbundanceSwarmPlot(CD274HighCountMat1, CD274lowCountMat1, groupsName = c("high", "low"), celltype = "Mono_CD57", marker = "CD274", savePath = paste0(savePath, "knn20_celluarPat/reclustering/"))
AbundanceSwarmPlot(CD274HighCountMat1, CD274lowCountMat1, groupsName = c("high", "low"), celltype = "Mono_CLEC9A", marker = "CD274", savePath = paste0(savePath, "knn20_celluarPat/reclustering/"))

#### k-means
AbundanceSwarmPlot(CD274HighCountMat2, CD274lowCountMat2, groupsName = c("high", "low"), celltype = "kmeans_7", marker = "CD274", savePath = paste0(savePath, "knn20_celluarPat/reclustering/"))

#### survival
SurvivalForPhenoAssoLabel(CD274HighCountMat2, CD274lowCountMat2, time = "RFS_time", status = "RFS_status", marker = "CD274", savePath = paste0(savePath, "knn20_celluarPat/reclustering/"))
