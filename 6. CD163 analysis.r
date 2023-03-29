# Fyrther investigate CD163+ Macrophage
library(SingleCellExperiment)
library(pheatmap)
library(survival)
library(survminer)
source("./spatial_analysis_functions.r")
source("./structural_analysis_functions.r")

## Load SCE
sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
savePath <- "/mnt/data/lyx/IMC/analysis/Macro_CD163/"
if (!dir.exists(savePath)) {
    dir.create(savePath, recursive = T)
}

## remain Tissue associated Tumor
sce <- sce[, sce$Tissue == "IM"]
# sce <- sce[, sce$Tissue == "IM" | sce$Tissue == "CT"]

## Load CNP
scimapResult <- read.csv("/mnt/data/lyx/IMC/analysis/spatial/cellular_neighbor_IM.csv")
scimapResult <- scimapResult[, -1]

colnames(scimapResult)

colName <- colnames(scimapResult)[14:ncol(scimapResult)]
sce <- BindResult(sce, scimapResult, colName)

## Differential expressiong markers in Macro_CD163
metaMarkers <- c(
    "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
    "VEGF", "CAIX", ## Hypoxia
    "CD279", "CD274", "CD127", ## Immune-checkpoint
    "CD27", "CD80" ## Immune-activation

    # "CD45", "CD20", "CD3", "CD8a", "CD4", "FoxP3", "CD57", ## Lymphocyte
    # "HLADR", "CD68", "CD14", "CD11c", "CD11b", "CD16", "CLEC9A", "CD169", "CD163", ## Myeloid
    # "CollagenI", "Vimentin", "AlphaSMA", "FAP", ## Stromal
    # "EpCAM", ## Tumor
    # "Ki67", "HK2", "FASN", "PRPS1", "GLUT1", ## Cell grouth
    # "VEGF", "CAIX", ## Hypoxia
    # "CD127", "CD27", "CD80", ## Immuno-activation
    # "CD274", "CD279", "TIGIT", "CD366" ## Immuno-checkpoint
)

## plot the differential expression genes
alltypes <- names(table(sce$SubType))

interstMajorTypes <- c(
    "Lymphocyte", "Lymphocyte", "Lymphocyte", "Lymphocyte", "Myeloid",
    "Lymphocyte", "Myeloid", "Myeloid", "Myeloid", "Myeloid",
    "Myeloid", "Myeloid", "Myeloid", "Lymphocyte", "Stromal",
    "Stromal", "Stromal", "Stromal", "Lymphocyte", "Tumor",
    "Tumor", "Tumor", "Tumor", "Lymphocyte"
)
interstSubTypes <- alltypes

## seperate major type
for (i in 1:(length(interstMajorTypes))) {
    interstMajorType <- interstMajorTypes[i]
    interstSubType <- interstSubTypes[i]

    sce_ <- sce[, sce$MajorType %in% interstMajorType]
    mat <- as.data.frame(t(assay(sce_)[metaMarkers, ]))
    mat$phenoLabel <- 0
    mat$phenoLabel[colData(sce_)[, "SubType"] == interstSubType] <- 1
    FCDF <- FCandPvalueCal(mat, xCol = c(1, length(metaMarkers)), yCol = (length(metaMarkers) + 1), need.sample = TRUE)
    VolcanoPlot(FCDF, pthreshold = 0.01, fcthreshold = 1.4, feature = "Phenotype-Associated cells", filename = paste0(savePath, interstSubType, "-associated differential markers in ", interstMajorType, ".pdf"))
}

## Same celltypes different in Relapse and non-relapse
for (i in 1:(length(interstMajorTypes))) {
    interstSubType <- interstSubTypes[i]

    sce_ <- sce[, sce$SubType %in% interstSubType]
    mat <- as.data.frame(t(assay(sce_)[metaMarkers, ]))
    mat$classLabel <- colData(sce_)[, "RFS_status"]
    FCDF <- FCandPvalueCal(mat, xCol = c(1, length(metaMarkers)), yCol = (length(metaMarkers) + 1), need.sample = TRUE)
    VolcanoPlot(FCDF, pthreshold = 0.01, fcthreshold = 1.4, feature = "Phenotype-Associated cells", filename = paste0(savePath, interstSubType, " differential markers in relapse.pdf"))
}

## Assign the new label
### VEGF
VEGF_Pos <- interstType
TempList <- AssignNewlabel(sce_, allROIs = names(table(sce$ID)), phenoLabel = "VEGF_Type", ReclusterName = ReclusterName, interstType = VEGF_Pos, cutoffType = "manual", cutoffValue = 6, numGroup = 2)
sce_Temp <- TempList[[1]]
VEGFhigh_ROI <- TempList[[2]]
VEGFlow_ROI <- TempList[[3]]

ResultPath <- paste0("/mnt/data/lyx/IMC/analysis/spatial/permutation_IM/")
celltypes <- names(table(sce$SubType))
celltypes

### VEGF high ROIs
list_ <- getResult(ResultPath, ROIs = VEGFhigh_ROI, celltypes)
MergeDF1 <- list_[[1]]
labelDF1 <- list_[[2]]
rm(list_)

### VEGF low ROIs
list_ <- getResult(ResultPath, ROIs = VEGFlow_ROI, celltypes)
MergeDF2 <- list_[[1]]
labelDF2 <- list_[[2]]
rm(list_)

DoubleHeat(MergeDF1, labelDF1,
    group1 = "VEGFpos high",
    MergeDF2, labelDF2, group2 = "VEGFpos low", plot = "groupheatmap",
    savePath = paste0(savePath, "/reclustering/", ReclusterName, "_VEGF Interaction Heatmap.pdf")
)

### plot the cell subtype number different
sce_Temp1 <- sce[, sce$ID %in% VEGFhigh_ROI]
VEGFHighCountMat1 <- GetAbundance(sce_Temp1, countcol = "SubType", is.fraction = T)
VEGFHighCountMat2 <- GetAbundance(sce_Temp1, countcol = structure, is.fraction = T)

sce_Temp2 <- sce[, sce$ID %in% VEGFlow_ROI]
VEGFlowCountMat1 <- GetAbundance(sce_Temp2, countcol = "SubType", is.fraction = T)
VEGFlowCountMat2 <- GetAbundance(sce_Temp2, countcol = structure, is.fraction = T)


#### celltype
celltypes <- names(table(sce$SubType))
savePathTemp <- paste0(savePath, "/reclustering/CellAbundance/")
if (!file.exists(savePathTemp)) {
    dir.create(savePathTemp, recursive = T)
}
for (celltype in celltypes) {
    AbundanceSwarmPlot(VEGFHighCountMat1, VEGFlowCountMat1, groupsName = c("high", "low"), celltype = celltype, marker = "VEGF", savePath = savePathTemp, numGroup = 2)
}

#### k-means
k_structures <- names(table(colData(sce)[, structure]))
savePathTemp <- paste0(savePath, "/reclustering/StructureAbundance/")
if (!file.exists(savePathTemp)) {
    dir.create(savePathTemp, recursive = T)
}
for (k_structure in k_structures) {
    AbundanceSwarmPlot(VEGFHighCountMat2, VEGFlowCountMat2, groupsName = c("high", "low"), celltype = k_structure, marker = "VEGF", savePath = savePathTemp, numGroup = 2)
}

#### survival
VEGFCountMat <- GetAbundance(sce_Temp, countcol = "VEGF_Type", is.fraction = T, is.reuturnMeans = T)

SurvivalForPhenoAssoLabel(VEGFCountMat, GroupCol = "Pheno_pos", time = "RFS_time", status = "RFS_status", marker = "VEGF", cutoffType = "median", savePath = paste0(savePath, "/reclustering/"))

