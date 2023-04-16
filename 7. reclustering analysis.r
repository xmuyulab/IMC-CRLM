# Fyrther investigate CD163+ Macrophage
library(SingleCellExperiment)
library(pheatmap)
library(survival)
library(survminer)
source("./spatial_analysis_functions.r")
source("./structural_analysis_functions.r")

## Load SCE
sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
savePath <- "/mnt/data/lyx/IMC/analysis/reclustering/"
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
    "CD274", "CD279", "TIGIT", "CD366", ## Immuno-checkpoint
    "CD127", "CD27", "CD80" ## Immuno-activation
)

## Same celltypes different in Relapse and non-relapse
HeatmapForMarkersInGroups(sce, metaMarkers, groupCol = "RFS_status", adjust.p = T, savePath)

## Treg & Mono_CD11c - Reclustering
anajorTypes <- c("Myeloid", "Lymphocyte")
anaclusterNames <- c("Mono_CD11c", "Treg")
structure <- "CNP20"

for (i in 1:length(anajorTypes)) {
    MajorName <- anajorTypes[i]
    SubName <- anaclusterNames[i]

    savePathTemp1 <- paste0(savePath, SubName, "/")
    if (!dir.exists(savePathTemp1)) {
        dir.create(savePathTemp1, recursive = T)
    }

    ## Reclustering
    sce_ <- Reclustering(sce, metaMarkers, ReMajorType = MajorName, SubName, ReSubType = SubName, PatternCol = structure, ncluster = 15, savePath = savePathTemp1)

    ## Cluster difference
    reclusterAbun <- GetAbundance(sce_, countcol = SubName, is.fraction = T)
    colnames(reclusterAbun)[1:15] <- sapply(colnames(reclusterAbun)[1:15], function(x) {
        return(paste0("rec_", x))
    })
    reclusterAbun$RFS_status <- as.factor(reclusterAbun$RFS_status)
    plotdf <- as.data.frame(matrix(data = NA, nrow = 15 * nrow(reclusterAbun), ncol = 3))
    plotdf[, 1] <- rep(colnames(reclusterAbun)[1:15], each = nrow(reclusterAbun))
    plotdf[, 2] <- as.numeric(as.matrix(reclusterAbun[, 1:15]))
    plotdf[, 3] <- rep(reclusterAbun$RFS_status, times = 15)

    colnames(plotdf) <- c("ReCluster", "Abundance", "Group")

    if (F) {
        plotdf <- plotdf[plotdf$ReCluster %in% c("rec_13", "rec_8", "rec_5", "rec_12"), ]
        rownames(plotdf) <- 1:nrow(plotdf)
    }

    p <- ggplot(plotdf, aes(x = Group, y = Abundance, fill = Group)) +
        geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
        scale_y_continuous(name = "Cell Abundance") +
        scale_x_discrete(name = "Cell Population") +
        scale_fill_manual(values = c("#56B4E9", "#D55E00")) +
        theme_bw() +
        theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            text = element_text(size = 12),
            axis.title = element_text(face = "bold"),
            axis.text.x = element_text(size = 11, angle = 90),
            legend.position = "none",
            strip.text = element_text(size = 12, face = "bold"),
            strip.background = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            panel.grid.major = element_line(color = "gray", size = 0.5),
            panel.grid.minor = element_blank()
        ) +
        stat_compare_means(aes(group = Group), label.y = .4, method = "t.test") +
        facet_wrap(~ReCluster, ncol = 4)

    # pdf(paste0(savePathTemp1, "recluster_abundance_analysis(part).pdf"), width = 8, height = 4)
    pdf(paste0(savePathTemp1, "recluster_abundance_analysis.pdf"), width = 8, height = 6)
    print(p)
    dev.off()

    saveRDS(sce_, paste0(savePathTemp1, SubName, "_recluster.rds"))
}

interstTypeList <- list()
interstTypeList[[1]] <- c("1", "7")
interstTypeList[[2]] <- c("5", "8")

for (i in 1:length(anajorTypes)) {
    MajorName <- anajorTypes[i]
    SubName <- anaclusterNames[i]

    savePathTemp1 <- paste0(savePath, SubName, "/")
    if (!dir.exists(savePathTemp1)) {
        dir.create(savePathTemp1, recursive = T)
    }
    sce_ <- readRDS(paste0("/mnt/data/lyx/IMC/analysis/reclustering/", SubName, "/", SubName, "_recluster.rds"))
    ## Assign the new label
    interstType <- interstTypeList[[i]]

    TempList <- AssignNewlabel(sce_, allROIs = names(table(sce$ID)), phenoLabel = paste0(SubName, "_Type"), ReclusterName = SubName, interstType = interstType, cutoffType = "mean", cutoffValue = 20, numGroup = 2)
    sce_Temp <- TempList[[1]]
    high_ROI <- TempList[[2]]
    low_ROI <- TempList[[3]]

    ## marker difference
    sce__ <- sce_[, sce_$SubType %in% SubName]
    mat <- as.data.frame(t(assay(sce_)[metaMarkers, ]))
    mat$phenoLabel <- 0
    mat$phenoLabel[colData(sce_)[, SubName] %in% interstType] <- 1
    FCDF <- FCandPvalueCal(mat, xCol = c(1, length(metaMarkers)), yCol = (length(metaMarkers) + 1), need.sample = TRUE)
    VolcanoPlot(FCDF, pthreshold = 0.01, fcthreshold = 3, feature = "Phenotype-Associated cells", filename = paste0(savePathTemp1, "Phenotype-associated differential markers in ", SubName, ".pdf"), Qvalue = TRUE)

    ## spatial difference
    ResultPath <- paste0("/mnt/data/lyx/IMC/analysis/spatial/permutation_IM/")
    celltypes <- names(table(sce$SubType))
    celltypes

    ### high ROIs
    list_ <- getResult(ResultPath, ROIs = high_ROI, celltypes)
    MergeDF1 <- list_[[1]]
    labelDF1 <- list_[[2]]
    rm(list_)

    ### low ROIs
    list_ <- getResult(ResultPath, ROIs = low_ROI, celltypes)
    MergeDF2 <- list_[[1]]
    labelDF2 <- list_[[2]]
    rm(list_)

    DoubleHeat(MergeDF1, labelDF1,
        group1 = paste0(SubName, " high"),
        MergeDF2, labelDF2, group2 = paste0(SubName, " low"), plot = "groupheatmap",
        savePath = paste0(savePathTemp1, SubName, " Interaction Heatmap.pdf")
    )
    ### Different
    getInteracDiff(ResultPath, sce, celltypes = celltypes, savepath = paste0(savePathTemp1, SubName, "_Spatial_Diff_.pdf"), IDs1 = high_ROI, IDs2 = low_ROI)

    ### plot the cell subtype number different
    sce_Temp1 <- sce[, sce$ID %in% high_ROI]
    HighCountMat1 <- GetAbundance(sce_Temp1, countcol = "SubType", is.fraction = T)
    HighCountMat2 <- GetAbundance(sce_Temp1, countcol = structure, is.fraction = T)

    sce_Temp2 <- sce[, sce$ID %in% low_ROI]
    lowCountMat1 <- GetAbundance(sce_Temp2, countcol = "SubType", is.fraction = T)
    lowCountMat2 <- GetAbundance(sce_Temp2, countcol = structure, is.fraction = T)


    #### celltype
    celltypes <- names(table(sce$SubType))
    savePathTemp <- paste0(savePathTemp1, "/CellAbundance/")
    if (!file.exists(savePathTemp)) {
        dir.create(savePathTemp, recursive = T)
    }
    for (celltype in celltypes) {
        AbundanceSwarmPlot(HighCountMat1, lowCountMat1, groupsName = c("high", "low"), celltype = celltype, marker = SubName, savePath = savePathTemp, numGroup = 2)
    }

    #### Multiple cell subpopulations abudnace plot
    PlotTypes <- c("CD8T", "Macro_HLADR", "Macro_CD11b", "Macro_CD163", "Macro_CD169", "Treg")
    p <- MultipleAbundanceSwarmPlot(HighCountMat1, lowCountMat1, groupsName = c("high", "low"), PlotTypes = PlotTypes, marker = SubName, numGroup = 2, style = "box")
    pdf(paste0(savePathTemp1, "ViolinPlot of multiple celltypes in ", marker, " group.pdf"), height = 10, width = 8)
    print(p)
    dev.off()

    #### CNPs
    k_structures <- names(table(colData(sce)[, structure]))
    savePathTemp <- paste0(savePathTemp1, "/StructureAbundance/")
    if (!file.exists(savePathTemp)) {
        dir.create(savePathTemp, recursive = T)
    }
    for (k_structure in k_structures) {
        AbundanceSwarmPlot(HighCountMat2, lowCountMat2, groupsName = c("high", "low"), celltype = k_structure, marker = SubName, savePath = savePathTemp, numGroup = 2)
    }

    #### survival
    CountMat <- GetAbundance(sce_Temp, countcol = paste0(SubName, "_Type"), is.fraction = T, is.reuturnMeans = F)

    SurvivalForPhenoAssoLabel(CountMat, GroupCol = "Pheno_pos", time = "RFS_time", status = "RFS_status", marker = SubName, cutoffType = "best", savePath = savePathTemp1)
}

## CD163+ Macrophage analysis
SubName <- "Macro_CD163"
interstType <- "Macro_CD163"

savePathTemp1 <- paste0(savePath, SubName, "/")
if (!dir.exists(savePathTemp1)) {
    dir.create(savePathTemp1, recursive = T)
}

## Assign the new label
TempList <- AssignNewlabel(sce, allROIs = names(table(sce$ID)), phenoLabel = paste0(SubName, "_Label"), ReclusterName = "SubType", interstType = interstType, cutoffType = "mean", cutoffValue = 20, numGroup = 2)
sceTemp <- TempList[[1]]
high_ROI <- TempList[[2]]
low_ROI <- TempList[[3]]

## spatial difference
ResultPath <- paste0("/mnt/data/lyx/IMC/analysis/spatial/permutation_IM/")
celltypes <- names(table(sce$SubType))
celltypes

### high ROIs
list_ <- getResult(ResultPath, ROIs = high_ROI, celltypes)
MergeDF1 <- list_[[1]]
labelDF1 <- list_[[2]]
rm(list_)

### low ROIs
list_ <- getResult(ResultPath, ROIs = low_ROI, celltypes)
MergeDF2 <- list_[[1]]
labelDF2 <- list_[[2]]
rm(list_)

DoubleHeat(MergeDF1, labelDF1,
    group1 = paste0(SubName, " high"),
    MergeDF2, labelDF2, group2 = paste0(SubName, " low"), plot = "groupheatmap",
    savePath = paste0(savePathTemp1, SubName, " Interaction Heatmap.pdf")
)
### Different
getInteracDiff(ResultPath, sce, celltypes = celltypes, savepath = paste0(savePathTemp1, SubName, "_Spatial_Diff_0.01.pdf"), IDs1 = high_ROI, IDs2 = low_ROI, Sig = 0.01)

PlotTypes <- c("B", "CD4T", "CD8T", "Treg", "Macro_CD163", "Macro_CD169", "Mono_Classic", "Mono_Intermediate")
p <- InterDiffInvertBarplot(
    ResultPath = ResultPath, sce = sce, PlotTypes = PlotTypes,
    IDs1 = high_ROI, IDs2 = low_ROI
)
pdf(paste0(savePathTemp1, SubName, "_Spatial_Diff_Barplot.pdf"), height = 6, width = 8)
print(p)
dev.off()

p <- InterAbundanceBarplot(
    ResultPath = ResultPath, sce = sce, PlotTypes = PlotTypes, groupsName = c("High", "Low"), IDs1 = high_ROI, IDs2 = low_ROI
)
pdf(paste0(savePathTemp1, SubName, "_SpatialAbundance.pdf"), height = 6, width = 10)
print(p)
dev.off()

### plot the cell subtype number different
sceTemp1 <- sce[, sce$ID %in% high_ROI]
HighCountMat1 <- GetAbundance(sceTemp1, countcol = "SubType", is.fraction = T)
HighCountMat2 <- GetAbundance(sceTemp1, countcol = structure, is.fraction = T)

sceTemp2 <- sce[, sce$ID %in% low_ROI]
lowCountMat1 <- GetAbundance(sceTemp2, countcol = "SubType", is.fraction = T)
lowCountMat2 <- GetAbundance(sceTemp2, countcol = structure, is.fraction = T)


#### celltype
celltypes <- names(table(sce$SubType))
savePathTemp <- paste0(savePathTemp1, "CellAbundance/")
if (!file.exists(savePathTemp)) {
    dir.create(savePathTemp, recursive = T)
}
for (celltype in celltypes) {
    AbundanceSwarmPlot(HighCountMat1, lowCountMat1, groupsName = c("high", "low"), celltype = celltype, marker = SubName, savePath = savePathTemp, numGroup = 2)
}

#### Multiple cell subpopulations abudnace plot
PlotTypes <- c("B", "CD8T", "Macro_CD11b", "Macro_CD163", "Macro_CD169", "Treg")
p <- MultipleAbundanceSwarmPlot(HighCountMat1, lowCountMat1, groupsName = c("high", "low"), PlotTypes = PlotTypes, marker = SubName, numGroup = 2, style = "violin")
pdf(paste0(savePathTemp1, "ViolinPlot of multiple celltypes in ", marker, " group.pdf"), height = 10, width = 8)
print(p)
dev.off()

#### CNPs
k_structures <- names(table(colData(sce)[, structure]))
savePathTemp <- paste0(savePathTemp1, "StructureAbundance/")
if (!file.exists(savePathTemp)) {
    dir.create(savePathTemp, recursive = T)
}
for (k_structure in k_structures) {
    AbundanceSwarmPlot(HighCountMat2, lowCountMat2, groupsName = c("high", "low"), celltype = k_structure, marker = SubName, savePath = savePathTemp, numGroup = 2)
}

#### survival
CountMat <- GetAbundance(sceTemp, countcol = paste0(SubName, "_Label"), is.fraction = T, is.reuturnMeans = F)

SurvivalForPhenoAssoLabel(CountMat, GroupCol = "Pheno_pos", time = "RFS_time", status = "RFS_status", marker = SubName, cutoffType = "best", savePath = savePathTemp1)
