#### Abundance analysis
library(ggplot2)
source("/home/lyx/project/IMC/abundance_functions.r")

## load data
# clusterResult <- read.csv("/mnt/data/lyx/IMC/clustering/qc5_norm_asin_Idenmarke_flowsom_pheno.csv")
clusterResult <- readRDS("/mnt/data/lyx/IMC/analysis/clustering/annotate_allcells.rds")

colnames(clusterResult)
table(clusterResult$SubType)
table(clusterResult$filelist)

## PID
clusterResult$PID <- sapply(clusterResult$filelist, function(x) {
    strsplit(x, "_")[[1]][1]
})

table(clusterResult$PID)

## clinical data
clinical <- read.csv("/mnt/data/lyx/IMC/clinical.csv")
head(clinical, 2)

## QC data (the tissue region and ROI)
qc <- read.csv("/mnt/data/lyx/IMC/IMC_CRC_QC.csv")
head(qc, 2)

## match the tissue region and cell
qc$filelist <- paste0(qc$sample, "_", qc$ROI)
qc$tissueRegion <- NA

for (i in 1:nrow(qc)) {
    if (!is.na(qc[i, 3])) {
        qc$tissueRegion[i] <- "TAT"
        next
    }
    if (!is.na(qc[i, 4])) {
        qc$tissueRegion[i] <- "CT"
        next
    }
    if (!is.na(qc[i, 5])) {
        qc$tissueRegion[i] <- "IM"
        next
    }
    if (!is.na(qc[i, 6])) {
        qc$tissueRegion[i] <- "TLS"
        next
    }
}

qc <- qc[!is.na(qc$QC), ]

clusterResult$Tissue <- qc$tissueRegion[match(clusterResult$ID, qc$filelist)]
table(clusterResult$Tissue)

## transform data into SingleCellExperiment object
print(colnames(clusterResult))
sce <- SCE_Transform(clusterResult, assay_col = c(1, 35), cellmeta_col = c(36, 46), clinical)
saveRDS(sce, paste0("/mnt/data/lyx/IMC/analysis/allsce.rds"))

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
sce <- sce[, sce$Tissue == "IM" | sce$Tissue == "TAT" | sce$Tissue == "CT"]
length(names(table(sce$PID)))

## plot cell marker expression value distribution
savePath <- "/mnt/data/lyx/IMC/analysis/abundance/"
table(sce$Batch)
PlotMarkerExpDistribution(sce, savePath)

## cell type fraction barplot
BarPlotForCelltypeFraction(sce, rowSep = "Tissue", colSep = "RFS_status", savePath)

## Density dotplot
savePathTemp <- "/mnt/data/lyx/IMC/analysis/clustering/"

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
sce <- sce[, sce$Tissue == "IM"]
rownames(sce)

PlotDensityDotplot(sce, marker1 = "CD20", marker2 = "CD3", MajorType = "Lymphocyte", sampleSize = 50000, savePathTemp)
PlotDensityDotplot(sce, marker1 = "EpCAM", marker2 = "CD274", MajorType = "Tumor", sampleSize = 50000, savePathTemp)
PlotDensityDotplot(sce, marker1 = "CD57", marker2 = "CD274", MajorType = "Monocyte", sampleSize = 50000, savePathTemp)
PlotDensityDotplot(sce, marker1 = "CD57", marker2 = "CD279", MajorType = "Monocyte", sampleSize = 50000, savePathTemp)


## downstream analysis
## Other clinical information subgroup abundance analysis
library(SingleCellExperiment)
source("/home/lyx/project/IMC/abundance_functions.r")
savePath <- "/mnt/data/lyx/IMC/analysis/abundance/"

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
print(colnames(colData(sce)))

ClassifyCF <- c(
    "Tissue","Prognosis", "RFS_status", "RFS_time_12", "zs_rec_riskmodel", "fong_score_2", "fong_score_4", "Gender", "Age_60", "KRAS_mutation",
    "TBS_3", "TBS_8", "CRLM_number_4", "CRLM_size_3", "CRLM_size_5", "Liver_involvement_num", "CEA_5","CA199_37","Pathology","Differential_grade","Lymph_positive"
)
countdf <- Transform_CellCountMat(sce, group = c("IM", "CT", "TAT"), clinicalFeatures = ClassifyCF, is.fraction = T)
colnames(countdf)

xCol <- c(1, 27)
clinicalFeatures <- ClassifyCF[-1]  ## -1: remove tissue
for (tissue in c("IM", "CT", "TAT")) {
    CountMAT <- abundanceDotPlotMat(countdf, xCol = xCol, clinicalFeatures = clinicalFeatures, tissue = tissue) 
    MultiCliDotplot(CountMAT, tissue = tissue, savePath)
}


## Abundance Boxplot of relapse
clinicalFeatures <- c("Tissue", "RFS_status")
countdf <- Transform_CellCountMat(sce, group = c("IM", "CT", "TAT"), clinicalFeatures = clinicalFeatures, is.fraction = T)

dim(countdf)
head(countdf)

celltypes2Plot <- c("B", "CD4T", "CD8T", "NK", "Macro_CD163", "Macro_CD11b", "Macro_CD169", "Macro_HLADR", "Mono_CD80")
#celltypes2Plot <- colnames(countdf)[xCol[1]:xCol[2]]
AbunBoxDF <- abundanceBoxplotMat(countdf, celltypes2Plot, MetaCol = c("Tissue", "RFS_status"), expCol = xCol)

AbunBoxDF$RFS_status <- as.factor(ifelse(AbunBoxDF$RFS_status == 1, "Relapse", "Non-Relapse"))

p <- ggplot(AbunBoxDF, aes(x = Tissue, y = Abundance, fill = RFS_status)) +
    geom_boxplot(alpha = 0.7) +
    scale_y_continuous(name = "Cell Abundance") +
    scale_x_discrete(name = "Cell Population") +
    theme_bw() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11, angle = 90)
    ) +
    facet_wrap(~Celltype, ncol = 3) +
    scale_fill_lancet() +
    stat_compare_means(aes(group = RFS_status), label.y = 0.4, method = "t.test")

pdf("/mnt/data/lyx/IMC/analysis/abundance/relapse_abundance_analysis.pdf", width = 15, height = 10)
print(p)
dev.off()

## Abundance Boxplot of KRAS mutation
if (F) {
    clinicalFeatures <- c("Tissue", "KRAS_mutation")
    countdf <- Transform_CellCountMat(sce, group = c("IM", "CT", "TAT"), clinicalFeatures = clinicalFeatures, is.fraction = T)

    dim(countdf)
    head(countdf)

    celltypes2Plot <- c("B", "CD4T", "CD8T", "NK", "Macro_CD163", "Macro_CD11b", "Macro_CD169", "Macro_HLADR", "Mono_CD80")
    AbunBoxDF <- abundanceBoxplotMat(countdf, celltypes2Plot, MetaCol = c("Tissue", "KRAS_mutation"), expCol = xCol)

    AbunBoxDF$KRAS_mutation <- as.factor(ifelse(AbunBoxDF$KRAS_mutation == 1, "KRAS Mutation", "WT"))

    p <- ggplot(AbunBoxDF, aes(x = Tissue, y = Abundance, fill = KRAS_mutation)) +
        geom_boxplot(alpha = 0.7) +
        scale_y_continuous(name = "Cell Abundance") +
        scale_x_discrete(name = "Cell Population") +
        theme_bw() +
        theme(
            plot.title = element_text(size = 14, face = "bold"),
            text = element_text(size = 12),
            axis.title = element_text(face = "bold"),
            axis.text.x = element_text(size = 11, angle = 90)
        ) +
        facet_wrap(~Celltype, ncol = 3) +
        scale_fill_lancet() +
        stat_compare_means(aes(group = KRAS_mutation), label.y = 0.5, method = "t.test")

    pdf("/mnt/data/lyx/IMC/analysis/abundance/kras_abundance_analysis.pdf", width = 15, height = 8)
    print(p)
    dev.off()
}

## KM curve
clinicalFeatures <- c("Tissue", "RFS_time", "RFS_status")
# celltypes2Plot <- c("DPT", "B", "Mono_CLEC9A", "Mono_Classic", "Mono_CD57", "Macro_Multi", "Macro_CD11b", "SC_Vimentin", "TC_Ki67", "TC_TIGHT")
celltypes2Plot <- colnames(countdf)[xCol[1]:xCol[2]]
print(celltypes2Plot)

countdf <- Transform_CellCountMat(sce, group = c("IM", "CT", "TAT"), clinicalFeatures = clinicalFeatures, is.fraction = T)
for (tissue in c("IM", "CT", "TAT")) {
    plotdf <- MergeCountsofROI(countdf, tissue = tissue, expCol = xCol, scale = F)

    for (feature in clinicalFeatures) {
        plotdf[, feature] <- clinical[match(rownames(plotdf), clinical$PID), feature]
    }
    savePathTemp <- paste0(savePath, "KM/", tissue, "/")
    if (!dir.exists(savePathTemp)) {
        dir.create(savePathTemp, recursive = T)
    }
    for (celltype in celltypes2Plot) {
        KMVisualize(plotdf, celltype, cutoff = "best", savePath = paste0(savePathTemp, tissue, "_", celltype, "_KMCurve.pdf"))
    }
}

## meta analysis
ContinuousCF <- c(
    "RFS_status","RFS_time", "fong_score", "Age", "TBS", "CRLM_number", "CRLM_size", "CEA", "CA199", "T_stage"
)
for(tissue in c("IM","CT","TAT")){
    plotdf <- MergeCountsofROI(countdf, tissue = tissue, expCol = xCol, scale = T)
    abundanceMetaAnalysis(plotdf, celltypes2Plot, clinical, features = ContinuousCF, tissue = tissue, savePath)
}


## correlation analysis
for (tissut in c("IM", "CT")) {
    plotdf <- subset(countdf, Tissue == tissut)
    plotdf <- plotdf[, c(xCol[1]:xCol[2])]
    CorrelationAnalysis(plotdf, savePath = paste0(savePath, tissut, "_abundance_correlation.pdf"))
}
