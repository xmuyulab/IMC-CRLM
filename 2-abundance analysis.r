#### Abundance analysis
library(ggplot2)
source("/home/lyx/project/IMC/abundance_functions.r")

## load data
# clusterResult <- read.csv("/mnt/data/lyx/IMC/prostate/analysis/clustering/annotate_allcells.csv")
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
clinical <- read.table("/mnt/data/lyx/IMC/clinical.csv", sep = ",", header = T)
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

PlotDensityDotplot(sce, marker1 = "CD20", marker2 = "CD3", MajorType = "Lymphocyte", sampleSize = 1000, savePathTemp)
PlotDensityDotplot(sce, marker1 = "EpCAM", marker2 = "CD274", MajorType = "Tumor", sampleSize = 10000, savePathTemp)
PlotDensityDotplot(sce, marker1 = "CD57", marker2 = "CD274", MajorType = "Monocyte", sampleSize = 10000, savePathTemp)
PlotDensityDotplot(sce, marker1 = "CD57", marker2 = "CD279", MajorType = "Monocyte", sampleSize = 10000, savePathTemp)


## downstream analysis
## Other clinical information subgroup abundance analysis
source("/home/lyx/project/IMC/abundance_functions.r")
savePath <- "/mnt/data/lyx/IMC/analysis/abundance/"

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
print(colnames(colData(sce)))

ClassifyCF <- c(
  "Tissue", "Prognosis", "RFS_status", "zs_rec_riskmodel", "fong_score_3", "Gender", "Age_60", "KRAS_mutation",
  "Liver_involvement_num", "CEA_5", "CA199_37", "Pathology", "Differential_grade", "Lymph_positive"
)

## Cell subpopulation counting
countdf <- Transform_CellCountMat(sceobj = sce, group = c("IM", "CT", "TAT"), clinicalFeatures = ClassifyCF, is.fraction = T)
colnames(countdf)

countdf <- countdf[, !(colnames(countdf) %in% "UNKNOWN")] ## Remove Unknown

xCol <- c(1, 20)
clinicalFeatures <- ClassifyCF[-1] ## -1: remove tissue
for (tissue in c("IM", "CT", "TAT")) {
  CountMAT <- abundanceDotPlotMat(countdf, xCol = xCol, clinicalFeatures = clinicalFeatures, tissue = tissue)
  MultiCliDotplot(CountMAT, tissue = tissue, savePath)
}

## Singletype abundance in multiple clincial groups
if (T) {
  ClassifyCF <- c("Tissue", "RFS_status", "Gender", "Age_60", "KRAS_mutation", "TBS_8", "Differential_grade")
  countdf <- Transform_CellCountMat(sce, group = c("IM", "CT", "TAT"), clinicalFeatures = ClassifyCF, is.fraction = T)
  dim(countdf)
  head(countdf)
  CrossBoxplotForAbundance(countdf, celltype = c("CD4T", "CD8T", "Treg"), TypeName = "T cell", clinicalFeatures = ClassifyCF, savePath = savePath)
  CrossBoxplotForAbundance(countdf, celltype = c("B", "CD4T", "CD8T", "NK"), TypeName = "Lymphocyte", clinicalFeatures = ClassifyCF, savePath = savePath)
  CrossBoxplotForAbundance(countdf, celltype = c("SC_aSMA", "SC_COLLAGEN", "SC_FAP", "SC_Vimentin"), TypeName = "Stromal", clinicalFeatures = ClassifyCF, savePath = savePath)
  CrossBoxplotForAbundance(countdf, celltype = c("Macro_CD163", "Macro_CD11b", "Macro_CD169", "Macro_HLADR"), TypeName = "Macrophage", clinicalFeatures = ClassifyCF, savePath = savePath)
  CrossBoxplotForAbundance(countdf, celltype = c("Macro_CD163", "Macro_CD169", "Treg"), TypeName = "Immune Inhibit", clinicalFeatures = ClassifyCF, savePath = savePath)
}


## Abundance Boxplot of relapse
clinicalFeatures <- c("Tissue", "RFS_status")
countdf <- Transform_CellCountMat(sce, group = c("IM", "CT", "TAT"), clinicalFeatures = clinicalFeatures, is.fraction = T)
countdf <- countdf[, !(colnames(countdf) %in% "UNKNOWN")] ## Remove Unknown

# celltypes2Plot , "Macro_CD163", "TC_Ki67")
celltypes2Plot <- colnames(countdf)[xCol[1]:(xCol[2])]

p <- AbundanceBoxPlot(
  sce = sce, countdf = countdf, celltypes2Plot = celltypes2Plot,
  expCol = xCol, tissueCol = "Tissue", clinicalGroupCol = "RFS_status", ClinicalGroupName = c("Relapse", "NonRelapse")
)
pdf("/mnt/data/lyx/IMC/analysis/abundance/relapse_abundance_analysis.pdf", width = 8, height = 6)
print(p)
dev.off()

### Abudance boxplot for selecting types
if (F) {
  # celltypes2Plot <- c("B", "CD4T", "CD8T", "Macro_CD163")
  celltypes2Plot <- names(table(sce$SubType))

  savedir <- "Abundance Boxplot/"
  savePathTemp_ <- paste0("/mnt/data/lyx/IMC/analysis/abundance/", savedir)
  if (!dir.exists(savePathTemp_)) {
    dir.create(savePathTemp_, recursive = T)
  }

  for (i in celltypes2Plot) {
    p <- AbundanceBoxPlot2(
      countdf = countdf, celltypes2Plot = i,
      expCol = xCol, tissueCol = "Tissue", tissue = c("IM", "TAT"), clinicalGroupCol = "RFS_status", ClinicalGroupName = c("Relapse", "NonRelapse")
    )
    pdf(paste0(savePathTemp_, "abudance boxplot of ", i, ".pdf"), width = 4, height = 3)
    print(p)
    dev.off()
  }
  rm(savedir, savePathTemp_)
}

## KM curve
clinicalFeatures <- c("Tissue", "RFS_time", "RFS_status")
countdf <- Transform_CellCountMat(sce, group = c("IM", "CT", "TAT"), clinicalFeatures = clinicalFeatures, is.fraction = T)

# celltypes2Plot <- c("CD4T", "CD8T", "B", "Macro_CD163", "Macro_CD169", "CD3T")
celltypes2Plot <- colnames(countdf)[xCol[1]:xCol[2]]
print(celltypes2Plot)

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

## COX
ContinuousCF <- c(
  "RFS_status", "RFS_time", "fong_score", "Age", "TBS", "CRLM_number", "CRLM_size", "CEA", "CA199", "T_stage"
)
for (tissue in c("IM", "CT", "TAT")) {
  plotdf <- MergeCountsofROI(countdf, tissue = tissue, expCol = xCol, scale = T)
  abundanceMetaAnalysis(plotdf, celltypes2Plot, clinical, features = ContinuousCF, tissue = tissue, savePath)
}


## correlation analysis
for (tissut in c("IM", "CT", "TAT")) {
  plotdf <- subset(countdf, Tissue == tissut)
  plotdf <- plotdf[, c(xCol[1]:xCol[2])]
  CorrelationAnalysis(plotdf, savePath = paste0(savePath, tissut, "_abundance_correlation.pdf"))
}

## Abundance fraction in different tissue
sce <- sce[, sce$Tissue %in% c("CT", "TAT", "IM")]
sce <- sce[, sce$MajorType != "UNKNOWN"]

BarPlotForCelltypeCounts(sce, tissueCol = "Tissue", groupCol = "RFS_status", savePath)
