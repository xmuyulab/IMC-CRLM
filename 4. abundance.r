# Abundance analysis
library(ggplot2)
source("/home/lyx/project/IMC/abundance_functions.r")

## load data
# clusterResult <- read.csv("/mnt/data/lyx/IMC/clustering/qc5_norm_asin_Idenmarke_flowsom_pheno.csv")
clusterResult <- readRDS("/mnt/data/lyx/IMC/analysis/clustering/annotate_allcells.rds")

print(colnames(clusterResult))
table(clusterResult$filelist)

### PID
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

### match the tissue region and cell
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

qc$tissueRegion[129] <- "IM"

clusterResult$Tissue <- qc$tissueRegion[match(clusterResult$filelist, qc$filelist)]
table(clusterResult$Tissue)

## transform data into SingleCellExperiment object
sce <- SCE_Transform(clusterResult, assay_col = c(1, 35), cellmeta_col = c(36, 49), clinical)
saveRDS(sce, paste0("/mnt/data/lyx/IMC/analysis/allsce.rds"))

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")

## abundance analysis
countdf <- Transform_CellCountMat(sce, group = c("IM", "CT", "TAT"),is.fraction = F)
abundanceVolcanoPlot(countdf, pthreshold = 0.05, fcthreshold = 1.2)

celltypes2Plot <- c("CD8T", "NK", "B", "Mono_cDC_ITGAX", "Mono_Multi",  "Macro_CD14", "Macro_CD11b", "SC_VEGF","TC_EPCAM")
abundanceBoxplot(countdf, celltypes2Plot, expCol = c(1, 24))

## meta analysis
features <- c(
    "zs_rec_riskmodel", "fong_score", "Age", "Gender", "KRAS_mutation", "BRAF_mutation", "mTBS",
    "CRLM_number", "CRLM_size", "Live_involvement_num", "Pathology", "Differential_grad", "T_grade"
)
IM_plotdf <- MergeCountsofROI(countdf,tissue = "IM", expCol = c(1,24), scale = T)
abundanceMetaAnalysis(IM_plotdf, celltypes2Plot, clinical, features)

## KM curve
IM_plotdf <- subset(countdf, Tissue == "IM")
for (celltype in celltypes2Plot) {
    KMVisualize(IM_plotdf, celltype, cutoff = "median", savePath = paste0("/mnt/data/lyx/IMC/analysis/abundance/IM_", celltype, "_KMCurve.pdf"))
}

## correlation analysis
IM_plotdf <- subset(countdf, Tissue == "IM")
IM_plotdf <- IM_plotdf[, c(1:24)]
CorrelationAnalysis(IM_plotdf, savePath = "/mnt/data/lyx/IMC/analysis/abundance/IM_abundance_correlation.pdf")
