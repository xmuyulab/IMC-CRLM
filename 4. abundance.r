#### Abundance analysis
library(ggplot2)
source("/home/lyx/project/IMC/abundance_functions.r")

## load data
#clusterResult <- read.csv("/mnt/data/lyx/IMC/clustering/qc5_norm_asin_Idenmarke_flowsom_pheno.csv")
clusterResult <- readRDS("/mnt/data/lyx/IMC/analysis/clustering/annotate_allcells.rds")

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
sce <- SCE_Transform(clusterResult, assay_col = c(1, 35), cellmeta_col = c(36, 44), clinical)
saveRDS(sce, paste0("/mnt/data/lyx/IMC/analysis/allsce.rds"))

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")

## plot cell marker expression value distribution
savePath <- "/mnt/data/lyx/IMC/analysis/abundance/"
table(sce$Batch)
PlotMarkerExpDistribution(sce, savePath)

## cell type fraction barplot
BarPlotForCelltypeFraction(sce, rowSep = "Tissue", colSep = "RFS_status", savePath)

## downstream analysis
## Other clinical information subgroup abundance analysis
library(SingleCellExperiment)
source("/home/lyx/project/IMC/abundance_functions.r")
savePath <- "/mnt/data/lyx/IMC/analysis/abundance/"

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
print(colnames(colData(sce)))

clinicalFeatures <- c(
    "Tissue", "Prognosis", "RFS_status", "RFS_time", "Recurrence_site", "fong_score", "Age", "Gender", "KRAS_mutation", "BRAF_mutation", "mTBS",
    "CRLM_number", "CRLM_size", "Live_involvement_num", "CEA", "CA199", "Pathology", "Differential_grad", "T_grade", "Lymph_grade"
)

countdf <- Transform_CellCountMat(sce, group = c("IM", "CT", "TAT"), clinicalFeatures = clinicalFeatures, is.fraction = T)

## Binary clinical information
## Good Prognosis: 0, Worse: 1
## Non-Relapse: 0, Relapse: 1
countdf$Recurrence_site <- ifelse(countdf$Recurrence_site == "Liver", 0, 1) ## Recurrence in Liver: 0, Extrahepatic: 1
countdf$fong_score <- ifelse(as.numeric(countdf$fong_score) <= 2, 0, 1) ## International Liver-recurrence score, 0: low-risk, 1: high-risk
countdf$Age <- ifelse(as.numeric(countdf$Age) >= 60, 1, 0) ## Age: >=60: 1, otherwise: 0
countdf$Gender <- ifelse(as.numeric(countdf$Gender) == 1, 0, 1) ## Gender, 0: Male, 1: Female
## KRAS mutation: 1: Mut, 0: WT
## BRAF mutation: 1: Mut, 0: WT
countdf$mTBS <- ifelse(as.numeric(countdf$mTBS) >= median(as.numeric(countdf$mTBS)), 1, 0) ## Tumor Burden Score, 1: high, 0: low
countdf$CRLM_number <- ifelse(as.numeric(countdf$CRLM_number) >= median(as.numeric(countdf$CRLM_number)), 1, 0) ## CRC-Liver metastasis number, 1: high, 0: low
countdf$CRLM_size <- ifelse(as.numeric(countdf$CRLM_size) >= median(as.numeric(countdf$CRLM_size)), 1, 0) ## CRC-Liver metastasis Size, 1: high, 0: low
countdf$countdf$Live_involvement_num <- ifelse(as.numeric(countdf$countdf$Live_involvement_num) == 1, 0, 1) ## Live involvement number, 0: 1, 1: 2
countdf$CEA <- ifelse(as.numeric(countdf$CEA) >= 5, 1, 0) ## Pre-operation value, 1: abnormal, 0: normal
countdf$CA199 <- ifelse(as.numeric(countdf$CA199) >= 37, 1, 0) ## Pre-operation value, 1: abnormal, 0: normal
## Pathology: 0: Adenocarcinoma, 1: Mucinous adenocarcinoma
## Differential grade: 1: highly-diffrenciatial, 2: mediate, 3: low
## Tumor grade
## Lymhonode grade: 1: positive, 0: negative

xCol <- c(1, 21)
clinicalFeatures <- c(
    "Prognosis", "RFS_status", "Recurrence_site", "fong_score", "Age", "Gender", "KRAS_mutation", "BRAF_mutation", "mTBS",
    "CRLM_number", "CRLM_size", "Live_involvement_num", "CEA", "CA199", "Pathology", "Lymph_grade"
)

for (tissue in c("IM", "CT", "TAT")) {
    CountMAT <- abundanceDotPlotMat(countdf, xCol = xCol, clinicalFeatures = clinicalFeatures, tissue = tissue)
    MultiCliDotplot(CountMAT, tissue = tissue, savePath)
}


## Abundance Boxplot of relapse
clinicalFeatures <- c("Tissue","RFS_status")
countdf <- Transform_CellCountMat(sce, group = c("IM", "CT", "TAT"), clinicalFeatures = clinicalFeatures, is.fraction = T)

dim(countdf)
head(countdf)

xCol <- c(1, 21)
celltypes2Plot <- c("DPT", "B", "Mono_CLEC9A", "Mono_Classic", "Mono_CD57", "Macro_Multi", "Macro_CD11b", "SC_Vimentin", "TC_Ki67", "TC_TIGHT")
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
    stat_compare_means(aes(group = RFS_status), label.y = 0.5, method = "t.test")

pdf("/mnt/data/lyx/IMC/analysis/abundance/relapse_abundance_analysis.pdf", width = 15, height = 8)
print(p)
dev.off()

## Abundance Boxplot of KRAS mutation
clinicalFeatures <- c("Tissue","KRAS_mutation")
countdf <- Transform_CellCountMat(sce, group = c("IM", "CT", "TAT"), clinicalFeatures = clinicalFeatures, is.fraction = T)

dim(countdf)
head(countdf)

xCol <- c(1, 21)
celltypes2Plot <- c("DPT", "B", "Mono_CLEC9A", "Mono_Classic", "Mono_CD57", "SC_Vimentin", "TC_CAIX", "TC_TIGHT","TC_COLLAGEN")
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


## KM curve
clinicalFeatures <- c("Tissue","RFS_time","RFS_status")
celltypes2Plot <- c("DPT", "B", "Mono_CLEC9A", "Mono_Classic", "Mono_CD57", "Macro_Multi", "Macro_CD11b", "SC_Vimentin", "TC_Ki67", "TC_TIGHT")
countdf <- Transform_CellCountMat(sce, group = c("IM", "CT", "TAT"), clinicalFeatures = clinicalFeatures, is.fraction = T)
for (tissue in c("IM", "CT", "TAT")) {
    plotdf <- MergeCountsofROI(countdf, tissue = tissue, expCol = xCol, scale = F)

    for (feature in clinicalFeatures) {
    plotdf[, feature] <- clinical[match(rownames(plotdf), clinical$PID), feature]
  }

    for (celltype in celltypes2Plot) {
        KMVisualize(plotdf, celltype, cutoff = "best", savePath = paste0(savePath, "KM/", tissue, "/",tissue, "_", celltype, "_KMCurve.pdf"))
    }
}

## meta analysis
if(F){
    features <- c(
    "zs_rec_riskmodel", "fong_score", "Age", "Gender", "KRAS_mutation", "BRAF_mutation", "mTBS",
    "CRLM_number", "CRLM_size", "Live_involvement_num", "Pathology", "Differential_grad", "T_grade"
)
IM_plotdf <- MergeCountsofROI(countdf, tissue = "TAT", expCol = xCol, scale = T)
abundanceMetaAnalysis(IM_plotdf, celltypes2Plot, clinical, features, savePath)
}

## correlation analysis
for(tissut in c("IM", "CT")){
    plotdf <- subset(countdf, Tissue == tissut)
plotdf <- plotdf[, c(xCol[1]:xCol[2])]
CorrelationAnalysis(plotdf, savePath = paste0(savePath, tissut,"_abundance_correlation.pdf"))
}
