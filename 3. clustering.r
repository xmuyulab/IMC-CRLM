# clustering
library(RColorBrewer)
library(pheatmap)
library(ComplexHeatmap)
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

## Loading protein expression matrix and meta data
raw_csv_dir <- "/mnt/data/lyx/IMC/IMCell_Output/denoise_result"
meta_csv_dir <- "/mnt/data/lyx/IMC/IMC_CRC_QC.csv"

total.res2.qc <- load_RawData(raw_csv_dir, meta_csv_dir)
print(dim(total.res2.qc))
total.res2.qcBack <- total.res2.qc

## Load Marker
panel <- read.csv("/mnt/data/lyx/IMC/IMC_CRC_panel_v2.csv", stringsAsFactors = FALSE)
MarkerList <- load_Marker(panel)

## data normalization
savePath <- "/mnt/data/lyx/IMC/analysis/clustering/"

colnames(total.res2.qc)[which(colnames(total.res2.qc) == "sample")] <- "filelist"

total.res2.qc.norm <- normData(total.res2.qc, marker_total = (MarkerList[["All_Marker"]]), censor_val = 0.99, arcsinh = FALSE, is.Batch = T, norm_method = "0-1")
total.res2.qc.norm.asin <- normData(total.res2.qc, marker_total = (MarkerList[["All_Marker"]]), censor_val = 0.99, arcsinh = TRUE, is.Batch = T, norm_method = "0-1")

saveRDS(total.res2.qc, paste0(savePath, "total.res2.qc.rds"))
saveRDS(total.res2.qc.norm, paste0(savePath, "total.res2.qc.norm.rds"))
saveRDS(total.res2.qc.norm.asin, paste0(savePath, "total.res2.qc.norm.asin.rds"))

total.res2.qc <- readRDS(paste0(savePath, "total.res2.qc.rds"))
total.res2.qc.norm <- readRDS(paste0(savePath, "total.res2.qc.norm.rds"))
total.res2.qc.norm.asin <- readRDS(paste0(savePath, "total.res2.qc.norm.asin.rds"))

# Flowsom clustering

## norm, Identification markers
markers <- MarkerList[["Major_Marker"]]
print((markers))

### paralle
if (F) {
    library(parallel)

    perFunc <- function(list_) {
        source("/home/lyx/project/IMC/clustering_functions.r")
        source("/home/lyx/project/IMC/FlowSOM_metaClustering.r")
        total.res2.qc <- list_[[1]]
        total.res2.qc.norm <- list_[[2]]
        markers <- list_[[3]]
        savePath <- list_[[4]]
        phenok <- list_[[5]]

        cat("Performing phenograph clustering with phenok = ", phenok, "\n")
        ## norm
        norm_exp <- FlowSOM_clustering(total.res2.qc, total.res2.qc.norm, markers,
            phenographOnly = F, xdim = 100, ydim = 100, method = paste0("qc5_norm_Idenmarker"), savePath = savePath, phenok = phenok, using.kmenas = T
        )
        saveRDS(norm_exp, paste0(savePath, "qc5_norm_Idenmarker_norm_exp_k", phenok, ".rds"))
        pdf(paste0(savePath, "qc5_norm_Idenmarker_norm_exp_k", phenok, ".pdf"), width = 15, height = 12)
        plotHeatmap(norm_exp, marker_total = markers, clustercol = "flowsom100pheno15")
        dev.off()

        ### norm.arcsin
        norm_exp <- FlowSOM_clustering(total.res2.qc, total.res2.qc.norm.asin, markers,
            phenographOnly = F, xdim = 100, ydim = 100, method = "qc5_norm_asin_Idenmarke", savePath = savePath, phenok = phenok, using.kmenas = T
        )
        saveRDS(norm_exp, paste0(savePath, "qc5_norm_asin_Idenmarke_norm_exp_k", phenok, ".rds"))
        pdf(paste0(savePath, "qc5_norm_asin_Idenmarke_norm_exp_k", phenok, ".pdf"), width = 15, height = 12)
        plotHeatmap(norm_exp, marker_total = markers, clustercol = "flowsom100pheno15")
        dev.off()

        return(NULL)
    }
    phenoks <- c(20, 10, 5, 1, 0.5, 0.1, 0.01)
    targets <- lapply(phenoks, FUN = function(x) {
        return(list(total.res2.qc, total.res2.qc.norm, markers, savePath, x))
    })
    cl <- makeCluster(7)
    results <- parLapply(cl, targets, perFunc)
    stopCluster(cl)
}

## only phenograph
if (F) {
    ### norm
    FlowSOM_clustering(total.res2.qc, total.res2.qc.norm, phenographOnly = T, markers = markers, k = 100, method = "qc5_norm_marker20_pheno100", savePath)
    ### norm.arcsin
    FlowSOM_clustering(total.res2.qc, total.res2.qc.norm.asin, phenographOnly = T, markers = markers, k = 100, method = "qc5_norm_asin_marker20_pheno100", savePath)
    norm_exp <- readRDS(paste0(savePath, "qc5_norm_Idenmarker_norm_exp_k20.rds"))
    pdf(paste0(savePath, "qc5_norm_Idenmarker_norm_exp_k20.pdf"), width = 15, height = 12)
    plotHeatmap(norm_exp, marker_total = markers, clustercol = "flowsom100pheno15")
    dev.off()
}


# plotFeaturePlot("/mnt/data/lyx/IMC/analysis/clustering/qc5_norm_Idenmarker_flowsom_pheno_tsne_k25.csv", markers = markers, savePath = "/mnt/data/lyx/IMC/analysis/clustering/")

# Annotation
### Tumor: EPCAM
### Lymphcyte:
### Stromal: Collagen I
### Macrophage: CD68
### Monocyte: CD14
### Unknown

MajorTypeList <- c(
    "Lymphocyte", "Lymphocyte", "Stromal", "Tumor", "Tumor",
    "Tumor", "Macrophage", "Monocyte", "Monocyte", "Stromal",
    "Unknown", "Tumor", "Stromal", "Tumor", "Macrophage",
    "Monocyte", "Stromal", "Macrophage", "Monocyte", "Monocyte",
    "Macrophage", "Lymphocyte", "Lymphocyte", "Lymphocyte"
)

## Lymphocyte
## 17: B cell
## 8: DPT
## 1：NK

## Monocyte
## 3: CD14+ CD16+ (Mono_Intermediate)
## 4: CD14+ CLEC9A (Mono_CLEC9A)
## 5: CD14+ SIGLEC1 (Mono_SIGLEC1)
## 6: CD14+ (Mono_Classic)
## 11: CD14+ CD11c (Mono_CD11c)
## 14: CD14+ CD11c (Mono_CD11c)

## Macrophage
## 2: CD68+ CD14+ CD16+ CD163+ HLADR+ (Macro_Multi)
## 15: CD68- CD14 CD16 CD11c CD11b CD163 CLEC9A HLADR+ (Macro_HLADR)
## 21: CD68+ CD11b (Macro_CD11b)


MinorTypeList <- c(
    "Treg", "NK", "SC_VEGF", "TC_CAIX", "TC_VEGF",
    "TC_Ki67", "Macro_Multi", "Mono_CD11c", "cDC1", "SC_aSMA",
    "Unknown", "TC_EpCAM", "SC_FAP", "TC_COLLAGEN", "Macro_CD11b",
    "Mono_CD169", "SC_COLLAGEN", "Macro_Vimentin", "Mono_Intermediate", "Mono_Classic",
    "Macro_CD14", "CD4T", "CD8T", "B"
)

# c(
#    "1", "2", "3", "4", "5",
#    "6", "7", "8", "9", "10",
#    "11", "12", "13", "14", "15",
#    "16", "17", "18", "19", "20",
#    "21", "22", "23", "24", "25",
#    "26", "27", "28", "29","30"
# )

k <- 20
Annotate(
    paste0(savePath, "qc5_norm_Idenmarker_flowsom_pheno_tsne_k", k, ".csv"), paste0(savePath, "qc5_norm_Idenmarker_flowsom_pheno_k", k, ".csv"),
    MajorTypeList, MinorTypeList, "/mnt/data/lyx/IMC/analysis/clustering/"
)

clustercsv <- readRDS(paste0(savePath, "annotate_allcells.rds"))
write.csv(clustercsv, file = paste0(savePath, "annotate_allcells.csv"), row.names = FALSE) ## save annotate result

TsnecsvPath <- read.csv(paste0(savePath, "qc5_norm_Idenmarker_flowsom_pheno_tsne_k", k, ".csv"))
head(TsnecsvPath, 2)
table(TsnecsvPath$metacluster)

markersList <- MarkerList[["All_Marker"]]
markersList <- c(
    "CD45", "CD20", "CD3", "CD8a", "CD4", "FoxP3", "CD57", ## Lymphocyte
    "HLADR", "CD68", "CD14", "CD11c", "CD11b", "CD16", "CLEC9A", "CD169", "CD163", ## Myeloid
    "CollagenI", "Vimentin", "AlphaSMA", "FAP", ## Stromal
    "EpCAM", ## Tumor
    "Ki67", "HK2", "FASN", "PRPS1", "GLUT1", ## Cell grouth
    "VEGF", "CAIX", ## Hypoxia
    "CD127", "CD27", "CD80", ## Immuno-activation
    "CD274", "CD279", "TIGIT", "CD366" ## Immuno-checkpoint
)

MarkerHeatmap(clustercsv, markersList, savePath)

## Clustering assement
clustercsv <- readRDS(paste0(savePath, "annotate_allcells.rds"))
markers <- MarkerList[["Iden_Marker"]]
clusterCol <- "SubType"

AssessClustering(clustercsv, markers = markers, clusterCol = clusterCol, tFraction = 0.85, sampleSize = 50000, savePath = savePath)


## Two-run clustering
savePath <- "/mnt/data/lyx/IMC/analysis/clustering/"

total.res2.qc <- readRDS(paste0(savePath, "total.res2.qc.rds"))
total.res2.qc.norm <- readRDS(paste0(savePath, "total.res2.qc.norm.rds"))
total.res2.qc.norm.asin <- readRDS(paste0(savePath, "total.res2.qc.norm.asin.rds"))

dim(total.res2.qc.norm)
colnames(total.res2.qc.norm)

## Major clustering
## Load Marker
panel <- read.csv("/mnt/data/lyx/IMC/IMC_CRC_panel_v2.csv", stringsAsFactors = FALSE)
MarkerList <- load_Marker(panel)

markers <- MarkerList[["Major_Marker"]]
print((markers))

### paralle
if (T) {
    library(parallel)

    perFunc <- function(list_) {
        source("/home/lyx/project/IMC/clustering_functions.r")
        source("/home/lyx/project/IMC/FlowSOM_metaClustering.r")
        total.res2.qc <- list_[[1]]
        total.res2.qc.norm <- list_[[2]]
        total.res2.qc.norm.asin <- list_[[3]]
        markers <- list_[[4]]
        Type <- list_[[5]]
        savePath <- list_[[6]]
        phenok <- list_[[7]]

        cat("Performing phenograph clustering with phenok = ", phenok, "\n")
        ## norm
        norm_exp <- FlowSOM_clustering(total.res2.qc, total.res2.qc.norm, markers,
            phenographOnly = F, xdim = 100, ydim = 100, Type = Type, method = paste0("qc5_norm_Idenmarker"), savePath = savePath, phenok = phenok,
        )
        # saveRDS(norm_exp, paste0(savePath, Type, "_qc5_norm_Idenmarker_norm_exp_k", phenok, ".rds"))
        pdf(paste0(savePath, Type, "_qc5_norm_Idenmarker_norm_exp_k", phenok, ".pdf"), width = 15, height = 12)
        plotHeatmap(norm_exp, marker_total = markers, clustercol = paste0(Type, "_flowsom100pheno15"))
        dev.off()

        ### norm.arcsin
        norm_exp <- FlowSOM_clustering(total.res2.qc, total.res2.qc.norm.asin, markers,
            phenographOnly = F, xdim = 100, ydim = 100, Type = Type, method = "qc5_norm_asin_Idenmarke", savePath = savePath, phenok = phenok
        )
        # saveRDS(norm_exp, paste0(savePath, Type, "_qc5_norm_asin_Idenmarke_norm_exp_k", phenok, ".rds"))
        pdf(paste0(savePath, Type, "_qc5_norm_asin_Idenmarke_norm_exp_k", phenok, ".pdf"), width = 15, height = 12)
        plotHeatmap(norm_exp, marker_total = markers, clustercol = paste0(Type, "_flowsom100pheno15"))
        dev.off()

        return(NULL)
    }
    phenoks <- c(30, 20, 15, 10, 5)
    Type <- "Major"
    targets <- lapply(phenoks, FUN = function(x) {
        return(list(total.res2.qc, total.res2.qc.norm, total.res2.qc.norm.asin, markers, Type, savePath, x))
    })
    cl <- makeCluster(32)
    results <- parLapply(cl, targets, perFunc)
    stopCluster(cl)
}

## Major annotation
bestnorm_expPath <- "/mnt/data/lyx/IMC/analysis/clustering/Major_qc5_norm_Idenmarker_flowsom_pheno_k10.csv"
bestnorm_exp <- read.csv(bestnorm_expPath)
table(bestnorm_exp$Major_flowsom100pheno15)

MajorTypeList <- c(
    "Tumor", "UNKNOWN", "Tumor", "Lymphocyte", "UNKNOWN",
    "UNKNOWN", "UNKNOWN", "Myeloid", "Tumor", "Myeloid",
    "Myeloid", "UNKNOWN", "Myeloid", "UNKNOWN", "Tumor",
    "UNKNOWN", "Myeloid", "Myeloid", "UNKNOWN", "Stromal",
    "UNKNOWN", "Tumor", "UNKNOWN", "Stromal", "Lymphocyte",
    "Myeloid", "Stromal", "UNKNOWN", "Lymphocyte", "Myeloid"
)
bestnorm_exp$MajorType <- NA
clusterids <- as.numeric(names(table(bestnorm_exp$Major_flowsom100pheno15)))

for (clusterid in clusterids) {
    bestnorm_exp[bestnorm_exp$Major_flowsom100pheno15 %in% clusterid, ]$MajorType <- MajorTypeList[clusterid]
}
table(bestnorm_exp$MajorType)
bestnorm_exp$CellID <- as.character(1:nrow(bestnorm_exp))
saveRDS(bestnorm_exp, "/mnt/data/lyx/IMC/analysis/clustering/Major_qc5_norm_Idenmarker_flowsom_pheno_k15.rds")

MajorTypes <- names(table(bestnorm_exp$MajorType))

## Minor clustering
cat("Loading Major annotation", "\n")
bestnorm_expPath <- "/mnt/data/lyx/IMC/analysis/clustering/Major_qc5_norm_Idenmarker_flowsom_pheno_k15.rds"
bestnorm_exp <- readRDS(bestnorm_expPath)
bestnorm_exp$SubType <- "UNKNOWN"

### paralle
if (T) {
    library(parallel)

    perFunc <- function(list_) {
        source("/home/lyx/project/IMC/clustering_functions.r")
        source("/home/lyx/project/IMC/FlowSOM_metaClustering.r")
        total.res2.qc <- list_[[1]]
        total.res2.qc.norm <- list_[[2]]
        total.res2.qc.norm.asin <- list_[[3]]
        markers <- list_[[4]]
        Type <- list_[[5]]
        savePath <- list_[[6]]
        phenok <- list_[[7]]

        cat("Performing phenograph clustering with phenok = ", phenok, "\n")
        ## norm
        norm_exp <- FlowSOM_clustering(total.res2.qc, total.res2.qc.norm, markers,
            phenographOnly = F, xdim = 100, ydim = 100, Type = Type, method = paste0("qc5_norm_Idenmarker"), savePath = savePath, phenok = phenok
        )
        # saveRDS(norm_exp, paste0(savePath, Type, "_qc5_norm_Idenmarker_norm_exp_k", phenok, ".rds"))
        pdf(paste0(savePath, Type, "_qc5_norm_Idenmarker_norm_exp_k", phenok, ".pdf"), width = 15, height = 12)
        plotHeatmap(norm_exp, marker_total = markers, clustercol = paste0(Type, "_flowsom100pheno15"))
        dev.off()

        ### norm.arcsin
        norm_exp <- FlowSOM_clustering(total.res2.qc, total.res2.qc.norm.asin, markers,
            phenographOnly = F, xdim = 100, ydim = 100, Type = Type, method = "qc5_norm_asin_Idenmarke", savePath = savePath, phenok = phenok
        )
        # saveRDS(norm_exp, paste0(savePath, Type, "_qc5_norm_asin_Idenmarke_norm_exp_k", phenok, ".rds"))
        pdf(paste0(savePath, Type, "_qc5_norm_asin_Idenmarke_norm_exp_k", phenok, ".pdf"), width = 15, height = 12)
        plotHeatmap(norm_exp, marker_total = markers, clustercol = paste0(Type, "_flowsom100pheno15"))
        dev.off()

        return(NULL)
    }

    ## form targets to multiple precess
    cat("Form the multi-process targets", "\n")
    phenoks <- c(20, 15, 10, 5)
    Type <- c("Lymphocyte", "Myeloid", "Stromal", "Tumor", "Unknown")

    targets <- list()
    z <- 1
    for (i in 1:length(Type)) {
        idx <- bestnorm_exp$MajorType %in% Type[i]
        for (j in phenoks) {
            cat("Precessing ", Type[i], " of k=", j, "\n")
            targets[[z]] <- list(bestnorm_exp[idx, ], bestnorm_exp[idx, ], total.res2.qc.norm.asin[idx, ], MarkerList[[i + 2]], Type[i], savePath, j)()
            z <- z + 1
        }
    }
    cat("Form the multi-process targets done.", "\n")
    cl <- makeCluster(32)
    cat("Start multi-process.", "\n")
    results <- parLapply(cl, targets, perFunc)
    stopCluster(cl)
}

############### Lymphocyte
Lymnorm_expPath <- "/mnt/data/lyx/IMC/analysis/clustering/Lymphocyte_qc5_norm_Idenmarker_flowsom_pheno_k15.csv"
Lymnorm_exp <- read.csv(Lymnorm_expPath)

dim(Lymnorm_exp)
colnames(Lymnorm_exp)
Lymnorm_exp$SubType <- NA

LymphTypeList <- c(
    "B", "B", "UNKNOWN", "UNKNOWN", "B",
    "CD4T", "UNKNOWN", "UNKNOWN", "UNKNOWN", "NK",
    "UNKNOWN", "CD8T", "CD4T", "CD4T", "CD4T",
    "CD4T", "UNKNOWN", "CD4T", "CD4T", "NK",
    "NK", "CD8T", "Treg", "NK", "CD4T",
    "Treg", "CD8T", "CD4T", "CD8T", "Treg"
)

clusterids <- as.numeric(names(table(Lymnorm_exp$Lymphocyte_flowsom100pheno15)))
for (clusterid in clusterids) {
    Lymnorm_exp[Lymnorm_exp$Lymphocyte_flowsom100pheno15 %in% clusterid, ]$SubType <- LymphTypeList[clusterid]
}
table(Lymnorm_exp$SubType)
colnames(Lymnorm_exp)
bestnorm_exp[match(Lymnorm_exp$CellID, bestnorm_exp$CellID), ]$SubType <- Lymnorm_exp$SubType

############ Myeloid
Myeloid_expPath <- "/mnt/data/lyx/IMC/analysis/clustering/Myeloid_qc5_norm_Idenmarker_flowsom_pheno_k15.csv"
Myeloid_exp <- read.csv(Myeloid_expPath)

dim(Myeloid_exp)
colnames(Myeloid_exp)
Myeloid_exp$SubType <- NA

MyeloidTypeList <- c(
    "Macro_CD169", "Macro_CD169", "Macro_CD163", "Macro_CD11b", "Macro_HLADR",
    "Mono_CD11c", "Macro_HLADR", "UNKNOWN", "UNKNOWN", "Mono_Classic",
    "UNKNOWN", "Mono_Classic", "UNKNOWN", "Mono_Intermediate", "Mono_Classic",
    "UNKNOWN", "Mono_Classic", "Mono_Intermediate", "Mono_Classic", "Mono_CD11c",
    "Mono_Intermediate", "Macro_CD163", "Mono_CD11c", "Macro_CD163", "Mono_CD11c",
    "Macro_CD163", "Macro_CD163", "Macro_CD163", "Macro_CD163", "Mono_Classic"
)
clusterids <- as.numeric(names(table(Myeloid_exp$Myeloid_flowsom100pheno15)))
for (clusterid in clusterids) {
    Myeloid_exp[Myeloid_exp$Myeloid_flowsom100pheno15 %in% clusterid, ]$SubType <- MyeloidTypeList[clusterid]
}
table(Myeloid_exp$SubType)
colnames(Myeloid_exp)
bestnorm_exp[match(Myeloid_exp$CellID, bestnorm_exp$CellID), ]$SubType <- Myeloid_exp$SubType
table(bestnorm_exp$SubType)

############ Stromal
Stromal_expPath <- "/mnt/data/lyx/IMC/analysis/clustering/Stromal_qc5_norm_Idenmarker_flowsom_pheno_k20.csv"
Stromal_exp <- read.csv(Stromal_expPath)

dim(Stromal_exp)
colnames(Stromal_exp)
Stromal_exp$SubType <- NA

StromalTypeList <- c(
    "SC_aSMA", "SC_aSMA", "SC_FAP", "SC_COLLAGEN", "SC_COLLAGEN",
    "SC_COLLAGEN", "SC_aSMA", "UNKNOWN", "SC_COLLAGEN", "SC_Vimentin",
    "UNKNOWN", "SC_COLLAGEN", "UNKNOWN", "UNKNOWN", "UNKNOWN",
    "SC_COLLAGEN", "SC_COLLAGEN", "SC_Vimentin", "SC_COLLAGEN", "UNKNOWN",
    "SC_Vimentin", "UNKNOWN", "UNKNOWN", "UNKNOWN", "UNKNOWN",
    "UNKNOWN", "UNKNOWN"
)
clusterids <- as.numeric(names(table(Stromal_exp$Stromal_flowsom100pheno15)))
for (clusterid in clusterids) {
    Stromal_exp[Stromal_exp$Stromal_flowsom100pheno15 %in% clusterid, ]$SubType <- StromalTypeList[clusterid]
}
table(Stromal_exp$SubType)
colnames(Stromal_exp)
bestnorm_exp[match(Stromal_exp$CellID, bestnorm_exp$CellID), ]$SubType <- Stromal_exp$SubType
table(bestnorm_exp$SubType)

############ Tumor
Tumor_expPath <- "/mnt/data/lyx/IMC/analysis/clustering/Tumor_qc5_norm_Idenmarker_flowsom_pheno_k10.csv"
Tumor_exp <- read.csv(Tumor_expPath)

dim(Tumor_exp)
colnames(Tumor_exp)
Tumor_exp$SubType <- NA

TumorTypeList <- c(
    "UNKNOWN", "UNKNOWN", "TC_CAIX", "TC_CAIX", "TC_CAIX",
    "TC_CAIX", "TC_CAIX", "TC_EpCAM", "TC_CAIX", "TC_VEGF",
    "TC_EpCAM", "UNKNOWN", "UNKNOWN", "TC_VEGF", "TC_EpCAM",
    "TC_EpCAM", "TC_EpCAM", "TC_Ki67", "TC_EpCAM", "TC_VEGF",
    "TC_EpCAM", "TC_EpCAM", "TC_EpCAM", "TC_Ki67", "TC_VEGF",
    "TC_EpCAM", "TC_EpCAM", "TC_Ki67", "TC_VEGF", "TC_Ki67",
    "TC_Ki67", "UNKNOWN", "UNKNOWN", "TC_Ki67"
)

clusterids <- as.numeric(names(table(Tumor_exp$Tumor_flowsom100pheno15)))
for (clusterid in clusterids) {
    Tumor_exp[Tumor_exp$Tumor_flowsom100pheno15 %in% clusterid, ]$SubType <- TumorTypeList[clusterid]
}
table(Tumor_exp$SubType)
colnames(Tumor_exp)
bestnorm_exp[match(Tumor_exp$CellID, bestnorm_exp$CellID), ]$SubType <- Tumor_exp$SubType
table(bestnorm_exp$SubType)
table(bestnorm_exp$MajorType)

bestnorm_exp[which(bestnorm_exp$SubType == "UNKNOWN"), ]$MajorType <- "UNKNOWN"
saveRDS(bestnorm_exp, paste0(savePath, "annotate_allcells.rds"))

### Further clustering Unknown cells
bestnorm_exp <- readRDS(paste0(savePath, "annotate_allcells.rds"))
colnames(bestnorm_exp)
table(bestnorm_exp$MajorType)
table(bestnorm_exp$SubType)

## clustering Unknown cells
markers <- unique(c(MarkerList[[2]], MarkerList[[3]], MarkerList[[4]], MarkerList[[5]], MarkerList[[6]]))
Type <- "UNKNOWN"

Unknown_exp <- bestnorm_exp[which(bestnorm_exp$MajorType == "UNKNOWN"), ]
for (phenok in c(10, 20, 30, 40)) {
    norm_exp <- FlowSOM_clustering(Unknown_exp, Unknown_exp[, markers], markers,
        phenographOnly = F, xdim = 100, ydim = 100, Type = Type, method = paste0("qc5_norm_Idenmarker"), savePath = savePath, phenok = phenok,
    )
    p <- plotHeatmap(norm_exp, marker_total = markers, clustercol = paste0(Type, "_flowsom100pheno15"))
    pdf(paste0(savePath, Type, "_qc5_norm_Idenmarker_norm_exp_k", phenok, ".pdf"), width = 15, height = 12)
    print(p)
    dev.off()
}

# saveRDS(norm_exp, paste0(savePath, Type, "_qc5_norm_Idenmarker_norm_exp_k", phenok, ".rds"))

Unknown_exp <- read.csv("/mnt/data/lyx/IMC/analysis/clustering/UNKNOWN_qc5_norm_Idenmarker_flowsom_pheno_k15.csv")
Unknown_exp$SubType <- "UNKNOWN"
UnknownMajor <- c(
    "Myeloid", "Tumor", "Tumor", "UNKNOWN", "UNKNOWN",
    "Myeloid", "Myeloid", "UNKNOWN", "Stromal", "Myeloid",
    "Stromal", "Myeloid", "Stromal", "Myeloid", "Lymphocyte",
    "Myeloid", "Tumor", "Tumor", "Myeloid", "Lymphocyte",
    "Stromal", "Myeloid", "Myeloid", "Lymphocyte", "UNKNOWN",
    "Myeloid", "Tumor", "Tumor", "Lymphocyte", "Lymphocyte"
)
UnknownMinor <- c(
    "Macro_CD11b", "TC_VEGF", "TC_EpCAM", "UNKNOWN", "UNKNOWN",
    "Mono_Classic", "Macro_CD163", "UNKNOWN", "SC_Vimentin", "Mono_CD11c",
    "SC_FAP", "Mono_Classic", "SC_COLLAGEN", "Mono_CD11c", "NK",
    "Mono_Classic", "TC_CAIX", "TC_EpCAM", "Macro_CD169", "Treg",
    "SC_aSMA", "Mono_Intermediate", "Mono_CD11c", "Treg", "UNKNOWN",
    "Macro_HLADR", "TC_Ki67", "TC_EpCAM", "CD8T", "B"
)
clusterids <- as.numeric(names(table(Unknown_exp$UNKNOWN_flowsom100pheno15)))
for (clusterid in clusterids) {
    Unknown_exp[Unknown_exp$UNKNOWN_flowsom100pheno15 %in% clusterid, ]$MajorType <- UnknownMajor[clusterid]
    Unknown_exp[Unknown_exp$UNKNOWN_flowsom100pheno15 %in% clusterid, ]$SubType <- UnknownMinor[clusterid]
}
table(Unknown_exp$MajorType)
table(Unknown_exp$SubType)

bestnorm_exp[match(Unknown_exp$CellID, bestnorm_exp$CellID), ]$MajorType <- Unknown_exp$MajorType
bestnorm_exp[match(Unknown_exp$CellID, bestnorm_exp$CellID), ]$SubType <- Unknown_exp$SubType
table(bestnorm_exp$MajorType)
table(bestnorm_exp$SubType)
bestnorm_exp[which(bestnorm_exp$SubType == "UNKNOWN"), ]$MajorType <- "UNKNOWN"

## slightly chage

markersList <- c(
    "CD45", "CD20", "CD3", "CD8a", "CD4", "FoxP3", "CD57", ## Lymphocyte
    "HLADR", "CD68", "CD14", "CD11c", "CD11b", "CD16", "CLEC9A", "CD169", "CD163", ## Myeloid
    "CollagenI", "Vimentin", "AlphaSMA", "FAP", ## Stromal
    "EpCAM", "Ki67", "HK2", "VEGF", "CAIX" ## Tumor
    # "FASN", "PRPS1", "GLUT1", ## Cell grouth
    # "VEGF", "CAIX", ## Hypoxia
    # "CD127", "CD27", "CD80", ## Immuno-activation
    # "CD274", "CD279", "TIGIT", "CD366" ## Immuno-checkpoint
)

MarkerHeatmap(bestnorm_exp, markersList, savePath)
saveRDS(bestnorm_exp, paste0(savePath, "annotate_allcells.rds"))
write.csv(bestnorm_exp, file = paste0(savePath, "annotate_allcells.csv"), row.names = FALSE)
