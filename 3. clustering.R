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
raw_csv_dir <- "/mnt/data/lyx/IMC/IMCell_Output/protein/"
meta_csv_dir <- "/mnt/data/lyx/IMC/IMC_CRC_QC.csv"

total.res2.qc <- load_RawData(raw_csv_dir, meta_csv_dir)
print(dim(total.res2.qc))
total.res2.qcBack <- total.res2.qc

## Load Marker
panel <- read.csv("/mnt/data/lyx/IMC/IMC_CRC_panel.csv", stringsAsFactors = FALSE)
MarkerList <- load_Marker(panel)

## data normalization
savePath <- "/mnt/data/lyx/IMC/analysis/clustering/"

colnames(total.res2.qc)[which(colnames(total.res2.qc) == "sample")] <- "filelist"

total.res2.qc.norm <- normData(total.res2.qc, marker_total = (MarkerList[["All_Marker"]]), censor_val = 0.99, arcsinh = FALSE, norm_method = "0-1")
total.res2.qc.norm.asin <- normData(total.res2.qc, marker_total = (MarkerList[["All_Marker"]]), censor_val = 0.99, arcsinh = TRUE, norm_method = "0-1")

saveRDS(total.res2.qc, paste0(savePath, "total.res2.qc.rds"))
saveRDS(total.res2.qc.norm, paste0(savePath, "total.res2.qc.norm.rds"))
saveRDS(total.res2.qc.norm.asin, paste0(savePath, "total.res2.qc.norm.asin.rds"))

total.res2.qc <- readRDS(paste0(savePath, "total.res2.qc.rds"))
total.res2.qc.norm <- readRDS(paste0(savePath, "total.res2.qc.norm.rds"))
total.res2.qc.norm.asin <- readRDS(paste0(savePath, "total.res2.qc.norm.asin.rds"))

# Flowsom clustering

## norm, Identification markers
markers <- MarkerList[["Iden_Marker"]]
print((markers))

### paralle
if (F) {
    library(parallel)
    perFunc <- function(list_) {
        total.res2.qc <- list_[[1]]
        total.res2.qc.norm <- list_[[2]]
        markers <- list_[[3]]
        savePath <- list_[[4]]
        phenok <- list_[[5]]

        cat("Performing phenograph clustering with phenok = ", phenok, "\n")
        ## norm
        norm_exp <- FlowSOM_clustering(total.res2.qc, total.res2.qc.norm, markers,
            phenographOnly = F, xdim = 100, ydim = 100, method = paste0("qc5_norm_Idenmarker"), savePath = savePath, phenok = phenok
        )
        saveRDS(norm_exp, paste0(savePath, "qc5_norm_Idenmarker_norm_exp_k", phenok, ".rds"))
        pdf(paste0(savePath, "qc5_norm_Idenmarker_norm_exp_k", phenok, ".pdf"), width = 15, height = 12)
        plotHeatmap(norm_exp, marker_total = markers, clustercol = "flowsom100pheno15")
        dev.off()

        ### norm.arcsin
        norm_exp <- FlowSOM_clustering(total.res2.qc, total.res2.qc.norm.asin, markers,
            phenographOnly = F, xdim = 100, ydim = 100, method = "qc5_norm_asin_Idenmarke", savePath = savePath, phenok = phenok
        )
        saveRDS(norm_exp, paste0(savePath, "qc5_norm_asin_Idenmarke_norm_exp_k", phenok, ".rds"))
        pdf(paste0(savePath, "qc5_norm_asin_Idenmarke_norm_exp_k", phenok, ".pdf"), width = 15, height = 12)
        plotHeatmap(norm_exp, marker_total = markers, clustercol = "flowsom100pheno15")
        dev.off()

        return(NULL)
    }
    phenoks <- c(10, 15, 20, 25, 30, 40)
    targets <- lapply(phenoks, FUN = function(x) {
        return(list(total.res2.qc, total.res2.qc.norm, markers, savePath, x))
    })
    cl <- makeCluster(8)
    results <- parLapply(cl, list_, perFunc)
    stopCluster(cl)
}

## only phenograph
### norm
FlowSOM_clustering(total.res2.qc, total.res2.qc.norm, phenographOnly = T, markers = markers, k = 100, method = "qc5_norm_marker20_pheno100", savePath)
### norm.arcsin
FlowSOM_clustering(total.res2.qc, total.res2.qc.norm.asin, phenographOnly = T, markers = markers, k = 100, method = "qc5_norm_asin_marker20_pheno100", savePath)


norm_exp <- readRDS(paste0(savePath, "qc5_norm_Idenmarker_norm_exp_k20.rds"))
pdf(paste0(savePath, "qc5_norm_Idenmarker_norm_exp_k20.pdf"), width = 15, height = 12)
plotHeatmap(norm_exp, marker_total = markers, clustercol = "flowsom100pheno15")
dev.off()

plotFeaturePlot("/mnt/data/lyx/IMC/analysis/clustering/qc5_norm_Idenmarker_flowsom_pheno_tsne_k25.csv", markers = markers, savePath = "/mnt/data/lyx/IMC/analysis/clustering/")

# Annotation
### Tumor: EPCAM
### Lymphcyte:
### Stromal: Collagen I
### Macrophage: CD68
### Monocyte: CD14

MajorTypeList <- c(
    "Monocyte", "Monocyte", "Monocyte", "Macrophage", "Tumor",
    "Tumor", "Tumor", "Monocyte", "Tumor", "Monocyte",
    "Lymphocyte", "Monocyte", "Tumor", "Monocyte", "Stromal",
    "Macrophage", "Stromal", "Tumor", "Macrophage", "Monocyte",
    "Macrophage", "Monocyte", "Stromal", "Tumor", "Lymphocyte",
    "Lymphocyte", "Lymphocyte"
)

## Lymphocyte
## T cell: CD3
## B cell: CD20
## NK: CD57

## Monocyte
## 1:  (Mono6 - Multi)
## 2: CD14+ CD16 CD11chi CLEC9A (Mono1 - CD11chi)
## 3: CD14+ CD16 CD11c CLEC9A (Mono2 - Classic Monocyte)
## 8: CD14+ CD16 CD11c+ CLEC9A+ (Mono3 - cDC1)
## 10: CD14+ CD16 CD11c CLEC9A CD169+ (Mono4 - SIGLEC1)
## 12: CD14+ CD16 CD11c CLEC9A (Mono2 - Classic Monocyte)
## 14: CD14+ CD16 CD11c CLEC9A (Mono2 - Classic Monocyte)
## 20: CD14+ CD16+ CD11c CLEC9A (Mono5 - Intermediate Monocyte)
## 22: CD14+ CD16 CD11c CLEC9A (Mono2 - Classic Monocyte)

## Macrophage
## 4: CD68+ CD14+ CD16+ CD11c CD11b CD163+ CLEC9A HLADR+ (Macro1 - Multi)
## 16: CD68+ CD14 CD16 CD11c CD11b CD163 CLEC9A HLADR+ (Macro2 - HLADR)
## 19: CD68+ CD14 CD16 CD11c CD11b+ CD163 CLEC9A HLADR (Macro3 - CD11b)
## 21: CD68+ CD14+ CD16 CD11c CD11b CD163 CLEC9A HLADR (Macro4 - CD14+)

## Mixture: CD14, CD11c, VEGF, CD4, FAP, FoxP3, CD57, CLEC9A, CD169hi

MinorTypeList <- c(
    "Mono_Multi", "Mono_CD11c", "Mono_Classic", "Macro_Multi", "TC_TIGHT",
    "TC_EPCAM", "TC_CAIX", "Mono_cDC_ITGAX", "TC_COLLAGEN", "Mono_SIGLEC1",
    "NK", "Mono_Classic", "TC_Ki67", "Mono_Classic", "SC_Vimentin",
    "Macro_HLADR", "SC_VEGF", "TC_HAVR2", "Macro_CD11b", "Mono_Intermediate",
    "Macro_CD14", "Mono_Classic", "SC_FAP", "TC_VEGF", "CD4T",
    "CD8T", "B"
)

# c(
#    "1", "2", "3", "4", "5",
#    "6", "7", "8", "9", "10",
#    "11", "12", "13", "14", "15",
#    "16", "17", "18", "19", "20",
#    "21", "22", "23", "24", "25",
#    "26", "27", "28", "29"
# )

Annotate(
    "/mnt/data/lyx/IMC/analysis/clustering/qc5_norm_Idenmarker_flowsom_pheno_tsne_k20.csv", "/mnt/data/lyx/IMC/analysis/clustering/qc5_norm_Idenmarker_flowsom_pheno_k20.csv",
    MajorTypeList, MinorTypeList, "/mnt/data/lyx/IMC/analysis/clustering/"
)

clustercsv <- readRDS("/mnt/data/lyx/IMC/analysis/clustering/annotate_allcells.rds")
write.csv(clustercsv, file = "/mnt/data/lyx/IMC/analysis/clustering/annotate_allcells.csv", row.names = FALSE) ## save annotate result

TsnecsvPath <- read.csv("/mnt/data/lyx/IMC/analysis/clustering/qc5_norm_Idenmarker_flowsom_pheno_tsne_k20.csv")
head(TsnecsvPath, 2)
table(TsnecsvPath$metacluster)

markersList <- MarkerList[["All_Marker"]]
markersList <- c(
    "CD45", "CD20", "CD3", "CD8a", "CD4", "TIGIT", "FoxP3", "CD57", "CD127", "CD366", ## Lymphocyte
    "HLADR", "CD68", "CD14", "CD11c", "CD11b", "CD16", "CLEC9A", "CD169", "CD163", ## Myeloid
    "CollagenI", "Vimentin", "AlphaSMA", "FAP", ## Stromal
    "EpCAM", "Ki67", "VEGF", "CAIX", "HK2", ## Tumor
    "FASN", "CD80", "CD274", "PRPS1", "CD279", "GLUT1", "CD27" ## Other
)

MarkerHeatmap(clustercsv, markersList)


## Density dotplot
## To perform this function, should run 4.abundance.r first
library(SingleCellExperiment)

savePath <- "/mnt/data/lyx/IMC/analysis/clustering/"

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
sce <- sce[,sce$Tissue == "IM"]
rownames(sce)

PlotDensityDotplot(sce, marker1 = "CD20", marker2 = "CD3", MajorType = "Lymphocyte", sampleSize = 50000, savePath)
PlotDensityDotplot(sce,marker1 = "EpCAM",marker2 = "CD274",MajorType = 'Tumor',sampleSize = 50000, savePath)
PlotDensityDotplot(sce,marker1 = "CD57",marker2 = "CD274",MajorType = 'Monocyte',sampleSize = 50000, savePath)
PlotDensityDotplot(sce,marker1 = "CD57",marker2 = "CD279",MajorType = 'Monocyte',sampleSize = 50000, savePath)
