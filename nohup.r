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
library(parallel)
source("/home/lyx/project/IMC/clustering_functions.r")
source("/home/lyx/project/IMC/FlowSOM_metaClustering.r")

savePath <- "/mnt/data/lyx/IMC/analysis/clustering/"

total.res2.qc <- readRDS(paste0(savePath, "total.res2.qc.rds"))
total.res2.qc.norm <- readRDS(paste0(savePath, "total.res2.qc.norm.rds"))
total.res2.qc.norm <- readRDS(paste0(savePath, "total.res2.qc.norm.rds"))

panel <- read.csv("/mnt/data/lyx/IMC/IMC_CRC_panel.csv", stringsAsFactors = FALSE)
MarkerList <- load_Marker(panel)
markers <- MarkerList[["Iden_Marker"]]
print((markers))

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

    return(NULL)
}

phenoks <- 20
list_ <- lapply(phenoks, FUN = function(x) {
    return(list(total.res2.qc, total.res2.qc.norm, markers, savePath, x))
})
lapply(list_, perFunc)
