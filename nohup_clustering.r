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

savePath <- "/mnt/data/lyx/IMC/analysis/clustering/"

## Load Marker
cat("Start Loading Marker", "\n")
panel <- read.csv("/mnt/data/lyx/IMC/IMC_CRC_panel_v2.csv", stringsAsFactors = FALSE)
MarkerList <- load_Marker(panel)



## Major
if (F) {
    ## Loading protein expression matrix and meta data
    raw_csv_dir <- "/mnt/data/lyx/IMC/IMCell_Output/denoise_result"
    meta_csv_dir <- "/mnt/data/lyx/IMC/IMC_CRC_QC.csv"

    cat("Start Loading Raw data", "\n")
    # total.res2.qc <- load_RawData(raw_csv_dir, meta_csv_dir)
    # print(dim(total.res2.qc))

    ## data normalization
    # colnames(total.res2.qc)[which(colnames(total.res2.qc) == "sample")] <- "filelist"

    cat("Start Normalizing Data", "\n")
    # total.res2.qc.norm <- normData(total.res2.qc, marker_total = (MarkerList[["All_Marker"]]), censor_val = 0.99, arcsinh = FALSE, is.Batch = F, norm_method = "0-1")
    # total.res2.qc.norm.asin <- normData(total.res2.qc, marker_total = (MarkerList[["All_Marker"]]), censor_val = 0.99, arcsinh = TRUE, is.Batch = F, norm_method = "0-1")

    # saveRDS(total.res2.qc, paste0(savePath, "total.res2.qc.rds"))
    # saveRDS(total.res2.qc.norm, paste0(savePath, "total.res2.qc.norm.rds"))
    # saveRDS(total.res2.qc.norm.asin, paste0(savePath, "total.res2.qc.norm.asin.rds"))

    cat("Start Loading Data", "\n")
    total.res2.qc <- readRDS(paste0(savePath, "total.res2.qc.rds"))
    total.res2.qc.norm <- readRDS(paste0(savePath, "total.res2.qc.norm.rds"))
    total.res2.qc.norm.asin <- readRDS(paste0(savePath, "total.res2.qc.norm.asin.rds"))
    cat("Loading Data done", "\n")

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
                phenographOnly = F, xdim = 100, ydim = 100, Type = Type, method = paste0("qc5_norm_Idenmarker"), savePath = savePath, phenok = phenok
            )
            # saveRDS(norm_exp, paste0(savePath, Type, "_qc5_norm_Idenmarker_norm_exp_k", phenok, ".rds"))
            pdf(paste0(savePath, Type, "_qc5_norm_Idenmarker_norm_exp_k", phenok, ".pdf"), width = 15, height = 12)
            plotHeatmap(norm_exp, marker_total = markers, clustercol = paste0(Type, "_flowsom100pheno15"))
            dev.off()

            return(NULL)
        }
        phenoks <- c(20, 15, 10)
        Type <- "Major"
        targets <- lapply(phenoks, FUN = function(x) {
            return(list(total.res2.qc, total.res2.qc.norm, total.res2.qc.norm.asin, markers, Type, savePath, x))
        })
        cl <- makeCluster(24)
        results <- parLapply(cl, targets, perFunc)
        stopCluster(cl)
    }
}

## Minor
if (F) {
    cat("Start Loading Data", "\n")
    total.res2.qc <- readRDS(paste0(savePath, "total.res2.qc.rds"))
    total.res2.qc.norm <- readRDS(paste0(savePath, "total.res2.qc.norm.rds"))
    total.res2.qc.norm.asin <- readRDS(paste0(savePath, "total.res2.qc.norm.asin.rds"))

    ## Minor clustering
    cat("Loading Major annotation", "\n")
    bestnorm_expPath <- "/mnt/data/lyx/IMC/analysis/clustering/Major_qc5_norm_Idenmarker_flowsom_pheno_k15.rds"
    bestnorm_exp <- readRDS(bestnorm_expPath)

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
            pdf(paste0(savePath, Type, "_qc5_norm_Idenmarker_norm_exp_k", phenok, ".pdf"), width = 15, height = 12)
            plotHeatmap(norm_exp, marker_total = markers, clustercol = paste0(Type, "_flowsom100pheno15"))
            dev.off()

            return(NULL)
        }

        ## form targets to multiple precess
        cat("Form the multi-process targets", "\n")
        phenoks <- c(25, 20, 15, 10)
        Type <- c("Lymphocyte", "Myeloid", "Stromal", "Tumor")

        targets <- list()
        z <- 1
        for (i in 1:length(Type)) {
            idx <- bestnorm_exp$MajorType %in% Type[i]
            for (j in phenoks) {
                cat("Precessing ", Type[i], " of k=", j, "\n")
                targets[[z]] <- list(bestnorm_exp[idx, ], bestnorm_exp[idx, ], NULL, MarkerList[[i + 2]], Type[i], savePath, j)
                z <- z + 1
            }
        }
        cat("Form the multi-process targets done.", "\n")
        cl <- makeCluster(24)
        cat("Start multi-process.", "\n")
        results <- parLapply(cl, targets, perFunc)
        stopCluster(cl)
    }
}

## Unknown
if (T) {
    ## Minor clustering
    cat("Loading Major annotation", "\n")
    bestnorm_expPath <- paste0(savePath, "annotate_allcells.rds")
    bestnorm_exp <- readRDS(bestnorm_expPath)

    bestnorm_exp[which(bestnorm_exp$SubType == "UNKNOWN"), ]$MajorType <- "UNKNOWN"
    markers <- unique(c(MarkerList[[2]], MarkerList[[3]], MarkerList[[4]], MarkerList[[5]], MarkerList[[6]]))

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
            pdf(paste0(savePath, Type, "_qc5_norm_Idenmarker_norm_exp_k", phenok, ".pdf"), width = 15, height = 12)
            plotHeatmap(norm_exp, marker_total = markers, clustercol = paste0(Type, "_flowsom100pheno15"))
            dev.off()

            return(NULL)
        }

        ## form targets to multiple precess
        cat("Form the multi-process targets", "\n")
        phenoks <- c(25, 20, 15, 10)
        Type <- c("UNKNOWN")

        targets <- list()
        z <- 1
        for (i in 1:length(Type)) {
            idx <- bestnorm_exp$MajorType %in% Type[i]
            for (j in phenoks) {
                cat("Precessing ", Type[i], " of k=", j, "\n")
                targets[[z]] <- list(bestnorm_exp[idx, ], bestnorm_exp[idx, ], NULL, markers, Type[i], savePath, j)
                z <- z + 1
            }
        }
        cat("Form the multi-process targets done.", "\n")
        cl <- makeCluster(24)
        cat("Start multi-process.", "\n")
        results <- parLapply(cl, targets, perFunc)
        stopCluster(cl)
    }
}
