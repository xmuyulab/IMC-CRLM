# Identify tumor microenviroment pattern
library(SingleCellExperiment)
library(pheatmap)
library(survival)
library(survminer)
source("./spatial_analysis_functions.r")
source("./structural_analysis_functions.r")

## Cellular neighbors analysis
sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
savePath <- "/mnt/data/lyx/IMC/analysis/structure/"
if (!dir.exists(savePath)) {
    dir.create(savePath, recursive = T)
}

## remain Tissue associated Tumor
sce <- sce[, sce$Tissue == "IM"]

scimapResult <- read.csv("/mnt/data/lyx/IMC/analysis/spatial/cellular_neighbor_IM.csv")
scimapResult <- scimapResult[, -1]

colnames(scimapResult)

colName <- colnames(scimapResult)[14:ncol(scimapResult)]
sce <- BindResult(sce, scimapResult, colName)

colnames(colData(sce))
table(colData(sce)$MajorType)

if (T) {
    for (structure in colName) {
        ## savepath
        savePathTemp1 <- paste0(savePath, structure, "/")
        if (!dir.exists(savePathTemp1)) {
            dir.create(savePathTemp1, recursive = T)
        }

        ## Cell subtype fraction in cellular neighbors pattern
        HeatmapForCelltypeInNeighbor(sce, "SubType", structure, savePathTemp1)

        ## Cellular pattern difference in Relaps and Non-Relaps
        CompareCellularPattern(sce, sep = "RFS_status", countcol = structure, n_cluster = 10, savePath = savePathTemp1)

        ## Celllular neighborhood pattern survival analysis
        CNP_countsDF <- GetAbundance(sce, countcol = structure, is.fraction = TRUE, is.reuturnMeans = T)
        CNPs <- names(table(colData(sce)[, structure]))
        if (!dir.exists(paste0(savePathTemp1, "KM/"))) {
            dir.create(paste0(savePathTemp1, "KM/"), recursive = T)
        }
        for (CNP in CNPs) {
            plotdf <- CNP_countsDF[, c(CNP, "RFS_status", "RFS_time")]
            KMForCNP(plotdf, CNP, savePath = paste0(savePathTemp1, "KM/", "Cellular Neighborhood pattern Suvival analysis of ", CNP, ".pdf"))
        }

        ## Cellular pattenr fraction
        CNP_countsDF <- GetAbundance(sce, countcol = structure, is.fraction = T, is.reuturnMeans = T)
        CNPFraction(CNP_countsDF, groupBy = "RFS_status", xCol = c(1, 10), savePath = savePathTemp1)

        ## CNP in image
        if (F) {
            SavePath1 <- paste0(savePathTemp1, "CNP_oncells/")
            if (!dir.exists(SavePath1)) {
                dir.create(SavePath1, recursive = T)
            }
            colData(sce)[, structure] <- as.factor(colData(sce)[, structure])

            ROIs <- names(table(colData(sce)$ID))
            for (ROI in ROIs) {
                PlotCelltypes(sce, ROI, TypeCol = structure, SavePath = paste0(SavePath1, ROI, "_"))
            }
        }

        ## Recluster
        ### clutering via metabolize molecular
        rownames(sce)
        reclusMarkers <- c(
            # "HLADR", "CD68", "CD14", "CD11c", "CD11b", "CD16", "CLEC9A", "CD169", "CD163", ## Cellidentity
            "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
            "VEGF", "CAIX", ## Hypoxia
            "CD279", "CD274", "CD127", ## Immune-checkpoint
            "CD27", "CD80" ## Immune-activation
        )

        ReMajorType <- c("Myeloid")
        ReclusterName <- "Myeloid"
        ReSubType <- NULL

        savePathTemp2 <- paste0(savePathTemp1, "reclustering/")
        if (!dir.exists(savePathTemp2)) {
            dir.create(savePathTemp2, recursive = T)
        }
        sce_ <- Reclustering(sce, reclusMarkers, ReMajorType, ReclusterName, ReSubType = ReSubType, PatternCol = structure, ncluster = 15, savePath = savePathTemp2)

        saveRDS(sce_, paste0(savePathTemp1, ReclusterName, "_recluster.rds"))
    }
}

## Certain reclustering types in cellular pattern (three group)
if (F) {
    ## Certain reclustering types in cellular pattern
    reclusMarkers <- c(
        # "HLADR", "CD68", "CD14", "CD11c", "CD11b", "CD16", "CLEC9A", "CD169", "CD163", ## Cellidentity
        "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
        "VEGF", "CAIX", ## Hypoxia
        "CD279", "CD274", "CD127", ## Immune-checkpoint
        "CD27", "CD80" ## Immune-activation
    )

    ReMajorType <- c("Myeloid")
    ReclusterName <- "Myeloid"

    for (structure in colName) {
        sce_ <- readRDS(paste0("/mnt/data/lyx/IMC/analysis/structure/", structure, "/Myeloid_recluster.rds"))
        interstType <- c("4")
        PlotCertainTypeinPattern(sce_, Col1 = ReclusterName, types1 = interstType, Col2 = structure, groupCol = "RFS_status", savePath = paste0(savePath, structure, "/reclustering/"))

        ##
        if (F) {
            interstType <- c("8")
            sce__ <- sce_[, sce_$Lymphocyte == interstType]
            table(sce__$MajorType)
            table(sce__$SubType)
            table(sce__$ID)[order(as.numeric(table(sce__$ID)), decreasing = T)]
            phenoLabelCountMat <- GetAbundance(sce_, countcol = "phenoLabel", is.fraction = FALSE, is.reuturnMeans = FALSE)
            phenoLabelCountMat <- TransformIntoPlotMat(phenoLabelCountMat, valueCol = c(1:2))
            head(phenoLabelCountMat)
            BoxPlotForPhenoAssCell(phenoLabelCountMat, savePath = paste0(savePath, structure, "/reclustering/"))
            SurvivalForPhenoAssCell(phenoLabelCountMat, savePath = paste0(savePath, structure, "/reclustering/"))
        }

        ## plot the differential expression genes
        sce__ <- sce_[, sce_$MajorType %in% ReclusterName]
        mat <- as.data.frame(t(assay(sce_)[reclusMarkers, ]))
        mat$phenoLabel <- 0
        mat$phenoLabel[colData(sce_)[, ReclusterName] %in% interstType] <- 1
        FCDF <- FCandPvalueCal(mat, xCol = c(1, length(reclusMarkers)), yCol = (length(reclusMarkers) + 1), need.sample = TRUE)
        VolcanoPlot(FCDF, pthreshold = 0.01, fcthreshold = 3, feature = "Phenotype-Associated cells", filename = paste0(savePath, structure, "/Phenotype-associated differential markers in ", ReMajorType, ".pdf"))


        ## Assign the new label
        ### CD274
        CD274_Pos <- interstType
        TempList <- AssignNewlabel(sce_, allROIs = names(table(sce$ID)), phenoLabel = "CD274_Type", ReclusterName = ReclusterName, interstType = CD274_Pos, cutoffType = "manual", cutoffValue = 6, numGroup = 3)
        sce_Temp <- TempList[[1]]
        CD274high_ROI <- TempList[[2]]
        CD274low_ROI <- TempList[[3]]
        CD274none_ROI <- TempList[[4]]

        ResultPath <- paste0("/mnt/data/lyx/IMC/analysis/spatial/permutation_IM/")
        celltypes <- names(table(sce$SubType))
        celltypes

        ### CD274 high ROIs
        list_ <- getResult(ResultPath, ROIs = CD274high_ROI, celltypes)
        MergeDF1 <- list_[[1]]
        labelDF1 <- list_[[2]]
        rm(list_)

        ### CD274 low ROIs
        list_ <- getResult(ResultPath, ROIs = CD274low_ROI, celltypes)
        MergeDF2 <- list_[[1]]
        labelDF2 <- list_[[2]]
        rm(list_)

        ### CD274 none ROIs
        list_ <- getResult(ResultPath, ROIs = CD274none_ROI, celltypes)
        MergeDF3 <- list_[[1]]
        labelDF3 <- list_[[2]]
        rm(list_)

        # DoubleHeat(MergeDF1, labelDF1, group1 = "CD274pos high", MergeDF2, labelDF2, group2 = "CD274pos low", plot = "groupheatmap", savePath = paste0(savePath, ReclusterName, "_CD274 Interaction Heatmap.pdf"))
        TribleHeat(
            MergeDF1, labelDF1,
            group1 = "CD274pos high",
            MergeDF2, labelDF2, group2 = "CD274pos low",
            MergeDF3, labelDF3, group3 = "CD274pos none",
            savePath = paste0(savePath, structure, "/reclustering/", ReclusterName, "_CD274 Interaction Heatmap.pdf")
        )

        ### plot the cell subtype number different
        sce_Temp1 <- sce[, sce$ID %in% CD274high_ROI]
        CD274HighCountMat1 <- GetAbundance(sce_Temp1, countcol = "SubType", is.fraction = T)
        CD274HighCountMat2 <- GetAbundance(sce_Temp1, countcol = structure, is.fraction = T)

        sce_Temp2 <- sce[, sce$ID %in% CD274low_ROI]
        CD274lowCountMat1 <- GetAbundance(sce_Temp2, countcol = "SubType", is.fraction = T)
        CD274lowCountMat2 <- GetAbundance(sce_Temp2, countcol = structure, is.fraction = T)

        sce_Temp3 <- sce[, sce$ID %in% CD274none_ROI]
        CD274noneCountMat1 <- GetAbundance(sce_Temp3, countcol = "SubType", is.fraction = T)
        CD274noneCountMat2 <- GetAbundance(sce_Temp3, countcol = structure, is.fraction = T)

        #### celltype
        celltypes <- names(table(sce$SubType))
        savePathTemp <- paste0(savePath, structure, "/reclustering/CellAbundance/")
        if (!file.exists(savePathTemp)) {
            dir.create(savePathTemp, recursive = T)
        }
        for (celltype in celltypes) {
            AbundanceSwarmPlot(CD274HighCountMat1, CD274lowCountMat1, CD274noneCountMat1, groupsName = c("high", "low", "none"), celltype = celltype, marker = "CD274", savePath = savePathTemp)
        }

        #### k-means
        k_structures <- names(table(colData(sce)[, structure]))
        savePathTemp <- paste0(savePath, structure, "/reclustering/StructureAbundance/")
        if (!file.exists(savePathTemp)) {
            dir.create(savePathTemp, recursive = T)
        }
        for (k_structure in k_structures) {
            AbundanceSwarmPlot(CD274HighCountMat2, CD274lowCountMat2, CD274noneCountMat2, groupsName = c("high", "low", "none"), celltype = k_structure, marker = "CD274", savePath = savePathTemp)
        }

        #### survival
        CD274CountMat <- GetAbundance(sce_Temp, countcol = "CD274_Type", is.fraction = T, is.reuturnMeans = T)

        SurvivalForPhenoAssoLabel(CD274CountMat, GroupCol = "Pheno_pos", time = "RFS_time", status = "RFS_status", marker = "CD274", cutoffType = "best", savePath = paste0(savePath, structure, "/reclustering/"))
    }
}
## Certain reclustering types in cellular pattern (two group)
if (T) {
    ## Certain reclustering types in cellular pattern
    reclusMarkers <- c(
        # "HLADR", "CD68", "CD14", "CD11c", "CD11b", "CD16", "CLEC9A", "CD169", "CD163", ## Cellidentity
        "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
        "VEGF", "CAIX", ## Hypoxia
        "CD279", "CD274", "CD127", ## Immune-checkpoint
        "CD27", "CD80" ## Immune-activation
    )

    ReMajorType <- c("Myeloid")
    ReclusterName <- "Myeloid"
    for (structure in colName) {
        sce_ <- readRDS(paste0("/mnt/data/lyx/IMC/analysis/structure/", structure, "/Myeloid_recluster.rds"))
        interstType <- c("4")
        PlotCertainTypeinPattern(sce_, Col1 = ReclusterName, types1 = interstType, Col2 = structure, groupCol = "RFS_status", savePath = paste0(savePath, structure, "/reclustering/"))

        ##
        if (F) {
            interstType <- c("4")
            sce__ <- sce_[, sce_$Lymphocyte == interstType]
            table(sce__$MajorType)
            table(sce__$SubType)
            table(sce__$ID)[order(as.numeric(table(sce__$ID)), decreasing = T)]
            phenoLabelCountMat <- GetAbundance(sce_, countcol = "phenoLabel", is.fraction = FALSE, is.reuturnMeans = FALSE)
            phenoLabelCountMat <- TransformIntoPlotMat(phenoLabelCountMat, valueCol = c(1:2))
            head(phenoLabelCountMat)
            BoxPlotForPhenoAssCell(phenoLabelCountMat, savePath = paste0(savePath, structure, "/reclustering/"))
            SurvivalForPhenoAssCell(phenoLabelCountMat, savePath = paste0(savePath, structure, "/reclustering/"))
        }

        ## plot the differential expression genes
        sce__ <- sce_[, sce_$MajorType %in% ReclusterName]
        mat <- as.data.frame(t(assay(sce_)[reclusMarkers, ]))
        mat$phenoLabel <- 0
        mat$phenoLabel[colData(sce_)[, ReclusterName] %in% interstType] <- 1
        FCDF <- FCandPvalueCal(mat, xCol = c(1, length(reclusMarkers)), yCol = (length(reclusMarkers) + 1), need.sample = TRUE)
        VolcanoPlot(FCDF, pthreshold = 0.01, fcthreshold = 3, feature = "Phenotype-Associated cells", filename = paste0(savePath, structure, "/Phenotype-associated differential markers in ", ReMajorType, ".pdf"))

        ## Assign the new label
        ### CD274
        CD274_Pos <- interstType
        TempList <- AssignNewlabel(sce_, allROIs = names(table(sce$ID)), phenoLabel = "CD274_Type", ReclusterName = ReclusterName, interstType = CD274_Pos, cutoffType = "manual", cutoffValue = 6, numGroup = 2)
        sce_Temp <- TempList[[1]]
        CD274high_ROI <- TempList[[2]]
        CD274low_ROI <- TempList[[3]]

        ResultPath <- paste0("/mnt/data/lyx/IMC/analysis/spatial/permutation_IM/")
        celltypes <- names(table(sce$SubType))
        celltypes

        ### CD274 high ROIs
        list_ <- getResult(ResultPath, ROIs = CD274high_ROI, celltypes)
        MergeDF1 <- list_[[1]]
        labelDF1 <- list_[[2]]
        rm(list_)

        ### CD274 low ROIs
        list_ <- getResult(ResultPath, ROIs = CD274low_ROI, celltypes)
        MergeDF2 <- list_[[1]]
        labelDF2 <- list_[[2]]
        rm(list_)

        DoubleHeat(MergeDF1, labelDF1,
            group1 = "CD274pos high",
            MergeDF2, labelDF2, group2 = "CD274pos low", plot = "groupheatmap",
            savePath = paste0(savePath, structure, "/reclustering/", ReclusterName, "_CD274 Interaction Heatmap.pdf")
        )

        ### plot the cell subtype number different
        sce_Temp1 <- sce[, sce$ID %in% CD274high_ROI]
        CD274HighCountMat1 <- GetAbundance(sce_Temp1, countcol = "SubType", is.fraction = T)
        CD274HighCountMat2 <- GetAbundance(sce_Temp1, countcol = structure, is.fraction = T)

        sce_Temp2 <- sce[, sce$ID %in% CD274low_ROI]
        CD274lowCountMat1 <- GetAbundance(sce_Temp2, countcol = "SubType", is.fraction = T)
        CD274lowCountMat2 <- GetAbundance(sce_Temp2, countcol = structure, is.fraction = T)


        #### celltype
        celltypes <- names(table(sce$SubType))
        savePathTemp <- paste0(savePath, structure, "/reclustering/CellAbundance/")
        if (!file.exists(savePathTemp)) {
            dir.create(savePathTemp, recursive = T)
        }
        for (celltype in celltypes) {
            AbundanceSwarmPlot(CD274HighCountMat1, CD274lowCountMat1, groupsName = c("high", "low"), celltype = celltype, marker = "CD274", savePath = savePathTemp, numGroup = 2)
        }

        #### k-means
        k_structures <- names(table(colData(sce)[, structure]))
        savePathTemp <- paste0(savePath, structure, "/reclustering/StructureAbundance/")
        if (!file.exists(savePathTemp)) {
            dir.create(savePathTemp, recursive = T)
        }
        for (k_structure in k_structures) {
            AbundanceSwarmPlot(CD274HighCountMat2, CD274lowCountMat2, groupsName = c("high", "low"), celltype = k_structure, marker = "CD274", savePath = savePathTemp, numGroup = 2)
        }

        #### survival
        CD274CountMat <- GetAbundance(sce_Temp, countcol = "CD274_Type", is.fraction = T, is.reuturnMeans = T)

        SurvivalForPhenoAssoLabel(CD274CountMat, GroupCol = "Pheno_pos", time = "RFS_time", status = "RFS_status", marker = "CD274", cutoffType = "best", savePath = paste0(savePath, structure, "/reclustering/"))
    }
}
