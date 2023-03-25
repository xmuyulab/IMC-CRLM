# Identify tumor microenviroment pattern
library(SingleCellExperiment)
library(pheatmap)
library(survival)
library(survminer)
source("./spatial_analysis_functions.r")
source("./structural_analysis_functions.r")

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
names(table(sce$Tissue))

## omit Tissue associated Tumor
# sce <- sce[, sce$Tissue == "CT" | sce$Tissue == "IM"]
sce <- sce[, sce$Tissue == "IM"]
sce

## TME archetype analysis
if (F) {
    ## load clinical information
    GroupInfo <- load_clinical(sce = sce, clinicalFilePath = "/mnt/data/lyx/IMC/clinical.csv")
    ## merge celltype abundance into a dataframe
    AbundanceDF <- MergeAbundanceResult(sce, return.fraction = T)

    selectCelltypes <- c("DPT", "Mono_CLEC9A", "Macro_Multi", "Macro_CD11b", "SC_Vimentin", "TC_Ki67") ## results from abundance analysis volcano plot

    celltypes <- names(table(sce$SubType))
    celltypes <- celltypes[!celltypes %in% selectCelltypes]

    ## cluster patients
    distMat <- dist(t(AbundanceDF[selectCelltypes, ]), method = "euclidean")
    colclust <- hclust(distMat, method = "complete")

    k <- 15
    TMEClusterID <- cutree(colclust, k = k)
    table(TMEClusterID)

    ## column annotation bar
    plotDF <- AbundanceDF[c(selectCelltypes, celltypes), ]

    annotationCol <- matrix(data = NA, nrow = ncol(plotDF), ncol = 3)
    annotationCol <- as.data.frame(annotationCol)

    rownames(annotationCol) <- colnames(plotDF)
    colnames(annotationCol) <- c("TME Archetypes", "RFSS", "KRAS Mutation")

    annotationCol$`TME Archetypes` <- as.factor(as.numeric(TMEClusterID))
    annotationCol$RFSS <- ifelse(GroupInfo$RFS_status == 1, "Relaps", "Non-Relaps")
    annotationCol$`KRAS Mutation` <- ifelse(GroupInfo$KRAS_mutation == 1, "Mutate", "WT")
    # annotationCol$`CRC Site` <- ifelse(GroupInfo$CRC_site == 1, "Right Colon", ifelse(GroupInfo$CRC_site == 2, "Left Colon", "Rectum"))

    p <- pheatmap(plotDF,
        scale = "column", gaps_row = length(selectCelltypes), cutree_cols = k,
        annotation_col = annotationCol, annotation_legend = TRUE,
        cluster_rows = FALSE, cluster_cols = colclust, clustering_distance_cols = "euclidean", clustering_method = "complete",
        show_rownames = TRUE, show_colnames = TRUE, fontsize_col = 6, fontsize_row = 6,
        color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
        angle_col = "270", cellwidth = 6, cellheight = 6
    )
    pdf("/mnt/data/lyx/IMC/analysis/spatial/TME Archetypes.pdf", height = 5, width = 16)
    print(p)
    dev.off()

    ## Assign TME Archetypes label (Voting)
    PlaeblDF <- as.data.frame(TMEClusterID)
    PlaeblDF$PID <- sapply(rownames(PlaeblDF), function(x) {
        strsplit(x, "_")[[1]][1]
    })
    PIDs <- names(table(PlaeblDF$PID))

    label <- c()
    for (i in PIDs) {
        PlaeblDFTemp <- subset(PlaeblDF, PID == i)
        label <- c(label, names(table(PlaeblDFTemp$TMEClusterID))[1])
    }
    names(label) <- PIDs

    ## survival analysis
    label <- as.data.frame(label)
    label$RFS_time <- GroupInfo$RFS_time[match(rownames(label), GroupInfo$PID)]
    label$RFS_status <- GroupInfo$RFS_status[match(rownames(label), GroupInfo$PID)]
    label <- label[-match("W21", rownames(label)), ]

    km <- survfit(Surv(RFS_time, RFS_status) ~ label, data = label)
    p <- ggsurvplot(km,
        data = label,
        linetype = c("solid", "solid"),
        surv.median.line = "hv", surv.scale = "percent",
        pval = T, risk.table = T,
        conf.int = T, conf.int.alpha = 0.1, conf.int.style = "ribbon",
        risk.table.y.text = T,
        palette = c("#3300CC", "#CC3300"),
        xlab = "Recurrence time"
    )

    pdf("/mnt/data/lyx/IMC/analysis/spatial/TME archetypes survival analysis.pdf", width = 8, height = 6)
    print(p)
    dev.off()

    ## visualize TME archetypes
    savePath <- "/mnt/data/lyx/IMC/analysis/spatial/archetypes/"
    if (!dir.exists(savePath)) {
        dir.create(savePath)
    }
    ROIs <- names(table(sce$ID))
    for (ROI in ROIs) {
        PlotCelltypes(sce, ROI, selectCelltypes, SavePath = paste0(savePath, ROI, " TME archetypes.pdf"))
    }
}

## Cellular neighbors analysis
sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
savePath <- "/mnt/data/lyx/IMC/analysis/structure/"
if (!dir.exists(savePath)) {
    dir.create(savePath, recursive = T)
}

## remain Tissue associated Tumor
sce <- sce[, sce$Tissue == "IM"]
# sce <- sce[, sce$Tissue == "IM" | sce$Tissue == "CT"]
scimapResult <- read.csv("/mnt/data/lyx/IMC/analysis/spatial/cellular_neighbor_IM.csv")
scimapResult <- scimapResult[, -1]

colnames(scimapResult)

colName <- colnames(scimapResult)[14:ncol(scimapResult)]
sce <- BindResult(sce, scimapResult, colName)

colnames(colData(sce))
table(colData(sce)$MajorType)

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
        "HLADR", "CD68", "CD14", "CD11c", "CD11b", "CD16", "CLEC9A", "CD169", "CD163", ## Cellidentity
        "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
        "VEGF", "CAIX", ## Hypoxia
        "CD279", "CD274", "CD127", ## Immune-checkpoint
        "CD27", "CD80" ## Immune-activation

        # "CD45", "CD20", "CD3", "CD8a", "CD4", "FoxP3", "CD57", ## Lymphocyte
        # "HLADR", "CD68", "CD14", "CD11c", "CD11b", "CD16", "CLEC9A", "CD169", "CD163", ## Myeloid
        # "CollagenI", "Vimentin", "AlphaSMA", "FAP", ## Stromal
        # "EpCAM", ## Tumor
        # "Ki67", "HK2", "FASN", "PRPS1", "GLUT1", ## Cell grouth
        # "VEGF", "CAIX", ## Hypoxia
        # "CD127", "CD27", "CD80", ## Immuno-activation
        # "CD274", "CD279", "TIGIT", "CD366" ## Immuno-checkpoint
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

## The difference between re-cluster ids
if (F) {
    reclusterAbun <- GetAbundance(sce_, countcol = ReclusterName, is.fraction = T)
    colnames(reclusterAbun)[1:15] <- sapply(colnames(reclusterAbun)[1:15], function(x) {
        return(paste0("rec_", x))
    })
    reclusterAbun$RFS_status <- as.factor(reclusterAbun$RFS_status)
    # reclusterAbun$rec_5_10 <- reclusterAbun$rec_5 + reclusterAbun$rec_10

    p <- ggplot(reclusterAbun, aes(x = RFS_status, y = rec_5_10, fill = RFS_status)) +
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
        stat_compare_means(aes(group = RFS_status), label.y = 0.4, method = "t.test")

    pdf(paste0(savePath, "knn20_celluarPat/reclustering/rec_", 510, "_abundance_analysis.pdf"), width = 8, height = 6)
    print(p)
    dev.off()
}

## Certain reclustering types in cellular pattern
for (structure in colName) {
    sce_ <- readRDS(paste0("/mnt/data/lyx/IMC/analysis/structure/", structure, "/Myeloid_recluster.rds"))
    interstType <- c("8", "10")
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
    ### CLEC9A
    CLEC9A_Pos <- interstType
    TempList <- AssignNewlabel(sce_, allROIs = names(table(sce$ID)), phenoLabel = "CLEC9A_Type", ReclusterName = ReclusterName, interstType = CLEC9A_Pos, cutoffType = "manual", cutoffValue = 6, numGroup = 3)
    sce_Temp <- TempList[[1]]
    CLEC9Ahigh_ROI <- TempList[[2]]
    CLEC9Alow_ROI <- TempList[[3]]
    CLEC9Anone_ROI <- TempList[[4]]

    ResultPath <- paste0("/mnt/data/lyx/IMC/analysis/spatial/permutation_IM/")
    celltypes <- names(table(sce$SubType))
    celltypes

    ### CLEC9A high ROIs
    list_ <- getResult(ResultPath, ROIs = CLEC9Ahigh_ROI, celltypes)
    MergeDF1 <- list_[[1]]
    labelDF1 <- list_[[2]]
    rm(list_)

    ### CLEC9A low ROIs
    list_ <- getResult(ResultPath, ROIs = CLEC9Alow_ROI, celltypes)
    MergeDF2 <- list_[[1]]
    labelDF2 <- list_[[2]]
    rm(list_)

    ### CLEC9A none ROIs
    list_ <- getResult(ResultPath, ROIs = CLEC9Anone_ROI, celltypes)
    MergeDF3 <- list_[[1]]
    labelDF3 <- list_[[2]]
    rm(list_)

    # DoubleHeat(MergeDF1, labelDF1, group1 = "CLEC9Apos high", MergeDF2, labelDF2, group2 = "CLEC9Apos low", plot = "groupheatmap", savePath = paste0(savePath, ReclusterName, "_CLEC9A Interaction Heatmap.pdf"))
    TribleHeat(
        MergeDF1, labelDF1,
        group1 = "CLEC9Apos high",
        MergeDF2, labelDF2, group2 = "CLEC9Apos low",
        MergeDF3, labelDF3, group3 = "CLEC9Apos none",
        savePath = paste0(savePath, structure, "/reclustering/", ReclusterName, "_CLEC9A Interaction Heatmap.pdf")
    )

    ### plot the cell subtype number different
    sce_Temp1 <- sce[, sce$ID %in% CLEC9Ahigh_ROI]
    CLEC9AHighCountMat1 <- GetAbundance(sce_Temp1, countcol = "SubType", is.fraction = T)
    CLEC9AHighCountMat2 <- GetAbundance(sce_Temp1, countcol = structure, is.fraction = T)

    sce_Temp2 <- sce[, sce$ID %in% CLEC9Alow_ROI]
    CLEC9AlowCountMat1 <- GetAbundance(sce_Temp2, countcol = "SubType", is.fraction = T)
    CLEC9AlowCountMat2 <- GetAbundance(sce_Temp2, countcol = structure, is.fraction = T)

    sce_Temp3 <- sce[, sce$ID %in% CLEC9Anone_ROI]
    CLEC9AnoneCountMat1 <- GetAbundance(sce_Temp3, countcol = "SubType", is.fraction = T)
    CLEC9AnoneCountMat2 <- GetAbundance(sce_Temp3, countcol = structure, is.fraction = T)

    #### celltype
    celltypes <- names(table(sce$SubType))
    savePathTemp <- paste0(savePath, structure, "/reclustering/CellAbundance/")
    if (!file.exists(savePathTemp)) {
        dir.create(savePathTemp, recursive = T)
    }
    for (celltype in celltypes) {
        AbundanceSwarmPlot(CLEC9AHighCountMat1, CLEC9AlowCountMat1, CLEC9AnoneCountMat1, groupsName = c("high", "low", "none"), celltype = celltype, marker = "CLEC9A", savePath = savePathTemp)
    }

    #### k-means
    k_structures <- names(table(colData(sce)[, structure]))
    savePathTemp <- paste0(savePath, structure, "/reclustering/StructureAbundance/")
    if (!file.exists(savePathTemp)) {
        dir.create(savePathTemp, recursive = T)
    }
    for (k_structure in k_structures) {
        AbundanceSwarmPlot(CLEC9AHighCountMat2, CLEC9AlowCountMat2, CLEC9AnoneCountMat2, groupsName = c("high", "low", "none"), celltype = k_structure, marker = "CLEC9A", savePath = savePathTemp)
    }

    #### survival
    CLEC9ACountMat <- GetAbundance(sce_Temp, countcol = "CLEC9A_Type", is.fraction = T, is.reuturnMeans = T)

    SurvivalForPhenoAssoLabel(CLEC9ACountMat, GroupCol = "Pheno_pos", time = "RFS_time", status = "RFS_status", marker = "CLEC9A", cutoffType = "best", savePath = paste0(savePath, structure, "/reclustering/"))
}
