# Fyrther investigate CD163+ Macrophage
library(SingleCellExperiment)

library(pheatmap)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(rstatix)
library(ggcor)

library(survival)
library(survminer)
library(pROC)

library(dplyr)
library(tidyr)

source("./spatial_analysis_functions.r")
source("./structural_analysis_functions.r")

## Load SCE
sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
savePath <- "/mnt/data/lyx/IMC/analysis/reclustering/"
if (!dir.exists(savePath)) {
    dir.create(savePath, recursive = T)
}

## Load clinical
clinical <- read.csv("/mnt/data/lyx/IMC/clinical.csv", header = T)

## remain Tissue associated Tumor
sce <- sce[, sce$Tissue == "IM"]
# sce <- sce[, sce$Tissue == "TAT"]

# sce <- sce[, sce$Tissue == "IM" | sce$Tissue == "CT"]

## Load CNP
scimapResult <- read.csv("/mnt/data/lyx/IMC/analysis/structure/cellular_neighbor_IM.csv")
scimapResult <- scimapResult[, -1]

colnames(scimapResult)

colName <- colnames(scimapResult)[14:ncol(scimapResult)]
sce <- BindResult(sce, scimapResult, colName)

## Differential expressiong markers in Macro_CD163
metaMarkers <- c(
    "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
    "VEGF", "CAIX", ## Hypoxia
    "CD274", "CD279", "TIGIT", "CD366", ## Immuno-checkpoint
    "CD127", "CD27", "CD80" ## Immuno-activation
)

## Same celltypes different in Relapse and non-relapse
HeatmapForMarkersInGroups(sce, metaMarkers, groupCol = "RFS_status", adjust.p = T, FCthreshold = 0, savePath)

## Compare certain types in different tissue
sce <- sce[, sce$Tissue %in% c("IM", "TAT")]
HeatmapForMarkersInGroups(sce, markers = metaMarkers, groupCol = "Tissue", adjust.p = T, FCthreshold = 1.2, savePath)
### Volcano
if (F) {
    sce_ <- sce[, sce$Tissue %in% c("IM", "TAT")]
    sce_ <- sce_[, sce_$SubType %in% c("Treg")]
    mat <- as.data.frame(t(assay(sce_)[metaMarkers, ]))
    mat$phenoLabel <- sce_$Tissue

    mat$phenoLabel <- ifelse(mat$phenoLabel %in% c("IM"), "1", "0")

    FCDF <- FCandPvalueCal(mat, xCol = c(1, length(metaMarkers)), yCol = (length(metaMarkers) + 1), need.sample = TRUE)
    VolcanoPlot(FCDF, pthreshold = 0.05, fcthreshold = 1.4, feature = "Treg in different tissue", filename = paste0("/mnt/data/lyx/IMC/analysis/reclustering/Treg markers diff in IM and TAT.pdf"), Qvalue = TRUE)
}

## Treg & Mono_CD11c - Reclustering
anajorTypes <- c("Myeloid", "Lymphocyte")
anaclusterNames <- c("Mono_CD11c", "Treg")
structure <- "CNP20"

for (i in 1:length(anajorTypes)) {
    MajorName <- anajorTypes[i]
    SubName <- anaclusterNames[i]

    savePathTemp1 <- paste0(savePath, SubName, "/")
    if (!dir.exists(savePathTemp1)) {
        dir.create(savePathTemp1, recursive = T)
    }

    ## Reclustering
    sce_ <- Reclustering(sce, markers = metaMarkers, ReMajorType = MajorName, ReclusterName = SubName, ReSubType = SubName, PatternCol = structure, ncluster = 15, savePath = savePathTemp1)

    ## Cluster difference
    reclusterAbun <- GetAbundance(sce_, countcol = SubName, is.fraction = T)
    colnames(reclusterAbun)[1:15] <- sapply(colnames(reclusterAbun)[1:15], function(x) {
        return(paste0("rec_", x))
    })
    reclusterAbun$RFS_status <- as.factor(reclusterAbun$RFS_status)
    plotdf <- as.data.frame(matrix(data = NA, nrow = 15 * nrow(reclusterAbun), ncol = 3))
    plotdf[, 1] <- rep(colnames(reclusterAbun)[1:15], each = nrow(reclusterAbun))
    plotdf[, 2] <- as.numeric(as.matrix(reclusterAbun[, 1:15]))
    plotdf[, 3] <- rep(reclusterAbun$RFS_status, times = 15)

    colnames(plotdf) <- c("ReCluster", "Abundance", "Group")

    if (F) {
        plotdf <- plotdf[plotdf$ReCluster %in% c("rec_13", "rec_8", "rec_5", "rec_12"), ]
        rownames(plotdf) <- 1:nrow(plotdf)
    }

    p <- ggplot(plotdf, aes(x = Group, y = Abundance, fill = Group)) +
        geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
        scale_y_continuous(name = "Cell Abundance") +
        scale_x_discrete(name = "Cell Population") +
        scale_fill_manual(values = c("#56B4E9", "#D55E00")) +
        theme_bw() +
        theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            text = element_text(size = 12),
            axis.title = element_text(face = "bold"),
            axis.text.x = element_text(size = 11, angle = 90),
            legend.position = "none",
            strip.text = element_text(size = 12, face = "bold"),
            strip.background = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            panel.grid.major = element_line(color = "gray", size = 0.5),
            panel.grid.minor = element_blank()
        ) +
        stat_compare_means(aes(group = Group), label.y = .4, method = "t.test") +
        facet_wrap(~ReCluster, ncol = 4)

    # pdf(paste0(savePathTemp1, "recluster_abundance_analysis(part).pdf"), width = 8, height = 4)
    pdf(paste0(savePathTemp1, "recluster_abundance_analysis.pdf"), width = 8, height = 6)
    print(p)
    dev.off()

    saveRDS(sce_, paste0(savePathTemp1, SubName, "_recluster.rds"))
}

interstTypeList <- list()
interstTypeList[[1]] <- c("1", "7") ## PD-L1+ CD11c+ Monocytes
interstTypeList[[2]] <- c("5", "8") ## PD-1+ Foxp3+ Lymphocytes

for (i in 1:length(anajorTypes)) {
    MajorName <- anajorTypes[i]
    SubName <- anaclusterNames[i]

    savePathTemp1 <- paste0(savePath, SubName, "/")
    if (!dir.exists(savePathTemp1)) {
        dir.create(savePathTemp1, recursive = T)
    }
    sce_ <- readRDS(paste0("/mnt/data/lyx/IMC/analysis/reclustering/", SubName, "/", SubName, "_recluster.rds"))
    ## Assign the new label
    interstType <- interstTypeList[[i]]

    TempList <- AssignNewlabel(sce_,
        allROIs = names(table(sce$ID)), phenoLabel = paste0(SubName, "_Type"), ReclusterName = SubName, interstType = interstType,
        clinicalFeatures = c("RFS_time", "RFS_status"), cutoffType = "mean", cutoffValue = 20, numGroup = 2
    )
    sce_Temp <- TempList[[1]]
    high_ROI <- TempList[[2]]
    low_ROI <- TempList[[3]]

    if (F) {
        ## spatial interaction analysis
        sce_Temp <- sce_Temp[, sce_Temp$Treg_Type == "Pheno_pos"]
        sce$SubType[match(sce_Temp$CellID, sce$CellID)] <- "Acitve Treg"

        library(spatstat)
        ROIs <- names(table(sce$ID))

        # Define the whole image size
        xlim <- c(0, 1000)
        ylim <- c(0, 1000)

        # Record the 10-nn cells of PD-L1+ CD11c+ Monocytes
        nn10_NeighborsList <- list()

        for (ROI in ROIs) {
            sce_ <- sce[, sce$ID == ROI]

            coor <- sce_$Position
            coor_df <- GetXandY(coor)
            coor_df$SubType <- sce_$SubType

            head(coor_df)

            # Convert your data frame to a ppp object
            points_ppp <- ppp(coor_df$cor_x, coor_df$cor_y, xrange = xlim, yrange = ylim, marks = coor_df$SubType)

            # Find nearest neighbors
            nn <- nncross(points_ppp, points_ppp, k = 2:11)

            nn10_neighbors <- nn[, 11:ncol(nn)]
            nn10_neighbors <- apply(nn10_neighbors, MARGIN = 2, function(x) {
                x <- as.numeric(x)
                return(sce_$CellID[x])
            })

            PDL1CD11c_idx <- (sce_$SubType == "Acitve Treg")
            nn10_neighbors_ <- nn10_neighbors[PDL1CD11c_idx, ]
            # table(nn10_neighbors_)

            nn10_NeighborsList[[ROI]] <- as.numeric(nn10_neighbors_)
        }
        Subtypes <- names(table(sce$SubType))
        InterDF <- as.data.frame(matrix(data = NA, nrow = length(nn10_NeighborsList), ncol = length(Subtypes)))
        colnames(InterDF) <- Subtypes

        for (i in 1:length(nn10_NeighborsList)) {
            CellIDTemp <- nn10_NeighborsList[[i]]
            CelltypeTemp <- sce$SubType[match(CellIDTemp, sce$CellID)]
            CelltypeTemp <- factor(CelltypeTemp, levels = Subtypes)
            CelltypeTemp <- as.data.frame(table(CelltypeTemp))
            InterDF[i, ] <- CelltypeTemp[, 2]
        }
        rownames(InterDF) <- names(nn10_NeighborsList)
        InterDF$PID <- sapply(rownames(InterDF), function(x) {
            return(strsplit(x, "_")[[1]][1])
        })
        InterDF$RFS_status <- sce$RFS_status[match(InterDF$PID, sce$PID)]

        plotdf <- tidyr::pivot_longer(InterDF, 1:22, names_to = "Neighbors", values_to = "Number")
        plotdf$RFS_status <- as.factor(plotdf$RFS_status)
        color <- ggsci::pal_aaas("default")(10)
        p <- ggplot(data = plotdf, aes(x = Neighbors, y = Number, fill = RFS_status)) +
            geom_boxplot(show.legend = FALSE, outlier.shape = NA, alpha = 0.8) +
            theme_bw() +
            labs(x = NULL, y = "Count", title = "Box Plot of Counts by Group and Cell Type") +
            stat_compare_means(
                aes(group = RFS_status),
                method = "t.test",
                label = "p.signif",
                hide.ns = T,
                size = 4
            ) +
            theme(
                plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                axis.title = element_text(size = 12, face = "bold"),
                axis.text = element_text(size = 10),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
                strip.text = element_text(size = 12, face = "bold"),
                strip.background = element_blank(),
                panel.border = element_rect(color = "black", fill = NA, size = 1),
                panel.grid.major = element_line(color = "gray", size = 0.5),
                panel.grid.minor = element_blank()
            )
        pdf("test.pdf", height = 6, width = 10)
        print(p)
        dev.off()
        if (F) {
            allNeighbors <- unlist(nn10_NeighborsList)
            allNeighbors_SubType <- sce$SubType[allNeighbors]

            ## Compare CD8T
            CD8T_idx <- sce$CellID[(sce$SubType == "CD8T")]
            CD8T_seu <- sce[, sce$CellID %in% CD8T_idx]

            Interact_CD8T_idx <- as.numeric(intersect(CD8T_idx, allNeighbors))
            Remain_CD8T_idx <- as.numeric(setdiff(CD8T_idx, Interact_CD8T_idx))

            CD8T_seu$label <- "0"
            CD8T_seu$label[match(Interact_CD8T_idx, CD8T_seu$CellID)] <- "1"

            mat <- as.data.frame(t(assay(CD8T_seu)))
            mat$label <- CD8T_seu$label

            xCol <- c(1, 35)
            yCol <- 36
            mat_foldchangeMat <- FCandPvalueCal(mat, xCol = xCol, yCol = yCol)
            mat_foldchangeMat$Q.value <- p.adjust(mat_foldchangeMat$P.value, method = "BH")
            mat_foldchangeMat$P.value <- mat_foldchangeMat$Q.value

            VolcanoPlot(mat_foldchangeMat, pthreshold = 0.05, fcthreshold = 1.2, feature = NULL, filename = "testC8T.pdf")
        }

        ## Take mean
        MeanDensity <- as.data.frame(matrix(data = NA, nrow = length(DensityList), ncol = 2))

        for (i in 1:length(DensityList)) {
            MeanDensity[i, ] <- apply(DensityList[[i]], MARGIN = 2, FUN = "mean")
        }
        rownames(MeanDensity) <- names(DensityList)
        colnames(MeanDensity) <- colnames(DensityList[[1]])
        MeanDensity$SL <- MeanDensity[, 1] / MeanDensity[, 2]

        MeanDensity$PID <- sapply(rownames(MeanDensity), function(x) {
            return(strsplit(x, "_")[[1]][1])
        })
        MeanDensity$RFS_time <- clinical$RFS_time[match(MeanDensity$PID, clinical$PID)]
        MeanDensity$RFS_status <- clinical$RFS_status[match(MeanDensity$PID, clinical$PID)]

        MeanDensity <- MeanDensity %>%
            group_by(PID) %>%
            summarise(across(c(1:5), mean, na.rm = TRUE))

        cutpoint <- surv_cutpoint(data = MeanDensity, time = "RFS_time", event = "RFS_status", variables = "SL")
        cutpoint <- summary(cutpoint)$cutpoint

        df <- MeanDensity[, c("SL", "RFS_status", "RFS_time")]
        df[, 1] <- ifelse(df[, 1] >= cutpoint, "high", "low")

        df$RFS_status <- as.numeric(df$RFS_status)
        df$RFS_time <- as.numeric(df$RFS_time)

        ## km curve
        fit <- survfit(Surv(RFS_time, RFS_status) ~ SL, data = df)
        p <- ggsurvplot(fit,
            data = df,
            linetype = c("solid", "solid"),
            surv.median.line = "hv", surv.scale = "percent",
            pval = T, risk.table = T,
            conf.int = T, conf.int.alpha = 0.1, conf.int.style = "ribbon",
            risk.table.y.text = T,
            palette = c("#3300CC", "#CC3300"),
            xlab = "Recurrence time"
        )

        pdf("KM for Stromal:Lymphocyte density (Patient).pdf", width = 8, height = 6)
        print(p)
        dev.off()
    }

    # colData(sce)[match(colData(sce_Temp)[colData(sce_Temp)$Mono_CD11c_Type=="Pheno_pos", ]$CellID, colData(sce)$CellID), "SubType"] <- "PDL1+ Mono_CD11c"
    # colData(sce)[match(colData(sce_Temp)[colData(sce_Temp)$Treg_Type=="Pheno_pos", ]$CellID, colData(sce)$CellID), "SubType"] <- "PD1+ Treg"
    # table(sce$SubType)
    ## saveRDS(sce,"/mnt/data/lyx/IMC/analysis/allsce(new).rds")

    ## marker difference
    sce__ <- sce_[, sce_$SubType %in% SubName]
    mat <- as.data.frame(t(assay(sce_)[metaMarkers, ]))
    mat$phenoLabel <- 0
    mat$phenoLabel[colData(sce_)[, SubName] %in% interstType] <- 1
    FCDF <- FCandPvalueCal(mat, xCol = c(1, length(metaMarkers)), yCol = (length(metaMarkers) + 1), need.sample = TRUE)
    VolcanoPlot(FCDF, pthreshold = 0.01, fcthreshold = 3, feature = "Phenotype-Associated cells", filename = paste0(savePathTemp1, "Phenotype-associated differential markers in ", SubName, ".pdf"), Qvalue = TRUE)

    ## spatial difference
    ResultPath <- paste0("/mnt/data/lyx/IMC/analysis/spatial/permutation_IM/")
    celltypes <- names(table(sce$SubType))
    celltypes

    ### high ROIs
    list_ <- getResult(ResultPath, ROIs = high_ROI, celltypes)
    MergeDF1 <- list_[[1]]
    labelDF1 <- list_[[2]]
    rm(list_)

    ### low ROIs
    list_ <- getResult(ResultPath, ROIs = low_ROI, celltypes)
    MergeDF2 <- list_[[1]]
    labelDF2 <- list_[[2]]
    rm(list_)

    DoubleHeat(MergeDF1, labelDF1,
        group1 = paste0(SubName, " high"),
        MergeDF2, labelDF2, group2 = paste0(SubName, " low"), plot = "groupheatmap",
        savePath = paste0(savePathTemp1, SubName, " Interaction Heatmap.pdf")
    )
    ### Different
    Sig <- 0.01
    getInteracDiff(ResultPath, sce, celltypes = celltypes, savepath = paste0(savePathTemp1, SubName, "_Spatial_Diff_", Sig, ".pdf"), IDs1 = high_ROI, IDs2 = low_ROI, Sig = Sig)

    ### plot the cell subtype number different
    sce_Temp1 <- sce[, sce$ID %in% high_ROI]
    HighCountMat1 <- GetAbundance(sce_Temp1, countcol = "SubType", is.fraction = T)
    HighCountMat2 <- GetAbundance(sce_Temp1, countcol = structure, is.fraction = T)

    sce_Temp2 <- sce[, sce$ID %in% low_ROI]
    lowCountMat1 <- GetAbundance(sce_Temp2, countcol = "SubType", is.fraction = T)
    lowCountMat2 <- GetAbundance(sce_Temp2, countcol = structure, is.fraction = T)


    #### celltype
    celltypes <- names(table(sce$SubType))
    savePathTemp <- paste0(savePathTemp1, "/CellAbundance/")
    if (!file.exists(savePathTemp)) {
        dir.create(savePathTemp, recursive = T)
    }
    for (celltype in celltypes) {
        AbundanceSwarmPlot(HighCountMat1, lowCountMat1, groupsName = c("high", "low"), celltype = celltype, marker = SubName, savePath = savePathTemp, numGroup = 2)
    }

    #### Multiple cell subpopulations abudnace plot
    PlotTypes <- c("CD8T", "Macro_HLADR", "Macro_CD11b", "Macro_CD163", "Macro_CD169", "Treg")
    p <- MultipleAbundanceSwarmPlot(HighCountMat1, lowCountMat1, groupsName = c("high", "low"), PlotTypes = PlotTypes, marker = SubName, numGroup = 2, style = "box")
    pdf(paste0(savePathTemp1, "ViolinPlot of multiple celltypes in ", marker, " group.pdf"), height = 10, width = 8)
    print(p)
    dev.off()

    #### CNPs
    k_structures <- names(table(colData(sce)[, structure]))
    savePathTemp <- paste0(savePathTemp1, "/StructureAbundance/")
    if (!file.exists(savePathTemp)) {
        dir.create(savePathTemp, recursive = T)
    }
    for (k_structure in k_structures) {
        AbundanceSwarmPlot(HighCountMat2, lowCountMat2, groupsName = c("high", "low"), celltype = k_structure, marker = SubName, savePath = savePathTemp, numGroup = 2)
    }

    #### survival
    CountMat <- GetAbundance(sce_Temp, countcol = paste0(SubName, "_Type"), is.fraction = T, is.reuturnMeans = F)

    SurvivalForPhenoAssoLabel(CountMat, GroupCol = "Pheno_pos", time = "RFS_time", status = "RFS_status", marker = SubName, cutoffType = "best", savePath = savePathTemp1)
}

## CD163+ Macrophage analysis
if (T) {
    SubName <- "Macro_CD163"
    interstType <- "Macro_CD163"

    savePathTemp1 <- paste0(savePath, SubName, "/")
    if (!dir.exists(savePathTemp1)) {
        dir.create(savePathTemp1, recursive = T)
    }

    ## Assign the new label
    TempList <- AssignNewlabel(sce, allROIs = names(table(sce$ID)), phenoLabel = paste0(SubName, "_Label"), ReclusterName = "SubType", interstType = interstType, cutoffType = "mean", cutoffValue = 20, numGroup = 2)
    sceTemp <- TempList[[1]]
    high_ROI <- TempList[[2]]
    low_ROI <- TempList[[3]]

    ## spatial difference
    ResultPath <- paste0("/mnt/data/lyx/IMC/analysis/spatial/permutation_IM/")
    celltypes <- names(table(sce$SubType))
    celltypes

    ### high ROIs
    list_ <- getResult(ResultPath, ROIs = high_ROI, celltypes)
    MergeDF1 <- list_[[1]]
    labelDF1 <- list_[[2]]
    rm(list_)

    ### low ROIs
    list_ <- getResult(ResultPath, ROIs = low_ROI, celltypes)
    MergeDF2 <- list_[[1]]
    labelDF2 <- list_[[2]]
    rm(list_)

    DoubleHeat(MergeDF1, labelDF1,
        group1 = paste0(SubName, " high"),
        MergeDF2, labelDF2, group2 = paste0(SubName, " low"), plot = "groupheatmap",
        savePath = paste0(savePathTemp1, SubName, " Interaction Heatmap.pdf")
    )
    ### Different
    getInteracDiff(ResultPath, sce, celltypes = celltypes, savepath = paste0(savePathTemp1, SubName, "_Spatial_Diff_0.01.pdf"), IDs1 = high_ROI, IDs2 = low_ROI, Sig = 0.01)

    PlotTypes <- c("B", "CD4T", "CD8T", "Treg", "Macro_CD163", "Macro_CD169", "Mono_Classic", "Mono_Intermediate")
    p <- InterDiffInvertBarplot(
        ResultPath = ResultPath, sce = sce, PlotTypes = PlotTypes,
        IDs1 = high_ROI, IDs2 = low_ROI
    )
    pdf(paste0(savePathTemp1, SubName, "_Spatial_Diff_Barplot.pdf"), height = 6, width = 8)
    print(p)
    dev.off()

    p <- InterAbundanceBarplot(
        ResultPath = ResultPath, sce = sce, PlotTypes = PlotTypes, groupsName = c("High", "Low"), IDs1 = high_ROI, IDs2 = low_ROI
    )
    pdf(paste0(savePathTemp1, SubName, "_SpatialAbundance.pdf"), height = 6, width = 10)
    print(p)
    dev.off()

    ### plot the cell subtype number different
    sceTemp1 <- sce[, sce$ID %in% high_ROI]
    HighCountMat1 <- GetAbundance(sceTemp1, countcol = "SubType", is.fraction = T)
    HighCountMat2 <- GetAbundance(sceTemp1, countcol = structure, is.fraction = T)

    sceTemp2 <- sce[, sce$ID %in% low_ROI]
    lowCountMat1 <- GetAbundance(sceTemp2, countcol = "SubType", is.fraction = T)
    lowCountMat2 <- GetAbundance(sceTemp2, countcol = structure, is.fraction = T)


    #### celltype
    celltypes <- names(table(sce$SubType))
    savePathTemp <- paste0(savePathTemp1, "CellAbundance/")
    if (!file.exists(savePathTemp)) {
        dir.create(savePathTemp, recursive = T)
    }
    for (celltype in celltypes) {
        AbundanceSwarmPlot(HighCountMat1, lowCountMat1, groupsName = c("high", "low"), celltype = celltype, marker = SubName, savePath = savePathTemp, numGroup = 2)
    }

    #### Multiple cell subpopulations abudnace plot
    PlotTypes <- c("B", "CD8T", "Macro_CD11b", "Macro_CD163", "Macro_CD169", "Treg")
    p <- MultipleAbundanceSwarmPlot(HighCountMat1, lowCountMat1, groupsName = c("high", "low"), PlotTypes = PlotTypes, marker = SubName, numGroup = 2, style = "violin")
    pdf(paste0(savePathTemp1, "ViolinPlot of multiple celltypes in ", marker, " group.pdf"), height = 10, width = 8)
    print(p)
    dev.off()

    #### CNPs
    k_structures <- names(table(colData(sce)[, structure]))
    savePathTemp <- paste0(savePathTemp1, "StructureAbundance/")
    if (!file.exists(savePathTemp)) {
        dir.create(savePathTemp, recursive = T)
    }
    for (k_structure in k_structures) {
        AbundanceSwarmPlot(HighCountMat2, lowCountMat2, groupsName = c("high", "low"), celltype = k_structure, marker = SubName, savePath = savePathTemp, numGroup = 2)
    }

    #### survival
    CountMat <- GetAbundance(sceTemp, countcol = paste0(SubName, "_Label"), is.fraction = T, is.reuturnMeans = F)

    SurvivalForPhenoAssoLabel(CountMat, GroupCol = "Pheno_pos", time = "RFS_time", status = "RFS_status", marker = SubName, cutoffType = "best", savePath = savePathTemp1)
}

## Treg analysis - 2
## Hypothesis: The PD-1+ Treg only in IM but not in TAT
## Treg Reclustering
if (T) {
    TargeType <- "Treg"
    savePathTemp1 <- paste0(savePath, TargeType, "/")
    if (!dir.exists(savePathTemp1)) {
        dir.create(savePathTemp1, recursive = T)
    }

    metaMarkers <- c(
        "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
        "VEGF", "CAIX", ## Hypoxia
        "CD279", "TIGIT", "CD366", ## Immuno-checkpoint
        "CD27" ## Immuno-activation
    )

    sce <- sce[, sce$Tissue %in% c("IM", "CT", "TAT")]
    sce_Type <- sce[, sce$SubType %in% TargeType]

    ## Dimention Redcution
    m <- DimenRec(sce_ = sce_Type, GroupFeature = "RFS_status", markers = metaMarkers, do.pca = FALSE, pca.dimen = 10, explaine.pca = F, embeeding.method = "tsne", num.threads = 8)

    ## calculate Fraction and Entropy via KNN graph
    k <- 20
    cordim <- 2

    resultdf <- KGandFracEntro(m, cordim = cordim, k = k, central.weight = 2, multi.cores = 8)

    ## Generate Fraction Null Distribution
    fraction_NullDistribution <- GeneNullDistri(m, cordim = cordim, sample.size = k, perm.times = 2000)

    ## Quantile of null distribution
    threshold <- 0.05
    low.95 <- as.numeric(quantile(fraction_NullDistribution, threshold))
    high.95 <- as.numeric(quantile(fraction_NullDistribution, (1 - threshold)))

    trueFraction <- resultdf[, 2]
    Class <- ifelse(trueFraction <= low.95, "Pheno_neg", ifelse(
        trueFraction >= high.95, "Pheno_pos", "Background"
    ))

    table(Class)

    ## Visualize class label and entropy
    PlotClassLabel(
        sce_ = sce_Type, Classlabels = Class, Entropylabels = resultdf[, 1], markers = metaMarkers, sample.size = 1e4, num.threads = 8,
        SavePath1 = paste0(savePathTemp1, "T-SNE of phenotype associated label.pdf"), SavePath2 = paste0(savePathTemp1, "T-SNE of phenotype associated entropy.pdf"),
        SavePath3 = savePathTemp1, seed = 619
    )

    ## Volcano for Pheno_pos and Pheno_neg cells

    ### Relapse associated
    mat <- as.data.frame(t(assay(sce_Type)))
    mat <- mat[, match(metaMarkers, colnames(mat))]
    mat$label <- Class

    # mat <- mat[!(mat$label%in%"Background"),]
    mat$label <- ifelse(mat$label == "Pheno_pos", 1, 0)

    mat_foldchangeMat <- FCandPvalueCal(mat = mat, xCol = c(1, length(metaMarkers)), yCol = (length(metaMarkers) + 1), need.sample = T, sample.size = 5e3)
    mat_foldchangeMat$Q.value <- p.adjust(mat_foldchangeMat$P.value, method = "BH")
    mat_foldchangeMat$P.value <- mat_foldchangeMat$Q.value

    VolcanoPlot(mat_foldchangeMat, pthreshold = 0.05, fcthreshold = 4, feature = NULL, filename = paste0(savePathTemp1, "DEGs of relapse associated Tregs.pdf"))

    ### Non relapse associated
    mat <- as.data.frame(t(assay(sce_Type)))
    mat <- mat[, match(metaMarkers, colnames(mat))]
    mat$label <- Class

    mat$label <- ifelse(mat$label == "Pheno_neg", 1, 0)

    mat_foldchangeMat <- FCandPvalueCal(mat = mat, xCol = c(1, length(metaMarkers)), yCol = (length(metaMarkers) + 1), need.sample = T, sample.size = 5e3)
    mat_foldchangeMat$Q.value <- p.adjust(mat_foldchangeMat$P.value, method = "BH")
    mat_foldchangeMat$P.value <- mat_foldchangeMat$Q.value

    VolcanoPlot(mat_foldchangeMat, pthreshold = 0.05, fcthreshold = 2, feature = NULL, filename = paste0(savePathTemp1, "DEGs of non-relapse associated Tregs.pdf"))

    sce_Type$Treg_Type <- Class

    saveRDS(sce_Type, paste0(savePathTemp1, "Treg_recluster.rds"))
}

## investigate the tissue specificity of PD-1+ CD27+ Treg
if (T) {
    sce_ <- readRDS(paste0(savePathTemp1, "Treg_recluster.rds"))

    ## Cluster difference in Tissue
    plotdf2 <- GetAbundance(sce_, countcol = "Treg_Type", clinicalFeatures = c("Tissue", "RFS_status"), is.fraction = T)
    plotdf2 <- plotdf2[, c("Pheno_pos", "PID", "Tissue", "RFS_status")]

    if (F) {
        write.table(plotdf2, "Active Treg Abundance for model construction.csv", sep = ",", row.names = T, col.names = T)
    }

    plotdf2 <- plotdf2[plotdf2$Tissue %in% c("IM", "CT", "TAT"), ]
    colnames(plotdf2)[1] <- "Abundance"

    p <- ggplot(plotdf2, aes(x = Tissue, y = Abundance, fill = "#56B4E9")) +
        geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
        scale_y_continuous(name = "Cell Abundance", limits = c(0, 0.6)) +
        scale_x_discrete(name = "Cell Population") +
        theme_bw() +
        theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            text = element_text(size = 12),
            axis.title = element_text(face = "bold"),
            axis.text.x = element_text(size = 11, angle = 90),
            legend.position = "none",
            strip.text = element_text(size = 12, face = "bold"),
            strip.background = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            panel.grid.major = element_line(color = "gray", size = 0.5),
            panel.grid.minor = element_blank()
        )

    stat.test <- plotdf2 %>%
        t_test(Abundance ~ Tissue) %>%
        add_significance() %>%
        add_xy_position(fun = "mean_sd", scales = "free_y")

    p1 <- p + stat_pvalue_manual(stat.test, hide.ns = FALSE, label = "{p}", step.increase = 0.08)

    pdf(paste0(savePathTemp1, "Pheno-associated Treg tissue distribution.pdf"), width = 4, height = 3)
    print(p1)
    dev.off()
}

## Plot the location of PD-1+ CD27+ Treg
if (T) {
    sce <- sce[, sce$Tissue == "IM"]
    sce_ <- readRDS(paste0(savePathTemp1, "Treg_recluster.rds"))

    head(colData(sce_))
    sce__ <- sce_[, sce_$Treg_Type == "Pheno_pos"]
    sce__ <- sce__[, sce__$Tissue == "IM"]

    sort(table(sce__$ID), decreasing = T)

    source("./spatial_analysis_functions.r")

    PD1TregID <- sce__$CellID
    sce[, match(PD1TregID, sce$CellID)]$SubType <- "CD27CD279+ Treg"

    ROI_num <- sort(table(sce__$ID), decreasing = T)
    ROIs <- names(ROI_num[ROI_num > 20])

    celltype <- c("CD27CD279+ Treg")
    for (ROI in ROIs) {
        VisualizeMarkerAndType(sce = sce, ROI = ROI, celltype = celltype, whichType = "SubType", marker = NULL, SavePath = savePathTemp1)
    }
}

## investigate the spatial interaction of PD-1+ CD27+ Treg
if (T) {
    library(SingleCellExperiment)
    library(spatstat)
    source("./structural_analysis_functions.r")

    TargeType <- "Treg"
    savePathTemp1 <- paste0(savePath, TargeType, "/")

    sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
    sce <- sce[, sce$Tissue %in% c("IM", "TAT")]

    sce_ <- readRDS(paste0("/mnt/data/lyx/IMC/analysis/reclustering/Treg/Treg_recluster.rds"))
    sce_ <- sce_[, sce_$Tissue %in% c("IM", "TAT")]
    sce_ <- sce_[, sce_$Treg_Type %in% c("Pheno_pos")]

    sce[, match(sce_$CellID, sce$CellID)]$SubType <- "PD-1+ CD27+ Treg"
    table(sce$SubType)

    metaMarkers <- c(
        "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
        "VEGF", "CAIX", ## Hypoxia
        "CD279", "TIGIT", "CD366", ## Immuno-checkpoint
        "CD27" ## Immuno-activation
    )

    ## spatial interaction analysis
    ROIs <- names(table(sce$ID))

    # Define the whole image size
    k_NeighborsList <- CalKNNBySpatstat(sce, ROIs = ROIs, CenType = "PD-1+ CD27+ Treg", NeighType = "SubType", k = 10, xlim = c(0, 1000), ylim = c(0, 1000), return.cellID = TRUE)

    ## Combine all neighbors
    neighbors <- c()

    for (i in k_NeighborsList) {
        neighbors <- c(neighbors, i)
    }

    neighbors <- neighbors[(neighbors %in% sce$CellID)]

    ## Calculate the neighbors type abudance
    sce_neighbor <- sce[, match(neighbors, sce$CellID)]
    neiborAbun <- as.data.frame(table(sce_neighbor$SubType))

    ### Visualize neighbors abudance
    if (T) {
        # Rename the columns
        colnames(neiborAbun) <- c("category", "counts")

        # Calculate the percentage for each category
        neiborAbun$percentage <- neiborAbun$counts / sum(neiborAbun$counts)

        # Define labels for the pie chart
        neiborAbun$label <- paste0(neiborAbun$category, "\n(", round(neiborAbun$percentage * 100, 2), "%)")

        # Sort the dataframe by percentage (in descending order)
        neiborAbun <- neiborAbun[order(-neiborAbun$percentage), ]

        # Select top k categories for labelling on the chart
        k <- 5
        neiborAbun$label[(k + 1):nrow(neiborAbun)] <- ""

        # Drop unused factor levels in the 'category' column
        neiborAbun$category <- droplevels(neiborAbun$category)

        # Create the pie chart
        p <- ggplot(neiborAbun, aes(x = "", y = percentage, fill = category)) +
            geom_bar(width = 1, stat = "identity", color = "black") +
            geom_label_repel(aes(label = ifelse(percentage > 0.08, label, "")), # labels for slices < 1% are set to ""
                nudge_y = 0.1, size = 4, show.legend = FALSE
            ) +
            coord_polar("y", start = 0) +
            theme_bw() +
            scale_fill_manual(values = ggsci::pal_ucscgb("default")(length(levels(neiborAbun$category)))) +
            theme(
                panel.grid.major = element_line(color = "gray", size = 0.2),
                panel.grid.minor = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                plot.title = element_text(size = 18, hjust = 0.5),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 12),
                legend.position = "bottom"
            )


        pdf(paste0(savePathTemp1, "CD27+ PD-1+ Treg neighbors abundance.pdf"), width = 10, height = 6)
        print(p)
        dev.off()
    }

    ## Calcualte the interaction strength
    if (T) {
        InterStrenDF <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))

        for (ROI in names(k_NeighborsList)) {
            sceTemp <- sce[, sce$ID == ROI]
            NeiborTemp <- k_NeighborsList[[ROI]]

            sce_neighborTemp <- sce[, match(NeiborTemp, sce$CellID)]
            neiborAbunTemp <- as.data.frame(table(sce_neighborTemp$SubType))

            neiborAbunTemp[, 2] <- neiborAbunTemp[, 2] / (length(NeiborTemp) / 10)
            neiborAbunTemp[, 3] <- rep(ROI, nrow(neiborAbunTemp))

            InterStrenDF <- rbind(InterStrenDF, neiborAbunTemp)
        }

        colnames(InterStrenDF) <- c("NeighborType", "InteractionStrenth", "ROI")
        InterStrenDF$Tissue <- colData(sce)[match(InterStrenDF$ROI, sce$ID), "Tissue"]

        ## Visualization
        if (T) {
            interstType <- c("Macro_CD163", "Mono_CD11c", "Treg", "NK", "CD8T", "B", "CD4T")

            plotdf <- InterStrenDF[InterStrenDF$NeighborType %in% interstType, ]

            p <- ggplot(plotdf, aes(x = NeighborType, y = InteractionStrenth, fill = Tissue)) +
                geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
                scale_y_continuous(name = "Interaction Strength", limits = c(0, 3.5)) +
                scale_x_discrete(name = "Neighbor Type") +
                scale_fill_lancet() +
                stat_compare_means(aes(group = Tissue),
                    method = "t.test",
                    hide.ns = TRUE, # Hide non-significant comparisons
                    label = "p.signif",
                    label.y.npc = "top"
                ) +
                theme_bw() +
                theme(
                    plot.title = element_text(size = 14, face = "bold"),
                    text = element_text(size = 12),
                    axis.title = element_text(face = "bold"),
                    axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5),
                    strip.background = element_blank(), # Remove facet grid background
                    panel.background = element_blank(), # Make the panel background empty
                    panel.grid.major = element_blank(), # Remove major grid lines
                    panel.grid.minor = element_blank(), # Remove minor grid lines
                    panel.border = element_blank() # Remove borders around the facets
                )


            pdf(paste0(savePathTemp1, "CD27+ PD-1+ Treg interaction strength between tissue.pdf"))
            print(p)
            dev.off()
        }
    }

    ## Active Treg and Normal Treg neighbors fraction comparison
    if (T) {
        ATreg_k_NeighborsList <- k_NeighborsList
        NTreg_k_NeighborsList <- CalKNNBySpatstat(sce, ROIs = ROIs, CenType = "Treg", NeighType = "SubType", k = 10, xlim = c(0, 1000), ylim = c(0, 1000), return.cellID = TRUE)

        ## Calcualte the interaction strength
        if (T) {
            ## Active Treg
            AInterStrenDF <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))

            for (ROI in names(ATreg_k_NeighborsList)) {
                sceTemp <- sce[, sce$ID == ROI]
                NeiborTemp <- ATreg_k_NeighborsList[[ROI]]

                sce_neighborTemp <- sce[, match(NeiborTemp, sce$CellID)]
                neiborAbunTemp <- as.data.frame(table(sce_neighborTemp$SubType))

                neiborAbunTemp[, 2] <- neiborAbunTemp[, 2] / (length(NeiborTemp) / 10)
                neiborAbunTemp[, 3] <- rep(ROI, nrow(neiborAbunTemp))

                AInterStrenDF <- rbind(AInterStrenDF, neiborAbunTemp)
            }


            ## Normal Treg
            NInterStrenDF <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))

            for (ROI in names(NTreg_k_NeighborsList)) {
                sceTemp <- sce[, sce$ID == ROI]
                NeiborTemp <- NTreg_k_NeighborsList[[ROI]]

                sce_neighborTemp <- sce[, match(NeiborTemp, sce$CellID)]
                neiborAbunTemp <- as.data.frame(table(sce_neighborTemp$SubType))

                neiborAbunTemp[, 2] <- neiborAbunTemp[, 2] / (length(NeiborTemp) / 10)
                neiborAbunTemp[, 3] <- rep(ROI, nrow(neiborAbunTemp))

                NInterStrenDF <- rbind(NInterStrenDF, neiborAbunTemp)
            }

            ANInterStrenDF <- rbind(AInterStrenDF, NInterStrenDF)
            colnames(ANInterStrenDF) <- c("NeighborType", "InteractionStrenth", "ROI")

            ANInterStrenDF$Group <- c(rep("ATreg", nrow(AInterStrenDF)), rep("NTreg", nrow(NInterStrenDF)))
        }

        ## Visualization
        if (T) {
            interstType <- c("Macro_CD163", "Mono_CD11c", "Treg", "NK", "CD8T", "B", "CD4T")
            interstType <- names(table(sce$SubType))

            plotdf <- ANInterStrenDF[ANInterStrenDF$NeighborType %in% interstType, ]

            p <- ggplot(plotdf, aes(x = NeighborType, y = InteractionStrenth, fill = Group)) +
                geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
                # geom_jitter(width = 0.5, size = 0.1, alpha = 0.2) +
                scale_y_continuous(name = "Interaction Strength", limits = c(0, 3)) +
                scale_x_discrete(name = "Neighbor Type") +
                scale_fill_jama() +
                stat_compare_means(aes(group = Group),
                    method = "t.test",
                    hide.ns = TRUE, # Hide non-significant comparisons
                    label = "p.signif",
                    label.y.npc = "top"
                ) +
                theme_bw() +
                theme(
                    plot.title = element_text(size = 14, face = "bold"),
                    text = element_text(size = 12),
                    axis.title = element_text(face = "bold"),
                    axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5),
                    strip.background = element_blank(), # Remove facet grid background
                    panel.background = element_blank(), # Make the panel background empty
                    panel.grid.major = element_blank(), # Remove major grid lines
                    panel.grid.minor = element_blank(), # Remove minor grid lines
                    panel.border = element_blank() # Remove borders around the facets
                )


            pdf(paste0(savePathTemp1, "Interaction strength difference between Active Treg and other Treg.pdf"))
            print(p)
            dev.off()
        }
    }

    ## Active Treg and Normal Treg neighbors expression comparison
    if (T) {
        neighbors_ <- unique(neighbors)

        ## Neigobrs SCE
        sce_neighbor <- sce[, match(neighbors_, sce$CellID)]
        sce_neighbor <- sce_neighbor[, sce_neighbor$Tissue == "IM"]
        table(sce_neighbor$SubType)

        PD1_Treg_associated_cell <- sce_neighbor$CellID

        ## Assign label
        sce$PD1Treg_Asso <- FALSE
        sce[, match(PD1_Treg_associated_cell, sce$CellID)]$PD1Treg_Asso <- TRUE

        ## Treg associated cells differential express analysis
        metaMarkers <- c(
            "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
            "VEGF", "CAIX" ## Hypoxia
        )

        immuneMarkers <- c(
            "CD279", "TIGIT", "CD366", ## Immuno-checkpoint
            "CD27" ## Immuno-activation
        )

        noimmuneMarkers <- c(
            "CD274", "CD80"
        )

        DEGsOfNeiborTypeDFImmune <- SpeciNeiDEGAnalysis(sce, IDcol = "PD1Treg_Asso", analysisType = c("B", "CD4T", "CD8T", "NK", "Treg"), metaMarkers = c(metaMarkers, immuneMarkers), sample.size = 1e4)
        DEGsOfNeiborTypeDFMye <- SpeciNeiDEGAnalysis(sce, IDcol = "PD1Treg_Asso", analysisType = c("Macro_CD163", "Macro_CD169", "Macro_HLADR", "Mono_CD11c", "Mono_Classic", "Mono_Intermediate"), metaMarkers = c(metaMarkers, noimmuneMarkers), sample.size = 1e4)
        DEGsOfNeiborTypeDF <- SpeciNeiDEGAnalysis(sce, IDcol = "PD1Treg_Asso", analysisType = c("SC_FAP", "SC_COLLAGEN", "SC_aSMA", "SC_Vimentin", "TC_CAIX", "TC_EpCAM", "TC_Ki67", "TC_VEGF"), metaMarkers = c(metaMarkers, noimmuneMarkers), sample.size = 1e4)

        ## DEGs visualization
        if (T) {
            visuaDF <- rbind(DEGsOfNeiborTypeDFImmune, DEGsOfNeiborTypeDFMye, DEGsOfNeiborTypeDF)

            colnames(visuaDF) <- c("NeighborType", "Marker", "FC", "P.value", "Q.value")
            visuaDF$MajorType <- colData(sce)[match(visuaDF$NeighborType, sce$SubType), "MajorType"]
            visuaDF$MajorType <- ifelse(visuaDF$MajorType %in% c("Lymphocyte"), "Immune", ifelse(visuaDF$MajorType %in% "Myeloid", "Myeloid", "Non-Immune"))

            visuaDF$FC <- as.numeric(visuaDF$FC)
            visuaDF$P.value <- as.numeric(visuaDF$P.value)
            visuaDF$Q.value <- as.numeric(visuaDF$Q.value)

            visuaDF$dir <- ""
            FC_threshold <- 2
            Q_threshold <- 0.05

            for (i in 1:nrow(visuaDF)) {
                if ((visuaDF[i, "FC"] >= FC_threshold) & (visuaDF[i, "Q.value"] <= Q_threshold)) {
                    visuaDF$dir[i] <- "up-regulated"
                }
                if ((visuaDF[i, "FC"] <= (1 / FC_threshold)) & (visuaDF[i, "Q.value"] <= Q_threshold)) {
                    visuaDF$dir[i] <- "down-regulated"
                }
            }

            visuaDF$label <- ifelse(visuaDF$dir != "", visuaDF$Marker, "")

            mycol <- c(ggsci::pal_npg("nrc", alpha = 0.8)(2), "gray")
            names(mycol) <- c("up-regulated", "down-regulated", "NOT")

            visuaDF$log2FC <- log2(visuaDF$FC)

            p1 <- ggplot(visuaDF, aes(x = NeighborType, y = log2FC, color = NeighborType)) +
                geom_jitter(aes(x = NeighborType, y = log2FC, color = dir), size = 0.2, width = 0.2) +
                theme_classic() +
                geom_text_repel(aes(label = label), size = 3) +
                scale_color_manual(values = mycol) +
                ylab("log2FC") +
                facet_grid(~MajorType, scales = "free") +
                theme(
                    legend.position = "none",
                    panel.border = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(size = 0.5, colour = "black"),
                    axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                    axis.text.y = element_text(colour = "black"),
                    axis.ticks = element_line(colour = "black"),
                    axis.title.y = element_text(size = 10),
                    plot.title = element_text(size = 10, hjust = 0.5)
                )

            pdf(paste0(savePathTemp1, "CD27+ PD-1+ Treg and other Treg neighrbos marker expression difference.pdf"), width = 10, height = 5)
            print(p1)
            dev.off()
        }
    }


    ## Are CD27, PD-1 expression of Treg correlated with RFS ? (No)
    if (F) {
        sceTreg <- sce[, sce$SubType %in% c("PD-1+ CD27+ Treg", "Treg")]
        sceTreg <- sceTreg[, sceTreg$Tissue == "IM"]
        TregExp <- as.data.frame(t(assay(sceTreg)))


        TregExp <- TregExp[, match(c("CD45", "CD4", "FoxP3", metaMarkers), colnames(TregExp))]
        dim(TregExp)

        p <- quickcor(TregExp, cluster = TRUE, type = "upper", cor.test = TRUE) +
            geom_colour(data = get_data(type = "upper")) +
            geom_mark(data = get_data(type = "upper"), size = 2.5, color = "black", fontface = 1) +
            scale_fill_gradientn(colours = c("#77C034", "white", "#C388FE")) +
            geom_panel_grid(colour = "white", size = 0.6)

        pdf(paste0(savePathTemp1, "Treg marker expression correlation.pdf"), width = 8, height = 8)
        print(p)
        dev.off()

        TregExp$PID <- sceTreg$PID

        ## Take mean
        TregExp <- TregExp %>%
            group_by(PID) %>%
            summarise(across(c(1:(ncol(TregExp) - 1)), mean, na.rm = TRUE))

        MeanDensity <- left_join(TregExp, clinical[, match(c("PID", "RFS_time", "RFS_status"), colnames(clinical))], by = "PID")

        cor(MeanDensity$CD279, MeanDensity$CD27)

        colnames(MeanDensity)

        cutpoint <- surv_cutpoint(data = MeanDensity, time = "RFS_time", event = "RFS_status", variables = "CD27")
        cutpoint <- summary(cutpoint)$cutpoint

        df <- MeanDensity[, c("CD27", "RFS_status", "RFS_time")]
        df[, 1] <- ifelse(df[, 1] >= cutpoint, "high", "low")

        df$RFS_status <- as.numeric(df$RFS_status)
        df$RFS_time <- as.numeric(df$RFS_time)

        ## km curve
        fit <- survfit(Surv(RFS_time, RFS_status) ~ CD27, data = df)
        p <- ggsurvplot(fit,
            data = df,
            linetype = c("solid", "solid"),
            surv.median.line = "hv", surv.scale = "percent",
            pval = T, risk.table = T,
            conf.int = T, conf.int.alpha = 0.1, conf.int.style = "ribbon",
            risk.table.y.text = T,
            palette = c("#3300CC", "#CC3300"),
            xlab = "Recurrence time"
        )

        pdf("KM for Treg CD27.pdf", width = 8, height = 6)
        print(p)
        dev.off()
    }
}

## Grouping IM ROIs into different abundance group
if (T) {
    sce
    structure <- "CNP20"

    sce_ <- readRDS(paste0(savePathTemp1, "Treg_recluster.rds"))
    sce_ <- sce_[, sce_$Tissue == "IM"]

    ## Assign the new label
    interstType <- c("Pheno_pos")

    TempList <- AssignNewlabel(
        sce_ = sce_,
        allROIs = names(table(sce_$ID)), phenoLabel = "Treg_Type2", ReclusterName = "Treg_Type", interstType = interstType,
        clinicalFeatures = c("RFS_time", "RFS_status"), cutoffType = "manual", cutoffValue = 10, numGroup = 2
    )
    sce_Temp <- TempList[[1]]
    high_ROI <- TempList[[2]]
    low_ROI <- TempList[[3]]

    # colData(sce)[match(colData(sce_Temp)[colData(sce_Temp)$Treg_Type=="Pheno_pos", ]$CellID, colData(sce)$CellID), "SubType"] <- "PD1+ Treg"
    # table(sce$SubType)
    ## saveRDS(sce,"/mnt/data/lyx/IMC/analysis/allsce(new).rds")

    ## spatial difference
    ResultPath <- paste0("/mnt/data/lyx/IMC/analysis/spatial/permutation_IM/")
    celltypes <- names(table(sce$SubType))
    celltypes <- celltypes[-12] ## remove "PD-1+ CD27+ Treg"

    ### high ROIs
    list_ <- getResult(ResultPath, ROIs = high_ROI, celltypes)
    MergeDF1 <- list_[[1]]
    labelDF1 <- list_[[2]]
    rm(list_)

    ### low ROIs
    list_ <- getResult(ResultPath, ROIs = low_ROI, celltypes)
    MergeDF2 <- list_[[1]]
    labelDF2 <- list_[[2]]
    rm(list_)

    DoubleHeat(MergeDF1, labelDF1,
        group1 = paste0("Active Treg high"),
        MergeDF2, labelDF2, group2 = paste0("Active Treg low"), plot = "groupheatmap",
        savePath = paste0(savePathTemp1, "Active Treg Interaction Heatmap.pdf")
    )

    ### plot the cell subtype number different
    sce_Temp1 <- sce[, sce$ID %in% high_ROI]
    HighCountMat1 <- GetAbundance(sce_Temp1, countcol = "SubType", clinicalFeatures = "RFS_status", is.fraction = T)
    HighCountMat2 <- GetAbundance(sce_Temp1, countcol = structure, clinicalFeatures = "RFS_status", is.fraction = T)

    sce_Temp2 <- sce[, sce$ID %in% low_ROI]
    lowCountMat1 <- GetAbundance(sce_Temp2, countcol = "SubType", clinicalFeatures = "RFS_status", is.fraction = T)
    lowCountMat2 <- GetAbundance(sce_Temp2, countcol = structure, clinicalFeatures = "RFS_status", is.fraction = T)


    #### celltype
    savePathTemp <- paste0(savePathTemp1, "CellAbundance/")
    if (!file.exists(savePathTemp)) {
        dir.create(savePathTemp, recursive = T)
    }
    for (celltype in celltypes) {
        AbundanceSwarmPlot(HighCountMat1, lowCountMat1, groupsName = c("high", "low"), celltype = celltype, marker = "Acitve Treg", savePath = savePathTemp, numGroup = 2)
    }

    #### Multiple cell subpopulations abudnace plot
    PlotTypes <- c("CD8T", "Macro_HLADR", "Macro_CD11b", "Macro_CD163", "Macro_CD169", "Treg")
    p <- MultipleAbundanceSwarmPlot(HighCountMat1, lowCountMat1, groupsName = c("high", "low"), PlotTypes = PlotTypes, marker = "Acitve Treg", numGroup = 2, style = "box")
    pdf(paste0(savePathTemp1, "ViolinPlot of multiple celltypes in ", "Acitve Treg", " group.pdf"), height = 10, width = 8)
    print(p)
    dev.off()

    #### CNPs
    k_structures <- names(table(colData(sce)[, structure]))
    savePathTemp <- paste0(savePathTemp1, "StructureAbundance/")
    if (!file.exists(savePathTemp)) {
        dir.create(savePathTemp, recursive = T)
    }
    for (k_structure in k_structures) {
        AbundanceSwarmPlot(HighCountMat2, lowCountMat2, groupsName = c("high", "low"), celltype = k_structure, marker = "Acitve Treg", savePath = savePathTemp, numGroup = 2)
    }

    #### survival
    CountMat <- GetAbundance(sce_Temp, countcol = paste0("Treg", "_Type"), clinicalFeatures = c("RFS_time", "RFS_status"), is.fraction = T, is.reuturnMeans = T)
    SurvivalForPhenoAssoLabel(CountMat, GroupCol = "Pheno_pos", time = "RFS_time", status = "RFS_status", marker = "Acitve Treg", cutoffType = "best", savePath = savePathTemp1)

    #### Save PD-1+ Treg abudance in IM to construct finaly model
    if (T) {
        CountMat <- GetAbundance(sce_Temp, countcol = paste0(SubName, "_Type"), clinicalFeatures = c("RFS_time", "RFS_status"), is.fraction = T, is.reuturnMeans = F)
        CountMat <- CountMat[, -1]
        write.table(CountMat, "PD-1+ Treg Abundance for model construction.csv", sep = ",", row.names = T, col.names = T)
    }
}


if (T) {
    sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
    sce <- sce[, sce$Tissue == "IM"]
    sce

    sce_ <- readRDS(paste0(savePathTemp1, "Treg_recluster.rds"))
    sce_ <- sce_[, sce_$Tissue == "IM"]
    sce_

    head(colData(sce_))
    head(colData(sce))

    colData(sce)[match(sce_[, sce_$Treg_Type == "Pheno_pos"]$CellID, sce$CellID), "SubType"] <- "Acitve Treg"

    ## The k nearest neighbours of Cells in TAT
    Minor_cell_neighbors_knn10 <- read.csv("/mnt/data/lyx/IMC/analysis/spatial/SubType_spatial_count_knn_10_IM.csv")
    Minor_cell_neighbors_knn10 <- Minor_cell_neighbors_knn10[, -1]
    head(Minor_cell_neighbors_knn10)

    Minor_cell_neighbors_knn10 <- Minor_cell_neighbors_knn10[, -21] ## Remove Unknown

    Minor_cell_neighbors_knn10 <- cbind(Minor_cell_neighbors_knn10, sce$ID)

    Minor_cell_neighbors_knn10[, 1:(ncol(Minor_cell_neighbors_knn10) - 1)] <- floor(Minor_cell_neighbors_knn10[, 1:(ncol(Minor_cell_neighbors_knn10) - 1)] * 10)
    colnames(Minor_cell_neighbors_knn10)[ncol(Minor_cell_neighbors_knn10)] <- c("ID")

    sourceTypeList <- c("Treg", "Acitve Treg")

    mat <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 22))

    for (sourceType in sourceTypeList) {
        idx <- (sce$SubType == sourceType)

        Minor_cell_neighbors_knn10Temp <- Minor_cell_neighbors_knn10[idx, ]
        Minor_cell_neighbors_knn10Temp$Label <- rep(sourceType, nrow(Minor_cell_neighbors_knn10Temp))
        mat <- rbind(mat, Minor_cell_neighbors_knn10Temp)
    }
    ## remove Treg subtypes
    mat <- mat[, -c(20)]

    xCol <- c(1, 19)
    yCol <- 21
    mat_foldchangeMat <- FCandPvalueCal(mat, xCol = xCol, yCol = yCol, groups = c("Treg", "Acitve Treg"))

    VolcanoPlot(mat_foldchangeMat, pthreshold = 0.05, fcthreshold = 1.4, feature = NULL, filename = "Volcano of neigobors in PD1+ Treg and Other Treg.pdf")
}

## PD-L1+ DC
## Mono_CD11c analysis - 2
## Hypothesis: The PD-1+ Mono_CD11c only in IM but not in TAT
## Mono_CD11c Reclustering
if (T) {
    TargeType <- "Mono_CD11c"
    savePathTemp1 <- paste0(savePath, TargeType, "/")
    if (!dir.exists(savePathTemp1)) {
        dir.create(savePathTemp1, recursive = T)
    }

    metaMarkers <- c(
        "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
        "VEGF", "CAIX", ## Hypoxia
        "CD274", "CD80" ## Immuno-checkpoint
    )

    sce <- sce[, sce$Tissue %in% c("IM", "CT", "TAT")]
    sce_Type <- sce[, sce$SubType %in% TargeType]

    ## Dimention Redcution
    m <- DimenRec(sce_ = sce_Type, GroupFeature = "RFS_status", markers = metaMarkers, do.pca = FALSE, pca.dimen = 10, explaine.pca = F, embeeding.method = "tsne", num.threads = 8)

    ## calculate Fraction and Entropy via KNN graph
    k <- 20
    cordim <- 2

    resultdf <- KGandFracEntro(m, cordim = cordim, k = k, central.weight = 2, multi.cores = 8)

    ## Generate Fraction Null Distribution
    fraction_NullDistribution <- GeneNullDistri(m, cordim = cordim, sample.size = k, perm.times = 2000)

    ## Quantile of null distribution
    threshold <- 0.05
    low.95 <- as.numeric(quantile(fraction_NullDistribution, threshold))
    high.95 <- as.numeric(quantile(fraction_NullDistribution, (1 - threshold)))

    trueFraction <- resultdf[, 2]
    Class <- ifelse(trueFraction <= low.95, "Pheno_neg", ifelse(
        trueFraction >= high.95, "Pheno_pos", "Background"
    ))

    table(Class)

    ## Visualize class label and entropy
    PlotClassLabel(
        sce_ = sce_Type, Classlabels = Class, Entropylabels = resultdf[, 1], markers = metaMarkers, sample.size = 1e4, num.threads = 8,
        SavePath1 = paste0(savePathTemp1, "T-SNE of phenotype associated label.pdf"), SavePath2 = paste0(savePathTemp1, "T-SNE of phenotype associated entropy.pdf"),
        SavePath3 = savePathTemp1, seed = 619
    )

    ## Volcano for Pheno_pos and Pheno_neg cells

    ### Relapse associated
    mat <- as.data.frame(t(assay(sce_Type)))
    mat <- mat[, match(metaMarkers, colnames(mat))]
    mat$label <- Class

    # mat <- mat[!(mat$label%in%"Background"),]
    mat$label <- ifelse(mat$label == "Pheno_pos", 1, 0)

    mat_foldchangeMat <- FCandPvalueCal(mat = mat, xCol = c(1, length(metaMarkers)), yCol = (length(metaMarkers) + 1), need.sample = T, sample.size = 5e3)
    mat_foldchangeMat$Q.value <- p.adjust(mat_foldchangeMat$P.value, method = "BH")
    mat_foldchangeMat$P.value <- mat_foldchangeMat$Q.value

    VolcanoPlot(mat_foldchangeMat, pthreshold = 0.05, fcthreshold = 1.5, feature = NULL, filename = paste0(savePathTemp1, "DEGs of relapse associated Mono_CD11cs.pdf"))

    ### Non relapse associated
    mat <- as.data.frame(t(assay(sce_Type)))
    mat <- mat[, match(metaMarkers, colnames(mat))]
    mat$label <- Class

    mat$label <- ifelse(mat$label == "Pheno_neg", 1, 0)

    mat_foldchangeMat <- FCandPvalueCal(mat = mat, xCol = c(1, length(metaMarkers)), yCol = (length(metaMarkers) + 1), need.sample = T, sample.size = 5e3)
    mat_foldchangeMat$Q.value <- p.adjust(mat_foldchangeMat$P.value, method = "BH")
    mat_foldchangeMat$P.value <- mat_foldchangeMat$Q.value

    VolcanoPlot(mat_foldchangeMat, pthreshold = 0.05, fcthreshold = 1.5, feature = NULL, filename = paste0(savePathTemp1, "DEGs of non-relapse associated Mono_CD11cs.pdf"))

    sce_Type$Mono_CD11c_Type <- Class

    saveRDS(sce_Type, paste0(savePathTemp1, "Mono_CD11c_recluster.rds"))
}

## investigate the tissue specificity of PD-1+ CD27+ Mono_CD11c
if (T) {
    sce_ <- readRDS(paste0(savePathTemp1, "Mono_CD11c_recluster.rds"))

    ## Cluster difference in Tissue
    plotdf2 <- GetAbundance(sce_, countcol = "Mono_CD11c_Type", clinicalFeatures = c("Tissue", "RFS_status"), is.fraction = T)
    plotdf2 <- plotdf2[, c("Pheno_pos", "PID", "Tissue", "RFS_status")]

    if (F) {
        write.table(plotdf2, "PD-L1+ DC Abundance for model construction.csv", sep = ",", row.names = T, col.names = T)
    }

    plotdf2 <- plotdf2[plotdf2$Tissue %in% c("IM", "CT", "TAT"), ]
    colnames(plotdf2)[1] <- "Abundance"

    p <- ggplot(plotdf2, aes(x = Tissue, y = Abundance, fill = "#56B4E9")) +
        geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
        scale_y_continuous(name = "Cell Abundance", limits = c(0, 0.6)) +
        scale_x_discrete(name = "Cell Population") +
        theme_bw() +
        theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            text = element_text(size = 12),
            axis.title = element_text(face = "bold"),
            axis.text.x = element_text(size = 11, angle = 90),
            legend.position = "none",
            strip.text = element_text(size = 12, face = "bold"),
            strip.background = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            panel.grid.major = element_line(color = "gray", size = 0.5),
            panel.grid.minor = element_blank()
        )

    stat.test <- plotdf2 %>%
        t_test(Abundance ~ Tissue) %>%
        add_significance() %>%
        add_xy_position(fun = "mean_sd", scales = "free_y")

    p1 <- p + stat_pvalue_manual(stat.test, hide.ns = FALSE, label = "{p}", step.increase = 0.08)

    pdf(paste0(savePathTemp1, "Pheno-associated Mono_CD11c tissue distribution.pdf"), width = 4, height = 3)
    print(p1)
    dev.off()
}

## Plot the location of PD-1+ CD27+ Mono_CD11c
if (T) {
    sce <- sce[, sce$Tissue == "IM"]
    sce_ <- readRDS(paste0(savePathTemp1, "Mono_CD11c_recluster.rds"))

    head(colData(sce_))
    sce__ <- sce_[, sce_$Mono_CD11c_Type == "Pheno_pos"]
    sce__ <- sce__[, sce__$Tissue == "IM"]

    sort(table(sce__$ID), decreasing = T)

    source("./spatial_analysis_functions.r")

    PD1Mono_CD11cID <- sce__$CellID
    sce[, match(PD1Mono_CD11cID, sce$CellID)]$SubType <- "CD27CD279+ Mono_CD11c"

    ROI_num <- sort(table(sce__$ID), decreasing = T)
    ROIs <- names(ROI_num[ROI_num > 20])

    celltype <- c("CD27CD279+ Mono_CD11c")
    for (ROI in ROIs) {
        VisualizeMarkerAndType(sce = sce, ROI = ROI, celltype = celltype, whichType = "SubType", marker = NULL, SavePath = savePathTemp1)
    }
}

## investigate the spatial interaction of PD-1+ CD27+ Mono_CD11c
if (T) {
    library(SingleCellExperiment)
    library(spatstat)
    source("./structural_analysis_functions.r")

    TargeType <- "Mono_CD11c"
    savePathTemp1 <- paste0(savePath, TargeType, "/")

    sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
    sce <- sce[, sce$Tissue %in% c("IM", "TAT")]

    sce_ <- readRDS(paste0("/mnt/data/lyx/IMC/analysis/reclustering/Mono_CD11c/Mono_CD11c_recluster.rds"))
    sce_ <- sce_[, sce_$Tissue %in% c("IM", "TAT")]
    sce_ <- sce_[, sce_$Mono_CD11c_Type %in% c("Pheno_pos")]

    sce[, match(sce_$CellID, sce$CellID)]$SubType <- "PD-L1+ Mono_CD11c"
    table(sce$SubType)

    ## spatial interaction analysis
    ROIs <- names(table(sce$ID))

    # Define the whole image size
    k_NeighborsList <- CalKNNBySpatstat(sce, ROIs = ROIs, CenType = "PD-L1+ Mono_CD11c", NeighType = "SubType", k = 10, xlim = c(0, 1000), ylim = c(0, 1000), return.cellID = TRUE)

    ## Combine all neighbors
    neighbors <- c()

    for (i in k_NeighborsList) {
        neighbors <- c(neighbors, i)
    }

    neighbors <- neighbors[(neighbors %in% sce$CellID)]

    ## Calculate the neighbors type abudance
    sce_neighbor <- sce[, match(neighbors, sce$CellID)]
    neiborAbun <- as.data.frame(table(sce_neighbor$SubType))

    ### Visualize neighbors abudance
    if (T) {
        # Rename the columns
        colnames(neiborAbun) <- c("category", "counts")

        # Calculate the percentage for each category
        neiborAbun$percentage <- neiborAbun$counts / sum(neiborAbun$counts)

        # Define labels for the pie chart
        neiborAbun$label <- paste0(neiborAbun$category, "\n(", round(neiborAbun$percentage * 100, 2), "%)")

        # Sort the dataframe by percentage (in descending order)
        neiborAbun <- neiborAbun[order(-neiborAbun$percentage), ]

        # Select top k categories for labelling on the chart
        k <- 7
        neiborAbun$label[(k + 1):nrow(neiborAbun)] <- ""

        # Drop unused factor levels in the 'category' column
        neiborAbun$category <- droplevels(neiborAbun$category)

        # Create the pie chart
        p <- ggplot(neiborAbun, aes(x = "", y = percentage, fill = category)) +
            geom_bar(width = 1, stat = "identity", color = "black") +
            geom_label_repel(aes(label = ifelse(percentage > 0.08, label, "")), # labels for slices < 1% are set to ""
                nudge_y = 0.1, size = 4, show.legend = FALSE
            ) +
            coord_polar("y", start = 0) +
            theme_bw() +
            scale_fill_manual(values = ggsci::pal_ucscgb("default")(length(levels(neiborAbun$category)))) +
            theme(
                panel.grid.major = element_line(color = "gray", size = 0.2),
                panel.grid.minor = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                plot.title = element_text(size = 18, hjust = 0.5),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 12),
                legend.position = "bottom"
            )


        pdf(paste0(savePathTemp1, "PD-L1+ Mono_CD11c neighbors abundance.pdf"), width = 10, height = 6)
        print(p)
        dev.off()
    }

    ## Calcualte the interaction strength
    if (T) {
        InterStrenDF <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))

        for (ROI in names(k_NeighborsList)) {
            sceTemp <- sce[, sce$ID == ROI]
            NeiborTemp <- k_NeighborsList[[ROI]]

            sce_neighborTemp <- sce[, match(NeiborTemp, sce$CellID)]
            neiborAbunTemp <- as.data.frame(table(sce_neighborTemp$SubType))

            neiborAbunTemp[, 2] <- neiborAbunTemp[, 2] / (length(NeiborTemp) / 10)
            neiborAbunTemp[, 3] <- rep(ROI, nrow(neiborAbunTemp))

            InterStrenDF <- rbind(InterStrenDF, neiborAbunTemp)
        }

        colnames(InterStrenDF) <- c("NeighborType", "InteractionStrenth", "ROI")
        InterStrenDF$Tissue <- colData(sce)[match(InterStrenDF$ROI, sce$ID), "Tissue"]

        ## Visualization
        if (T) {
            interstType <- c("Macro_CD163", "Mono_CD11c", "Mono_CD11c", "NK", "CD8T", "B", "CD4T")

            plotdf <- InterStrenDF[InterStrenDF$NeighborType %in% interstType, ]

            p <- ggplot(plotdf, aes(x = NeighborType, y = InteractionStrenth, fill = Tissue)) +
                geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
                scale_y_continuous(name = "Interaction Strength", limits = c(0, 3.5)) +
                scale_x_discrete(name = "Neighbor Type") +
                scale_fill_lancet() +
                stat_compare_means(aes(group = Tissue),
                    method = "t.test",
                    hide.ns = TRUE, # Hide non-significant comparisons
                    label = "p.signif",
                    label.y.npc = "top"
                ) +
                theme_bw() +
                theme(
                    plot.title = element_text(size = 14, face = "bold"),
                    text = element_text(size = 12),
                    axis.title = element_text(face = "bold"),
                    axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5),
                    strip.background = element_blank(), # Remove facet grid background
                    panel.background = element_blank(), # Make the panel background empty
                    panel.grid.major = element_blank(), # Remove major grid lines
                    panel.grid.minor = element_blank(), # Remove minor grid lines
                    panel.border = element_blank() # Remove borders around the facets
                )


            pdf(paste0(savePathTemp1, "PD-L1+ Mono_CD11c interaction strength between tissue.pdf"))
            print(p)
            dev.off()
        }
    }

    ## PDL1+ Mono_CD11c and Normal Mono_CD11c neighbors fraction comparison
    if (T) {
        AMono_CD11c_k_NeighborsList <- k_NeighborsList
        NMono_CD11c_k_NeighborsList <- CalKNNBySpatstat(sce, ROIs = ROIs, CenType = "Mono_CD11c", NeighType = "SubType", k = 10, xlim = c(0, 1000), ylim = c(0, 1000), return.cellID = TRUE)

        ## Calcualte the interaction strength
        if (T) {
            ## PDL1+ Mono_CD11c
            AInterStrenDF <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))

            for (ROI in names(AMono_CD11c_k_NeighborsList)) {
                sceTemp <- sce[, sce$ID == ROI]
                NeiborTemp <- AMono_CD11c_k_NeighborsList[[ROI]]

                sce_neighborTemp <- sce[, match(NeiborTemp, sce$CellID)]
                neiborAbunTemp <- as.data.frame(table(sce_neighborTemp$SubType))

                neiborAbunTemp[, 2] <- neiborAbunTemp[, 2] / (length(NeiborTemp) / 10)
                neiborAbunTemp[, 3] <- rep(ROI, nrow(neiborAbunTemp))

                AInterStrenDF <- rbind(AInterStrenDF, neiborAbunTemp)
            }


            ## Normal Mono_CD11c
            NInterStrenDF <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))

            for (ROI in names(NMono_CD11c_k_NeighborsList)) {
                sceTemp <- sce[, sce$ID == ROI]
                NeiborTemp <- NMono_CD11c_k_NeighborsList[[ROI]]

                sce_neighborTemp <- sce[, match(NeiborTemp, sce$CellID)]
                neiborAbunTemp <- as.data.frame(table(sce_neighborTemp$SubType))

                neiborAbunTemp[, 2] <- neiborAbunTemp[, 2] / (length(NeiborTemp) / 10)
                neiborAbunTemp[, 3] <- rep(ROI, nrow(neiborAbunTemp))

                NInterStrenDF <- rbind(NInterStrenDF, neiborAbunTemp)
            }

            ANInterStrenDF <- rbind(AInterStrenDF, NInterStrenDF)
            colnames(ANInterStrenDF) <- c("NeighborType", "InteractionStrenth", "ROI")

            ANInterStrenDF$Group <- c(rep("AMono_CD11c", nrow(AInterStrenDF)), rep("NMono_CD11c", nrow(NInterStrenDF)))
        }

        ## Visualization
        if (T) {
            interstType <- c("Macro_CD163", "Mono_CD11c", "Mono_CD11c", "NK", "CD8T", "B", "CD4T")
            interstType <- names(table(sce$SubType))

            plotdf <- ANInterStrenDF[ANInterStrenDF$NeighborType %in% interstType, ]

            p <- ggplot(plotdf, aes(x = NeighborType, y = InteractionStrenth, fill = Group)) +
                geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
                # geom_jitter(width = 0.5, size = 0.1, alpha = 0.2) +
                scale_y_continuous(name = "Interaction Strength", limits = c(0, 3)) +
                scale_x_discrete(name = "Neighbor Type") +
                scale_fill_jama() +
                stat_compare_means(aes(group = Group),
                    method = "t.test",
                    hide.ns = TRUE, # Hide non-significant comparisons
                    label = "p.signif",
                    label.y.npc = "top"
                ) +
                theme_bw() +
                theme(
                    plot.title = element_text(size = 14, face = "bold"),
                    text = element_text(size = 12),
                    axis.title = element_text(face = "bold"),
                    axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5),
                    strip.background = element_blank(), # Remove facet grid background
                    panel.background = element_blank(), # Make the panel background empty
                    panel.grid.major = element_blank(), # Remove major grid lines
                    panel.grid.minor = element_blank(), # Remove minor grid lines
                    panel.border = element_blank() # Remove borders around the facets
                )


            pdf(paste0(savePathTemp1, "Interaction strength difference between PDL1+ Mono_CD11c and other Mono_CD11c.pdf"))
            print(p)
            dev.off()
        }
    }

    ## PDL1+ Mono_CD11c and Normal Mono_CD11c neighbors expression comparison
    if (T) {
        neighbors_ <- unique(neighbors)

        ## Neigobrs SCE
        sce_neighbor <- sce[, match(neighbors_, sce$CellID)]
        sce_neighbor <- sce_neighbor[, sce_neighbor$Tissue == "IM"]
        table(sce_neighbor$SubType)

        PD1_Mono_CD11c_associated_cell <- sce_neighbor$CellID

        ## Assign label
        sce$PDL1Mono_CD11c_Asso <- FALSE
        sce[, match(PD1_Mono_CD11c_associated_cell, sce$CellID)]$PDL1Mono_CD11c_Asso <- TRUE

        ## Mono_CD11c associated cells differential express analysis
        metaMarkers <- c(
            "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
            "VEGF", "CAIX" ## Hypoxia
        )

        immuneMarkers <- c(
            "CD279", "TIGIT", "CD366", ## Immuno-checkpoint
            "CD27" ## Immuno-activation
        )

        noimmuneMarkers <- c(
            "CD274", "CD80"
        )

        DEGsOfNeiborTypeDFImmune <- SpeciNeiDEGAnalysis(sce, IDcol = "PDL1Mono_CD11c_Asso", analysisType = c("B", "CD4T", "CD8T", "NK", "Treg"), metaMarkers = c(metaMarkers, immuneMarkers), sample.size = 1e4)
        DEGsOfNeiborTypeDFMye <- SpeciNeiDEGAnalysis(sce, IDcol = "PDL1Mono_CD11c_Asso", analysisType = c("Macro_CD163", "Macro_CD169", "Macro_HLADR", "Mono_CD11c", "Mono_Classic", "Mono_Intermediate"), metaMarkers = c(metaMarkers, noimmuneMarkers), sample.size = 1e4)
        DEGsOfNeiborTypeDF <- SpeciNeiDEGAnalysis(sce, IDcol = "PDL1Mono_CD11c_Asso", analysisType = c("SC_FAP", "SC_COLLAGEN", "SC_aSMA", "SC_Vimentin", "TC_CAIX", "TC_EpCAM", "TC_Ki67", "TC_VEGF"), metaMarkers = c(metaMarkers, noimmuneMarkers), sample.size = 1e4)

        ## DEGs visualization
        if (T) {
            visuaDF <- rbind(DEGsOfNeiborTypeDFImmune, DEGsOfNeiborTypeDFMye, DEGsOfNeiborTypeDF)

            colnames(visuaDF) <- c("NeighborType", "Marker", "FC", "P.value", "Q.value")
            visuaDF$MajorType <- colData(sce)[match(visuaDF$NeighborType, sce$SubType), "MajorType"]
            visuaDF$MajorType <- ifelse(visuaDF$MajorType %in% c("Lymphocyte"), "Immune", ifelse(visuaDF$MajorType %in% "Myeloid", "Myeloid", "Non-Immune"))

            visuaDF$FC <- as.numeric(visuaDF$FC)
            visuaDF$P.value <- as.numeric(visuaDF$P.value)
            visuaDF$Q.value <- as.numeric(visuaDF$Q.value)

            visuaDF$dir <- ""
            FC_threshold <- 1.5
            Q_threshold <- 0.05

            for (i in 1:nrow(visuaDF)) {
                if ((visuaDF[i, "FC"] >= FC_threshold) & (visuaDF[i, "Q.value"] <= Q_threshold)) {
                    visuaDF$dir[i] <- "up-regulated"
                }
                if ((visuaDF[i, "FC"] <= (1 / FC_threshold)) & (visuaDF[i, "Q.value"] <= Q_threshold)) {
                    visuaDF$dir[i] <- "down-regulated"
                }
            }

            visuaDF$label <- ifelse(visuaDF$dir != "", visuaDF$Marker, "")

            mycol <- c(ggsci::pal_npg("nrc", alpha = 0.8)(2), "gray")
            names(mycol) <- c("up-regulated", "down-regulated", "NOT")

            visuaDF$log2FC <- log2(visuaDF$FC)

            p1 <- ggplot(visuaDF, aes(x = NeighborType, y = log2FC, color = NeighborType)) +
                geom_jitter(aes(x = NeighborType, y = log2FC, color = dir), size = 0.2, width = 0.2) +
                theme_classic() +
                geom_text_repel(aes(label = label), size = 3) +
                scale_color_manual(values = mycol) +
                ylab("log2FC") +
                facet_grid(~MajorType, scales = "free") +
                theme(
                    legend.position = "none",
                    panel.border = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(size = 0.5, colour = "black"),
                    axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                    axis.text.y = element_text(colour = "black"),
                    axis.ticks = element_line(colour = "black"),
                    axis.title.y = element_text(size = 10),
                    plot.title = element_text(size = 10, hjust = 0.5)
                )

            pdf(paste0(savePathTemp1, "PD-L1+ Mono_CD11c and other Mono_CD11c neighrbos marker expression difference.pdf"), width = 10, height = 5)
            print(p1)
            dev.off()
        }
    }


    ## Are CD27, PD-1 expression of Mono_CD11c correlated with RFS ? (No)
    if (F) {
        sceMono_CD11c <- sce[, sce$SubType %in% c("PD-L1+ Mono_CD11c", "Mono_CD11c")]
        sceMono_CD11c <- sceMono_CD11c[, sceMono_CD11c$Tissue == "IM"]
        Mono_CD11cExp <- as.data.frame(t(assay(sceMono_CD11c)))

        Mono_CD11cExp <- Mono_CD11cExp[, match(c("CD11c", metaMarkers, noimmuneMarkers), colnames(Mono_CD11cExp))]
        dim(Mono_CD11cExp)

        Mono_CD11cExp$PID <- sceMono_CD11c$PID

        ## Take mean
        Mono_CD11cExp <- Mono_CD11cExp %>%
            group_by(PID) %>%
            summarise(across(c(1:(ncol(Mono_CD11cExp) - 1)), mean, na.rm = TRUE))

        MeanDensity <- left_join(Mono_CD11cExp, clinical[, match(c("PID", "RFS_time", "RFS_status"), colnames(clinical))], by = "PID")

        colnames(MeanDensity)
        cutpoint <- surv_cutpoint(data = MeanDensity, time = "RFS_time", event = "RFS_status", variables = "CD274")
        cutpoint <- summary(cutpoint)$cutpoint

        df <- MeanDensity[, c("CD274", "RFS_status", "RFS_time")]
        df[, 1] <- ifelse(df[, 1] >= cutpoint, "high", "low")

        df$RFS_status <- as.numeric(df$RFS_status)
        df$RFS_time <- as.numeric(df$RFS_time)

        ## km curve
        fit <- survfit(Surv(RFS_time, RFS_status) ~ CD274, data = df)
        p <- ggsurvplot(fit,
            data = df,
            linetype = c("solid", "solid"),
            surv.median.line = "hv", surv.scale = "percent",
            pval = T, risk.table = T,
            conf.int = T, conf.int.alpha = 0.1, conf.int.style = "ribbon",
            risk.table.y.text = T,
            palette = c("#3300CC", "#CC3300"),
            xlab = "Recurrence time"
        )

        pdf("KM for Mono_CD11c CD274.pdf", width = 8, height = 6)
        print(p)
        dev.off()
    }
}

## Grouping IM ROIs into different abundance group
if (T) {
    sce
    structure <- "CNP20"

    sce_ <- readRDS(paste0(savePathTemp1, "Mono_CD11c_recluster.rds"))
    sce_ <- sce_[, sce_$Tissue == "IM"]

    ## Assign the new label
    interstType <- c("Pheno_pos")

    TempList <- AssignNewlabel(
        sce_ = sce_,
        allROIs = names(table(sce_$ID)), phenoLabel = "Mono_CD11c_Type2", ReclusterName = "Mono_CD11c_Type", interstType = interstType,
        clinicalFeatures = c("RFS_time", "RFS_status"), cutoffType = "manual", cutoffValue = 40, numGroup = 2
    )
    sce_Temp <- TempList[[1]]
    high_ROI <- TempList[[2]]
    low_ROI <- TempList[[3]]

    # colData(sce)[match(colData(sce_Temp)[colData(sce_Temp)$Mono_CD11c_Type=="Pheno_pos", ]$CellID, colData(sce)$CellID), "SubType"] <- "PD1+ Mono_CD11c"
    # table(sce$SubType)
    ## saveRDS(sce,"/mnt/data/lyx/IMC/analysis/allsce(new).rds")

    ## spatial difference
    ResultPath <- paste0("/mnt/data/lyx/IMC/analysis/spatial/permutation_IM/")
    celltypes <- names(table(sce$SubType))
    celltypes <- celltypes[-21] ## remove "PD-L1+ Mono_CD11c"

    ### high ROIs
    list_ <- getResult(ResultPath, ROIs = high_ROI, celltypes)
    MergeDF1 <- list_[[1]]
    labelDF1 <- list_[[2]]
    rm(list_)

    ### low ROIs
    list_ <- getResult(ResultPath, ROIs = low_ROI, celltypes)
    MergeDF2 <- list_[[1]]
    labelDF2 <- list_[[2]]
    rm(list_)

    DoubleHeat(MergeDF1, labelDF1,
        group1 = paste0("PDL1+ Mono_CD11c high"),
        MergeDF2, labelDF2, group2 = paste0("PDL1+ Mono_CD11c low"), plot = "groupheatmap",
        savePath = paste0(savePathTemp1, "PDL1+ Mono_CD11c Interaction Heatmap.pdf")
    )

    ### plot the cell subtype number different
    sce_Temp1 <- sce[, sce$ID %in% high_ROI]
    HighCountMat1 <- GetAbundance(sce_Temp1, countcol = "SubType", clinicalFeatures = "RFS_status", is.fraction = T)
    HighCountMat2 <- GetAbundance(sce_Temp1, countcol = structure, clinicalFeatures = "RFS_status", is.fraction = T)

    sce_Temp2 <- sce[, sce$ID %in% low_ROI]
    lowCountMat1 <- GetAbundance(sce_Temp2, countcol = "SubType", clinicalFeatures = "RFS_status", is.fraction = T)
    lowCountMat2 <- GetAbundance(sce_Temp2, countcol = structure, clinicalFeatures = "RFS_status", is.fraction = T)


    #### celltype
    savePathTemp <- paste0(savePathTemp1, "CellAbundance/")
    if (!file.exists(savePathTemp)) {
        dir.create(savePathTemp, recursive = T)
    }
    for (celltype in celltypes) {
        AbundanceSwarmPlot(HighCountMat1, lowCountMat1, groupsName = c("high", "low"), celltype = celltype, marker = "PDL1+ Mono_CD11c", savePath = savePathTemp, numGroup = 2)
    }

    #### Multiple cell subpopulations abudnace plot
    PlotTypes <- c("CD8T", "Macro_HLADR", "Macro_CD11b", "Macro_CD163", "Macro_CD169", "Mono_CD11c")
    p <- MultipleAbundanceSwarmPlot(HighCountMat1, lowCountMat1, groupsName = c("high", "low"), PlotTypes = PlotTypes, marker = "PDL1+ Mono_CD11c", numGroup = 2, style = "box")
    pdf(paste0(savePathTemp1, "ViolinPlot of multiple celltypes in ", "PDL1+ Mono_CD11c", " group.pdf"), height = 10, width = 8)
    print(p)
    dev.off()

    #### CNPs
    k_structures <- names(table(colData(sce)[, structure]))
    savePathTemp <- paste0(savePathTemp1, "StructureAbundance/")
    if (!file.exists(savePathTemp)) {
        dir.create(savePathTemp, recursive = T)
    }
    for (k_structure in k_structures) {
        AbundanceSwarmPlot(HighCountMat2, lowCountMat2, groupsName = c("high", "low"), celltype = k_structure, marker = "PDL1+ Mono_CD11c", savePath = savePathTemp, numGroup = 2)
    }

    #### survival
    CountMat <- GetAbundance(sce_Temp, countcol = paste0("Mono_CD11c", "_Type"), clinicalFeatures = c("RFS_time", "RFS_status"), is.fraction = T, is.reuturnMeans = T)
    SurvivalForPhenoAssoLabel(CountMat, GroupCol = "Pheno_pos", time = "RFS_time", status = "RFS_status", marker = "PDL1+ Mono_CD11c", cutoffType = "best", savePath = savePathTemp1)

    #### Save PD-1+ Mono_CD11c abudance in IM to construct finaly model
    if (T) {
        CountMat <- GetAbundance(sce_Temp, countcol = paste0(SubName, "_Type"), clinicalFeatures = c("RFS_time", "RFS_status"), is.fraction = T, is.reuturnMeans = F)
        CountMat <- CountMat[, -1]
        write.table(CountMat, "PD-1+ Mono_CD11c Abundance for model construction.csv", sep = ",", row.names = T, col.names = T)
    }
}
