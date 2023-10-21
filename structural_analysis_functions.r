# Functions for structural analysis
library(SingleCellExperiment)
library(parallel)

library(survival)
library(survminer)
library(Hmisc)

library(pheatmap)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(ggrepel)
library(ggcor)
library(ggsci)
library(RColorBrewer)
library(viridis)
library(ggunchained)
library(ggalluvial)

library(stats)
library(Rtsne)
library(FNN)
library(igraph)
library(dplyr)
library(tidyr)

## clinical information
load_clinical <- function(sce, clinicalFilePath) {
    clinical <- read.csv(clinicalFilePath)
    IDs <- names(table(sce$ID))

    ## match the IDs and response level
    GroupInfo <- data.frame(row.names = IDs)

    ### Patient ID
    GroupInfo["PID"] <- sapply(rownames(GroupInfo), function(x) {
        strsplit(x, split = "_")[[1]][1]
    })
    ### RFSS
    GroupInfo["RFS_status"] <- sapply(GroupInfo$PID, function(x) {
        clinical[clinical$PID == x, "RFS_status"]
    })
    GroupInfo["RFS_time"] <- sapply(GroupInfo$PID, function(x) {
        clinical[clinical$PID == x, "RFS_time"]
    })
    ### Mutation
    GroupInfo["KRAS_mutation"] <- sapply(GroupInfo$PID, function(x) {
        clinical[clinical$PID == x, "KRAS_mutation"]
    })
    GroupInfo["Pathology"] <- sapply(GroupInfo$PID, function(x) {
        clinical[clinical$PID == x, "Pathology"]
    })
    GroupInfo["CRC_site"] <- sapply(GroupInfo$PID, function(x) {
        clinical[clinical$PID == x, "CRC_site"]
    })

    return(GroupInfo)
}

## Merge abundance information
MergeAbundanceResult <- function(sce, return.fraction = F) {
    ## celltypes and ROIs
    celltypes <- names(table(sce$SubType))
    ROIs <- names(table(sce$ID))

    AbundanceDF <- matrix(data = 0, nrow = length(celltypes), ncol = length(ROIs))
    AbundanceDF <- as.data.frame(AbundanceDF)
    rownames(AbundanceDF) <- celltypes
    colnames(AbundanceDF) <- ROIs

    for (i in 1:length(ROIs)) {
        ROI <- ROIs[i]
        sceTemp <- sce[, sce$ID == ROI]
        abundanceTemp <- as.data.frame(table(sceTemp$SubType))

        for (j in 1:nrow(abundanceTemp)) {
            rowTemp <- as.character(abundanceTemp[j, 1])
            AbundanceDF[rowTemp, i] <- as.numeric(abundanceTemp[j, 2])
        }
    }
    if (return.fraction) {
        numCell <- apply(AbundanceDF, MARGIN = 2, FUN = "sum")
        numCell <- 1 / numCell
        for (i in 1:ncol(AbundanceDF)) {
            AbundanceDF[, i] <- AbundanceDF[, i] * numCell[i]
        }
    }

    return(AbundanceDF)
}

## Bind the cellular neighborhoods clustering result with sce
BindResult <- function(sce, mat, colName) {
    for (i in colName) {
        colData(sce)[i] <- mat[, i]
    }
    return(sce)
}

## Heatmap for celltypes in cellular neighborhoods
GetCelltype2NeighborMat <- function(mat, colname1, colname2) {
    vec1 <- as.vector(mat[colname1][, 1])
    vec2 <- as.vector(mat[colname2][, 1])

    VecDF <- cbind(vec1, vec2)

    names1 <- names(table(vec1))
    names2 <- names(table(vec2))

    plotdf <- matrix(data = NA, nrow = length(names2), ncol = length(names1))
    rownames(plotdf) <- names2
    colnames(plotdf) <- names1

    for (i in 1:nrow(plotdf)) {
        VecDFTemp <- subset(VecDF, vec2 == rownames(plotdf)[i])
        TableTemp <- as.data.frame(table(VecDFTemp[, "vec1"]))
        rownames(TableTemp) <- TableTemp[, 1]
        TableTemp <- TableTemp[match(colnames(plotdf), rownames(TableTemp)), ]
        TableTemp <- TableTemp[, -1]
        TableTemp <- ifelse(is.na(TableTemp), 0, TableTemp)
        plotdf[i, ] <- TableTemp
    }

    return(as.data.frame(plotdf))
}


HeatmapForCelltypeInNeighbor <- function(sce, colname1, colname2, savePath, scale = "column") {
    ## transfer into heatmap matrix
    plotdf <- GetCelltype2NeighborMat(colData(sce), colname1, colname2)

    ## heatmap
    color <- colorRampPalette(c("#436eee", "white", "#EE0000"))(100)
    p <- pheatmap(plotdf,
        color = color, scale = scale,
        cluster_rows = F, cluster_cols = T,
        # legend_labels = c("Abundance high", "Abundance low"), legend = T,
        show_rownames = T, show_colnames = T
    )

    pdf(paste0(savePath, "Cellular Neighbors celltype fraction heatmap.pdf"), width = 8, height = 6)
    print(p)
    dev.off()
    return(NULL)
}

## comapre and visualize the CNP into tissue
BarPlotForCelltypeCounts <- function(sce, tissueCol, groupCol, typeCol, savePath) {
    countdf <- GetAbundance(sceobj = sce, countcol = typeCol, clinicalFeatures = c("ID", tissueCol, groupCol), is.fraction = T)

    plotdf <- pivot_longer(countdf, cols = 1:10, names_to = "CNP", values_to = "Fraction")

    plotdf <- plotdf[, -c(1:2)]
    plotdf2 <- plotdf %>%
        group_by(Tissue, RFS_status, CNP) %>%
        summarise(across(c(1:(ncol(plotdf) - 3)), mean, na.rm = TRUE))

    # Create a list of colors
    colors <- ggsci::pal_npg("nrc")(10)

    plotdf2 <- as.data.frame(plotdf2)

    plotdf2[, 1] <- as.character(plotdf2[, 1])
    plotdf2[, 2] <- as.factor(plotdf2[, 2])
    plotdf2[, 3] <- as.character(plotdf2[, 3])
    plotdf2[, 4] <- as.numeric(plotdf2[, 4])

    # Create the bar plot
    p <- ggplot(plotdf2, aes(x = RFS_status, y = Fraction, fill = CNP)) +
        # Add bars
        geom_bar(stat = "identity") +
        # Use different sets of colors for different MajorTypes
        scale_fill_manual(values = colors) +
        # Rotate the x-axis labels to make them readable
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        # Add title and labels
        labs(
            title = "Fraction of each CNP by RFS_status and Tissue",
            x = "RFS Status", y = "Fraction", fill = "CNP"
        ) +
        # Separate the plot by Tissue and MajorType
        facet_grid(~Tissue, scales = "free") +
        # Hide the borders of the facet grid, make the background of the facet to be empty
        theme(
            strip.background = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank()
        )

    pdf(paste0(savePath, "Abundance Barplot of CNP.pdf"), height = 4, width = 8)
    print(p)
    dev.off()

    plotdf <- as.data.frame(plotdf)
    plotdf[, 2] <- as.factor(plotdf[, 2])

    plotdf <- plotdf[plotdf$Tissue != "TLS", ]

    ## Boxplot
    p <- ggplot(plotdf, aes(x = CNP, y = Fraction, fill = RFS_status)) +
        # Add bars
        geom_boxplot() +
        # Use different sets of colors for different MajorTypes
        scale_fill_manual(values = colors) +
        # Rotate the x-axis labels to make them readable
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        # Add title and labels
        labs(
            title = "Fraction of each CNP by RFS_status and Tissue",
            x = "CNP", y = "Fraction", fill = "RFS Status"
        ) +
        # Separate the plot by Tissue and MajorType
        facet_grid(~Tissue, scales = "free") +
        # Hide the borders of the facet grid, make the background of the facet to be empty
        theme(
            strip.background = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank()
        ) +
        stat_compare_means(aes(group = RFS_status), method = "t.test", hide.ns = T, label = "p.signif")


    pdf(paste0(savePath, "Abundance Boxplot of CNP.pdf"), height = 4, width = 8)
    print(p)
    dev.off()
    return(NULL)
}

## Comapre the Celluar pattern into groups
CompareCellularPattern <- function(sce_, sep, countcol, n_cluster, clinicalFeatures, savePath) {
    groups <- names(table(colData(sce_)[, sep]))
    cat("The category is: ", groups, " in ", sep, "\n")

    sce1 <- sce_[, colData(sce_)[, sep] == groups[1]]
    sce2 <- sce_[, colData(sce_)[, sep] == groups[2]]

    ## ROI-level Boxplot
    abundance1 <- GetAbundance(sceobj = sce1, countcol = countcol, clinicalFeatures = clinicalFeatures, is.reuturnMeans = F)
    abundance2 <- GetAbundance(sceobj = sce2, countcol = countcol, clinicalFeatures = clinicalFeatures, is.reuturnMeans = F)

    BoxPlotForCellular(abundance1, abundance2, sep = sep, groups = groups, valueCol = c(1:n_cluster), savePath)

    return(NULL)
}

## Cellular pattern fraction
CNPFraction <- function(countDF, groupBy, xCol, savePath) {
    groups <- names(table(countDF[, groupBy]))

    DFg1 <- countDF[which(countDF[, groupBy] == groups[1]), xCol[1]:xCol[2]]
    DFg2 <- countDF[which(countDF[, groupBy] == groups[2]), xCol[1]:xCol[2]]

    DFg1 <- apply(DFg1, MARGIN = 2, FUN = "sum")
    DFg2 <- apply(DFg2, MARGIN = 2, FUN = "sum")

    DFg1 <- DFg1 / sum(DFg1)
    DFg2 <- DFg2 / sum(DFg2)

    ValueVec <- c(as.numeric(DFg1), as.numeric(DFg2))
    CNPVec <- rep(names(DFg1), times = 2)
    GroupVec <- rep(groups, each = length(DFg1))

    plotdf <- cbind.data.frame(ValueVec, CNPVec, GroupVec)

    p <- ggplot(data = plotdf, mapping = aes(x = "Content", y = ValueVec, fill = CNPVec)) +
        geom_bar(stat = "identity", position = "stack", width = 1) +
        coord_polar(theta = "y") +
        labs(x = "", y = "", title = "") +
        theme(axis.text = element_blank()) +
        theme(axis.ticks = element_blank()) +
        theme_bw() +
        facet_grid(GroupVec ~ .)

    pdf(paste0(savePath, "CNP fraction in ", groupBy, ".pdf"), height = 6, width = 8)
    print(p)
    dev.off()

    return(NULL)
}

## get abundace
GetAbundance <- function(sceobj, countcol, clinicalFeatures, is.fraction = TRUE, is.reuturnMeans = FALSE) {
    cellMeta <- colData(sceobj)

    ## ROI, major celltype and cell subtype names and other clinical information
    ROIs <- names(table(cellMeta$ID))

    SubTypes <- names(table(cellMeta[, countcol]))
    alltypes <- unique(c(SubTypes))

    CellCountMat <- matrix(data = NA, nrow = length(ROIs), ncol = (length(alltypes)))
    CellCountMat <- as.data.frame(CellCountMat)
    rownames(CellCountMat) <- ROIs
    colnames(CellCountMat) <- c(alltypes)

    for (ROI in ROIs) {
        sceTemp <- sceobj[, sceobj$ID == ROI]

        coldataTemp <- colData(sceTemp)
        cellnum <- nrow(coldataTemp)

        ## count cells
        SubTem <- as.data.frame(t(table(coldataTemp[, countcol])))

        if (is.fraction) {
            CellCountMat[match(ROI, rownames(CellCountMat)), match(SubTem$Var2, colnames(CellCountMat))] <- SubTem$Freq / cellnum
        }
        if (!is.fraction) {
            CellCountMat[match(ROI, rownames(CellCountMat)), match(SubTem$Var2, colnames(CellCountMat))] <- SubTem$Freq
        }
    }

    ## match the row of clinical and plotdf
    CellCountMat$PID <- as.vector(sapply(rownames(CellCountMat), function(x) {
        strsplit(x, "_")[[1]][1]
    }))

    for (i in clinicalFeatures) {
        CellCountMat[, i] <- cellMeta[match(rownames(CellCountMat), cellMeta$ID), ][, i]
    }

    for (i in 1:ncol(CellCountMat)) {
        CellCountMat[, i][is.na(CellCountMat[, i])] <- 0
    }

    if (!is.reuturnMeans) {
        return(CellCountMat)
    }
    if (is.reuturnMeans) {
        PIDs <- names(table(CellCountMat$PID))
        CellCountMat2 <- matrix(data = 0, nrow = length(PIDs), ncol = ncol(CellCountMat))
        rownames(CellCountMat2) <- PIDs
        colnames(CellCountMat2) <- colnames(CellCountMat)
        CellCountMat2 <- as.data.frame(CellCountMat2)

        for (i in PIDs) {
            Temp <- subset(CellCountMat, PID == i)
            numMetaFeature <- length(clinicalFeatures) + 1
            expTemp <- Temp[, 1:(ncol(Temp) - numMetaFeature)]
            MeanexpTemp <- apply(expTemp, MARGIN = 2, FUN = "mean")
            CellCountMat2[i, 1:(ncol(Temp) - numMetaFeature)] <- MeanexpTemp
            CellCountMat2[i, (numMetaFeature + 1):ncol(Temp)] <- Temp[1, (numMetaFeature + 1):ncol(Temp)]
        }
        return(CellCountMat2)
    }
}

## Boxplot For cellular pattern
BoxPlotForCellular <- function(mat1, mat2, sep, groups, valueCol, savePath) {
    plotdf1 <- pivot_longer(mat1, cols = valueCol, values_to = "Abundance", names_to = "Pattern")
    plotdf2 <- pivot_longer(mat2, cols = valueCol, values_to = "Abundance", names_to = "Pattern")

    plotdf <- rbind(plotdf1, plotdf2)
    plotdf$Group <- c(rep(groups[1], nrow(plotdf1)), rep(groups[2], nrow(plotdf2)))

    p <- ggplot(data = plotdf, aes(x = Pattern, y = Abundance, fill = Group)) +
        geom_boxplot(alpha = 0.7) +
        scale_y_continuous(name = "Abundance") +
        scale_x_discrete(name = "Cell Neighborhood Pattern") +
        theme_bw() +
        theme(
            plot.title = element_text(size = 14, face = "bold"),
            text = element_text(size = 12),
            axis.title = element_text(face = "bold"),
            axis.text.x = element_text(size = 11, angle = 90)
        ) +
        scale_fill_manual(values = c("#5494cc", "#e18283")) +
        stat_compare_means(aes(group = Group), label.y = 0.8, method = "t.test", label = "p.format")

    pdf(paste0(savePath, "Cellular Neighborhood pattern difference of ", sep, ".pdf"), height = 6, width = 12)
    print(p)
    dev.off()

    return(NULL)
}

## KM curve For cellular pattern
KMForCNP <- function(df, cluster, savePath) {
    cutpoint <- surv_cutpoint(data = df, time = "RFS_time", event = "RFS_status", variables = cluster)
    cutpoint <- summary(cutpoint)$cutpoint
    df[, 1] <- ifelse(df[, 1] >= cutpoint, "high", "low")

    colnames(df) <- c("Abundance", "RFS_status", "RFS_time")
    fit <- survfit(Surv(RFS_time, RFS_status) ~ Abundance, data = df)
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


    pdf(savePath, width = 8, height = 6)
    print(p)
    dev.off()

    return(NULL)
}

## clutering via certain markers
Reclustering <- function(sce, markers, ReMajorType, ReclusterName, ReSubType = NULL, PatternCol = NULL, ncluster = 10, savePath) {
    ## extract major types
    if (is.null(ReSubType)) {
        sce_ <- sce[, colData(sce)$MajorType %in% ReMajorType]
    } else {
        sce_ <- sce[, colData(sce)$SubType %in% ReSubType]
    }

    exp <- assay(sce_)
    exp <- exp[markers, ]

    ## K-means clustering
    exp <- t(exp) ## row should be sample
    set.seed(619)
    fit <- kmeans(exp, centers = ncluster, nstart = 50, iter.max = 100000)

    table(fit$cluster)
    clusters <- fit$cluster

    colData(sce_)[, ReclusterName] <- clusters

    ## T-sne visualization
    if (nrow(exp) <= 15000) {
        sampleidx <- c(1:nrow(exp))
    } else {
        sampleidx <- sample(1:nrow(exp), size = 15000, replace = F)
    } ### sample 15k cells to visualize}

    exp_sample <- exp[sampleidx, ]
    tsne <- Rtsne(exp_sample, dims = 2, PCA = F, verbose = F, max_iter = 500, check_duplicates = F)
    tsne_coor <- data.frame(tSNE1 = tsne$Y[, 1], tSNE2 = tsne$Y[, 2])
    tsne_coor$cluster <- as.factor(fit$cluster[sampleidx])
    tsne_coor$group <- ifelse(sce_$RFS_status[sampleidx] == 0, "Non-Relapse", "Relapse")

    colour <- c(brewer.pal(8, "Set2"), brewer.pal(9, "Set1"), brewer.pal(10, "Set3"))

    centers <- tsne_coor[, c("tSNE1", "tSNE2", "cluster")] %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(x = median(x = tSNE1), y = median(x = tSNE2))

    p <- ggplot(tsne_coor, aes(tSNE1, tSNE2)) +
        geom_point(aes(color = cluster), size = 0.5) +
        scale_fill_manual(values = colour) +
        guides(color = guide_legend(override.aes = list(size = 8, alpha = 1))) +
        theme_classic() +
        geom_text(data = centers, aes(x, y, label = cluster)) +
        facet_grid(~group)

    pdf(paste0(savePath, "tSNE reclustering of ", ReclusterName, ".pdf"), height = 4, width = 8)
    print(p)
    dev.off()

    ## Plot the marker of each cluster
    BubbleForcluterMarker(sce_, ReclusterName, markers, savePath)

    ## plot all markers expression value
    if (!dir.exists(paste0(savePath, "marker TSNE of ", ReclusterName, "/"))) {
        dir.create(paste0(savePath, "marker TSNE of ", ReclusterName, "/"))
    }
    PlotMarkerOnTSNE(exp_sample, tsne_coor, ReclusterName, paste0(savePath, "marker TSNE of ", ReclusterName, "/"))

    ## The relationship between re-clustering, origin cell subtype and Cellular pattern
    if (!is.null(PatternCol)) {
        SubtypeInReclustering(sce_, reclusteringCol = ReclusterName, OrigintypeCol = "SubType", PatternCol = PatternCol, savePath)
    }
    return(sce_)
}

## Bubble plot for visualize the reclustering markers
BubbleForcluterMarker <- function(sce_, colname1, markers, savePath) {
    exp <- assay(sce_)
    exp <- exp[markers, ]

    labels <- names(table(colData(sce_)[, colname1]))
    plotdf <- matrix(data = NA, nrow = length(labels), ncol = length(markers))
    rownames(plotdf) <- labels
    colnames(plotdf) <- markers

    for (i in 1:nrow(plotdf)) {
        cluterTemp <- rownames(plotdf)[i]
        idxTemp <- colData(sce_)[, colname1] == cluterTemp
        expTemp <- exp[markers, idxTemp]
        expTemp <- apply(expTemp, MARGIN = 1, FUN = "mean")

        plotdf[i, ] <- expTemp
    }

    MarkerIntensity <- as.numeric(plotdf)
    cluterID <- rep(rownames(plotdf), times = ncol(plotdf))
    Marker <- rep(colnames(plotdf), each = nrow(plotdf))

    plotdf2 <- data.frame("MarkerIntensity" = MarkerIntensity, "cluterID" = cluterID, "Marker" = Marker)

    p <- ggplot(plotdf2, aes(x = Marker, y = cluterID, size = MarkerIntensity, color = cluterID)) +
        geom_point() +
        theme(
            panel.background = element_blank(),
            panel.grid.major = element_line(colour = "white"),
            panel.border = element_rect(colour = "white", fill = NA),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
        ) +
        guides(color = guide_legend(override.aes = list(size = 8, alpha = 1)))
    pdf(paste0(savePath, "Bubble plot of ", colname1, ".pdf"), height = 6, width = 8)
    print(p)
    dev.off()

    return(NULL)
}

## 0-1 normlization
zero2oneNor <- function(vec) {
    vec <- as.numeric(vec)
    max_ <- max(vec)
    min_ <- min(vec)

    return((vec - min_) / (max_ - min_))
}

## Plot the marker expression on T-sne
PlotMarkerOnTSNE <- function(expDF, tsneDF, ReclusterName, savePath) {
    markers <- colnames(expDF)
    # plotdf <- as.data.frame(matrix(data = NA, nrow = 0, ncol = (ncol(tsneDF) + 1)))

    for (marker in markers) {
        plotdfTemp <- cbind(tsneDF, expDF[, marker])
        colnames(plotdfTemp) <- c("tSNE1", "tSNE2", "cluster", "group", "Intensity")
        plotdfTemp[, "Intensity"] <- as.numeric(plotdfTemp[, "Intensity"])
        plotdfTemp$Marker <- rep(marker, times = nrow(plotdfTemp))
        # plotdf <- rbind(plotdf, plotdfTemp)


        # p <- ggplot(plotdf, aes(tSNE1, tSNE2)) +
        p <- ggplot(plotdfTemp, aes(tSNE1, tSNE2)) +
            geom_point(aes(color = Intensity), size = 0.5) +
            scale_colour_gradient(low = "grey", high = "#EE0000") +
            theme_classic() +
            facet_grid(Marker ~ group)

        pdf(paste0(savePath, marker, " expression on tSNE reclustering.pdf"), height = 3, width = 6)
        print(p)
        dev.off()
    }
}

## Plot the contribution of cell subtype and cellular pattern in re-clustering
SubtypeInReclustering <- function(sce_, reclusteringCol, OrigintypeCol, PatternCol, savePath) {
    recluster <- colData(sce_)[, reclusteringCol]
    oritype <- colData(sce_)[, OrigintypeCol]
    pattern <- colData(sce_)[, PatternCol]

    ## Subtype in re-clustering
    reclusters <- names(table(recluster))
    cluster2subtypeDF <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))
    for (clu in reclusters) {
        idxTemp <- recluster %in% clu
        oritypeTemp <- as.data.frame(table(oritype[idxTemp]))
        oritypeTemp <- cbind(Recluster = rep(clu, nrow(oritypeTemp)), oritypeTemp)
        cluster2subtypeDF <- rbind(cluster2subtypeDF, oritypeTemp)
    }
    colnames(cluster2subtypeDF) <- c("Recluster", "CellSubtype", "Counts")

    color <- c(brewer.pal(n = 8, "Set1"), brewer.pal(n = 8, "Set2"), brewer.pal(n = 8, "Set3"))
    p <- ggplot(data = cluster2subtypeDF, aes(x = Recluster, y = Counts)) +
        geom_bar(aes(fill = CellSubtype), stat = "identity", width = 0.9) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.title = element_text(size = 12, face = "bold"),
            axis.text.x = element_text(angle = 90, hjust = 1),
            plot.margin = unit(rep(3, 4), "lines")
        ) +
        scale_fill_manual("CellSubtype", values = color) +
        coord_flip()
    pdf(paste0(savePath, "Subtypes in ", reclusteringCol, ".pdf"), height = 8, width = 6)
    print(p)
    dev.off()

    ## re-clustering in pattern
    patterns <- names(table(pattern))
    pattern2reclusterDF <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))
    for (pat in patterns) {
        idxTemp <- pattern %in% pat
        reclusterTemp <- as.data.frame(table(recluster[idxTemp]))
        reclusterTemp <- cbind("Pattern" = rep(pat, nrow(reclusterTemp)), reclusterTemp)
        pattern2reclusterDF <- rbind(pattern2reclusterDF, reclusterTemp)
    }
    colnames(pattern2reclusterDF) <- c("Pattern", "Recluster", "Counts")

    p <- ggplot(data = pattern2reclusterDF, aes(x = Pattern, y = Counts)) +
        geom_bar(aes(fill = Recluster), stat = "identity", width = 0.9) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.title = element_text(size = 12, face = "bold"),
            axis.text.x = element_text(angle = 90, hjust = 1),
            plot.margin = unit(rep(3, 4), "lines")
        ) +
        scale_fill_manual("Recluster", values = color) +
        coord_flip()
    pdf(paste0(savePath, reclusteringCol, " in ", PatternCol, ".pdf"), height = 8, width = 6)
    print(p)
    dev.off()

    return(NULL)
}

## Plot certain certain reclustering types in cellular pattern
PlotCertainTypeinPattern <- function(sce_, Col1, types1, Col2, groupCol, savePath) {
    Vec1 <- as.character(colData(sce_)[, Col1])
    Vec2 <- as.character(colData(sce_)[, Col2])
    VecGroup <- as.character(colData(sce_)[, groupCol])

    idx <- Vec1 %in% as.character(types1)

    Vec2 <- Vec2[idx]
    VecGroup <- VecGroup[idx]

    ## Count
    Vec2Names <- names(table(Vec2))
    VecGroupNames <- names(table(VecGroup))

    plotdf <- matrix(data = NA, nrow = 0, ncol = 3)
    plotdf <- as.data.frame(plotdf)

    for (name1 in Vec2Names) {
        for (name2 in VecGroupNames) {
            idx1 <- Vec2 %in% name1
            idx2 <- VecGroup %in% name2
            idx_ <- idx1 & idx2
            idx_ <- sum(as.numeric(idx_))
            VecTemp <- c(idx_, name1, name2)
            plotdf <- rbind(plotdf, VecTemp)
        }
    }
    colnames(plotdf) <- c("Counts", "CellularPattern", "Relapse")
    plotdf$Counts <- as.numeric(plotdf$Counts)

    ## plot
    p <- ggplot(data = plotdf, aes(x = CellularPattern, y = Counts)) +
        geom_bar(aes(fill = Relapse), stat = "identity", width = 0.9, position = "dodge") +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.title = element_text(size = 12, face = "bold"),
            axis.text.x = element_text(angle = 90, hjust = 1),
            plot.margin = unit(rep(3, 4), "lines"),
            legend.position = "bottom", legend.box = "horizontal"
        ) +
        scale_fill_brewer(palette = "Paired") +
        coord_flip()
    pdf(paste0(savePath, Col1, " ", as.character(types1), " in ", Col2, ".pdf"), height = 8, width = 6)
    print(p)
    dev.off()

    return(NULL)
}

## Plot CNP on cell-level
## get cell coordinate
getCoordinate <- function(sce_) {
    if (!"Position" %in% colnames(colData(sce_))) {
        cat("There is not column named Position!", "\n")
        return(NULL)
    }

    position <- sapply(sce_$Position, function(a) {
        strsplit(a, ",")
    })
    names(position) <- NULL

    x <- lapply(position, function(b) {
        return(as.numeric(strsplit(b[1], "\\(")[[1]][2]))
    })
    y <- lapply(position, function(b) {
        return(as.numeric(strsplit(b[2], "\\)")[[1]][1]))
    })

    x <- unlist(x)
    y <- unlist(y)

    return(list(x, y))
}

## plot celltypes on cell-level
PlotCNP <- function(sce, ROI, TypeCol, SavePath) {
    sce_ <- sce[, sce$ID == ROI]

    coorList <- getCoordinate(sce_)

    celltypes <- colData(sce_)[, TypeCol]

    plotdf <- as.data.frame(matrix(data = NA, nrow = length(celltypes), ncol = 3))
    colnames(plotdf) <- c("x", "y", "CNP")
    plotdf["x"] <- coorList[[1]]
    plotdf["y"] <- coorList[[2]]
    plotdf["CNP"] <- celltypes

    myPalette <- brewer.pal(12, "Set3")

    p <- ggplot(plotdf, aes(x = x, y = y)) +
        geom_point(aes(color = CNP), size = 1) +
        scale_colour_manual(values = myPalette, drop = F) +
        labs(title = paste0(ROI)) +
        theme_test()

    pdf(paste0(SavePath, "_CNP.pdf"), height = 6, width = 8)
    print(p)
    dev.off()

    return(NULL)
}

## Phenotype-associated cell label Boxplot
BoxPlotForPhenoAssCell <- function(plotdf, savePath) {
    plotdf$Relapse <- ifelse(plotdf$Relapse == 0, "Non-Relapse", "Relapse")

    p <- ggplot(data = plotdf, aes(x = Pattern, y = Abundance, fill = Relapse)) +
        geom_boxplot(alpha = 0.7) +
        scale_y_continuous(name = "Abundance") +
        scale_x_discrete(name = "Cell Neighborhood Pattern") +
        theme_bw() +
        theme(
            plot.title = element_text(size = 14, face = "bold"),
            text = element_text(size = 12),
            axis.title = element_text(face = "bold"),
            axis.text.x = element_text(size = 11, angle = 90)
        ) +
        scale_fill_manual(values = c("#5494cc", "#e18283")) +
        stat_compare_means(aes(group = Relapse), label.y = max(plotdf$Abundance), method = "t.test", label = "p.signif")

    pdf(paste0(savePath, "Boxplot of Pheno-associated celllabel.pdf"), height = 6, width = 8)
    print(p)
    dev.off()

    return(NULL)
}

## Phenotype-associated cell label KM curve
SurvivalForPhenoAssCell <- function(plotdf, savePath) {
    ## merge patients
    plotdf$PID <- sapply(as.character(plotdf$ROI), function(x) {
        return(strsplit(x, split = "_")[[1]][1])
    })

    PIDs <- unique(plotdf$PID)
    plotdf2 <- as.data.frame(matrix(data = NA, nrow = length(PIDs), ncol = 3))
    rownames(plotdf2) <- PIDs

    plotdf <- subset(plotdf, Pattern == "Pheno_plus")

    for (pid in PIDs) {
        temp <- subset(plotdf, PID == pid)
        plotdf2[pid, ] <- c(mean(as.numeric(temp$Abundance)), temp[1, 4], temp[1, 5])
    }
    colnames(plotdf2) <- c("Abundance", "RFSS", "RFS_time")

    plotdf2$Abundance <- ifelse(plotdf2$Abundance >= median(plotdf2$Abundance), "High", "Low")

    plotdf2$RFSS <- as.numeric(plotdf2$RFSS)
    plotdf2$RFS_time <- as.numeric(plotdf2$RFS_time)

    fit <- survfit(Surv(RFS_time, RFSS) ~ Abundance, data = plotdf2)
    p <- ggsurvplot(fit,
        data = plotdf2,
        linetype = c("solid", "solid"),
        surv.median.line = "hv", surv.scale = "percent",
        pval = T, risk.table = T,
        conf.int = T, conf.int.alpha = 0.1, conf.int.style = "ribbon",
        risk.table.y.text = T,
        palette = c("#3300CC", "#CC3300"),
        xlab = "Recurrence time"
    )

    pdf(paste0(savePath, "Suvival analysis for phenotype-associated cluster.pdf"), height = 6, width = 8)
    print(p)
    dev.off()

    return(NULL)
}


## calcualte FC
FCandPvalueCal <- function(mat, xCol, yCol, groups = NULL, need.sample = FALSE, sample.size = 2000) {
    set.seed(619)
    if (need.sample) {
        if (nrow(mat) > sample.size) {
            idx <- sample(1:nrow(mat), size = sample.size, replace = FALSE)
            mat <- mat[idx, ]
        }
    }
    if (is.null(groups)) {
        groups <- names(table(mat[, yCol]))
        groups <- as.character(sort(groups, decreasing = F))
    }
    cat("The Group in DEGs are", groups[2], "versus", groups[1], "\n")
    if (length(groups) < 2) {
        return(0)
    }

    returnMat <- matrix(data = NA, nrow = 0, ncol = 3)
    returnMat <- as.data.frame(returnMat)

    group1mat <- mat[which(mat[, yCol] == groups[1]), ]
    group2mat <- mat[which(mat[, yCol] == groups[2]), ]

    for (i in xCol[1]:xCol[2]) {
        typeTemp <- colnames(mat)[i]

        v1 <- group1mat[, i]
        v2 <- group2mat[, i]

        ## relaps versus no relaps
        foldchange <- mean(v2) / mean(v1)
        pvalue <- t.test(v2, v1)$p.value

        returnMat <- rbind(returnMat, c(typeTemp, foldchange, pvalue))
    }

    colnames(returnMat) <- c("Celltype", "Foldchange", "P.value")
    return(returnMat)
}

## Plot volcano plot
VolcanoPlot <- function(df, pthreshold = 0.05, fcthreshold = 1.4, feature, filename = NULL, Qvalue = F) {
    ## Fold change
    df$Foldchange <- as.numeric(df$Foldchange)
    df$P.value <- as.numeric(df$P.value)

    ## p value
    if (Qvalue) {
        df$Q.value <- p.adjust(df$P.value, method = "BH")
        df$change <- as.factor(ifelse(df$Q.value < pthreshold & abs(log2(df$Foldchange)) > log2(fcthreshold),
            ifelse(log2(df$Foldchange) > log2(fcthreshold), "Up-regulate", "Down-regulate"), "Non-significant"
        ))

        ## label
        df$label <- ifelse(df[, "Q.value"] < pthreshold & abs(log2(df$Foldchange)) > log2(fcthreshold), as.character(df[, 1]), "")

        ## plot
        p.vol <- ggplot(
            data = df,
            aes(x = log2(Foldchange), y = -log10(Q.value), colour = change, fill = change)
        ) +
            scale_color_manual(values = c("Down-regulate" = "blue", "Non-significant" = "grey", "Up-regulate" = "red")) +
            geom_point(alpha = 0.4, size = 3.5) +
            geom_text_repel(aes(x = log2(Foldchange), y = -log10(Q.value), label = label),
                size = 3,
                box.padding = unit(0.6, "lines"), point.padding = unit(0.7, "lines"),
                segment.color = "black", show.legend = FALSE
            ) +
            geom_vline(xintercept = c(-(log2(fcthreshold)), (log2(fcthreshold))), lty = 4, col = "black", lwd = 0.8) +
            geom_hline(yintercept = -log10(pthreshold), lty = 4, col = "black", lwd = 0.8) +
            theme_bw() +
            labs(x = "log2(Fold Change)", y = "-log10(Q value)", title = paste0("Volcano Plot of Different Expression Markers in ", feature)) +
            theme(
                axis.text = element_text(size = 11), axis.title = element_text(size = 13),
                plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
                legend.text = element_text(size = 11), legend.title = element_text(size = 13),
                plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
            )
    } else {
        df$change <- as.factor(ifelse(df$P.value < pthreshold & abs(log2(df$Foldchange)) > log2(fcthreshold),
            ifelse(log2(df$Foldchange) > log2(fcthreshold), "Up-regulate", "Down-regulate"), "Non-significant"
        ))

        ## label
        df$label <- ifelse(df[, 3] < pthreshold & abs(log2(df$Foldchange)) > log2(fcthreshold), as.character(df[, 1]), "")

        ## plot
        p.vol <- ggplot(
            data = df,
            aes(x = log2(Foldchange), y = -log10(P.value), colour = change, fill = change)
        ) +
            scale_color_manual(values = c("Down-regulate" = "blue", "Non-significant" = "grey", "Up-regulate" = "red")) +
            geom_point(alpha = 0.4, size = 3.5) +
            geom_text_repel(aes(x = log2(Foldchange), y = -log10(P.value), label = label),
                size = 3,
                box.padding = unit(0.6, "lines"), point.padding = unit(0.7, "lines"),
                segment.color = "black", show.legend = FALSE
            ) +
            geom_vline(xintercept = c(-(log2(fcthreshold)), (log2(fcthreshold))), lty = 4, col = "black", lwd = 0.8) +
            geom_hline(yintercept = -log10(pthreshold), lty = 4, col = "black", lwd = 0.8) +
            theme_bw() +
            labs(x = "log2(Fold Change)", y = "-log10(P value)", title = paste0("Volcano Plot of Different Expression Markers in ", feature)) +
            theme(
                axis.text = element_text(size = 11), axis.title = element_text(size = 13),
                plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
                legend.text = element_text(size = 11), legend.title = element_text(size = 13),
                plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
            )
    }
    ggsave(p.vol, filename = filename, width = 7, height = 7)

    return(NULL)
}

## Assign new label
AssignNewlabel <- function(sce_, allROIs, phenoLabel, ReclusterName, interstType, clinicalFeatures, cutoffType, cutoffValue = 10, numGroup = 2, is.reuturnsce = FALSE) {
    colData(sce_)[, phenoLabel] <- "Pheno_neg"
    colData(sce_)[, phenoLabel][which(colData(sce_)[, ReclusterName] %in% interstType)] <- "Pheno_pos"

    if (!is.reuturnsce) {
        CountMat <- GetAbundance(sce_, countcol = phenoLabel, clinicalFeatures = clinicalFeatures, is.fraction = FALSE, is.reuturnMeans = FALSE)
        if (cutoffType == "median") {
            Cutoff <- median(CountMat[, "Pheno_pos"])
            cat("The median cutoff is ", Cutoff, "\n")
        }
        if (cutoffType == "mean") {
            Cutoff <- mean(CountMat[, "Pheno_pos"])
            cat("The mean cutoff is ", Cutoff, "\n")
        }
        if (cutoffType == "manual") {
            Cutoff <- cutoffValue
            cat("The manual cutoff is ", Cutoff, "\n")
        }
        if (cutoffType == "best") {
            df <- CountMat[, c("Pheno_pos", "RFS_time", "RFS_status")]
            cutpoint <- surv_cutpoint(data = df, time = "RFS_time", event = "RFS_status", variables = "Pheno_pos")
            Cutoff <- summary(cutpoint)$cutpoint
            cat("The best cutoff is ", Cutoff, "\n")
        }

        RList <- list()
        RList[["sce"]] <- sce_
        if (numGroup == 2) {
            RList[["high_ROI"]] <- rownames(CountMat[which(CountMat[, "Pheno_pos"] >= Cutoff), ])
            RList[["low_ROI"]] <- allROIs[!allROIs %in% c(RList[["high_ROI"]])]
            if (length(RList[["low_ROI"]]) <= 5 | length(RList[["high_ROI"]]) <= 5) {
                Cutoff <- mean(CountMat[, "Pheno_pos"])
                RList[["high_ROI"]] <- rownames(CountMat[which(CountMat[, "Pheno_pos"] >= Cutoff), ])
                RList[["low_ROI"]] <- allROIs[!allROIs %in% c(RList[["high_ROI"]])]
                cat("Using best cutoff cause unbalance, shift to using mean as cutoff", "\n")
            }
        }
        if (numGroup == 3) {
            RList[["high_ROI"]] <- rownames(CountMat[which(CountMat[, "Pheno_pos"] >= Cutoff), ])
            RList[["low_ROI"]] <- rownames(CountMat[which(CountMat[, "Pheno_pos"] < Cutoff & CountMat[, "Pheno_pos"] > 0), ])
            RList[["none_ROI"]] <- allROIs[!allROIs %in% c(RList[["high_ROI"]], RList[["low_ROI"]])]
        }
        return(RList)
    } else {
        return(sce_)
    }
}

### get the results of permutation test
getResult <- function(ResultPath, ROIs, celltypes, p_threshold = 0.05) {
    ## Interaction
    numtypes <- length(celltypes)

    interDF <- matrix(data = 0, nrow = numtypes, ncol = numtypes)
    interDF <- as.data.frame(interDF)
    rownames(interDF) <- celltypes
    colnames(interDF) <- celltypes

    numID <- length(ROIs)

    for (ID in ROIs) {
        filepathTemp <- paste0(ResultPath, ID, "/ClosePvalue/ClosePvalue.csv")
        csvTemp <- read.csv(filepathTemp)
        rownames(csvTemp) <- csvTemp[, 1]
        csvTemp <- csvTemp[, -1]

        csvTemp <- ifelse(csvTemp <= p_threshold, 1, 0)

        colname_ <- colnames(csvTemp)
        colname_ <- sapply(colname_, function(x) {
            gsub(pattern = "\\.\\.", replacement = "+ ", x)
        })
        colname_ <- sapply(colname_, function(x) {
            gsub(pattern = "\\.", replacement = " ", x)
        })
        colname_ <- unname(colname_)
        rowname_ <- colname_

        for (i in 1:nrow(csvTemp)) {
            for (j in 1:ncol(csvTemp)) {
                if (csvTemp[i, j] == 1) {
                    interDF[rowname_[i], colname_[j]] <- interDF[rowname_[i], colname_[j]] + 1
                }
            }
        }
    }
    interDF <- round(interDF / numID, 4)

    ## Avoideness
    avoidDF <- matrix(data = 0, nrow = numtypes, ncol = numtypes)
    avoidDF <- as.data.frame(avoidDF)
    rownames(avoidDF) <- celltypes
    colnames(avoidDF) <- celltypes

    for (ID in ROIs) {
        filepathTemp <- paste0(ResultPath, ID, "/AvoidPvalue/AvoidPvalue.csv")
        csvTemp <- read.csv(filepathTemp)
        rownames(csvTemp) <- csvTemp[, 1]
        csvTemp <- csvTemp[, -1]

        csvTemp <- ifelse(csvTemp <= p_threshold, 1, 0)

        colname_ <- colnames(csvTemp)
        colname_ <- sapply(colname_, function(x) {
            gsub(pattern = "\\.\\.", replacement = "+ ", x)
        })
        colname_ <- sapply(colname_, function(x) {
            gsub(pattern = "\\.", replacement = " ", x)
        })
        rowname_ <- colname_

        for (i in 1:nrow(csvTemp)) {
            for (j in 1:ncol(csvTemp)) {
                if (csvTemp[i, j] == 1) {
                    avoidDF[rowname_[i], colname_[j]] <- avoidDF[rowname_[i], colname_[j]] + 1
                }
            }
        }
    }

    avoidDF <- round(avoidDF / numID, 4)

    ## Merge
    labelDF <- matrix(data = 0, nrow = numtypes, ncol = numtypes)
    labelDF <- as.data.frame(labelDF)
    rownames(labelDF) <- celltypes
    colnames(labelDF) <- celltypes

    MergeDF <- labelDF

    for (i in 1:nrow(labelDF)) {
        for (j in 1:ncol(labelDF)) {
            if (interDF[i, j] > avoidDF[i, j]) {
                MergeDF[i, j] <- interDF[i, j]
                labelDF[i, j] <- 1
            } else {
                MergeDF[i, j] <- avoidDF[i, j]
                labelDF[i, j] <- -1
            }
        }
    }

    return(list("MergeDF" = MergeDF, "LabelDF" = labelDF))
}

## Double heatmap
DoubleHeat <- function(MergeDF1, labelDF1, group1, MergeDF2, labelDF2, group2, plot = "circle", savePath) {
    DF1 <- MergeDF1 * labelDF1
    DF2 <- MergeDF2 * labelDF2

    plotdf <- matrix(data = 0, nrow = nrow(DF1) * ncol(DF1), ncol = 4)
    plotdf <- as.data.frame(plotdf)

    plotdf[, 1] <- rep(rownames(DF1), times = ncol(DF1))
    plotdf[, 2] <- rep(colnames(DF1), each = nrow(DF1))
    plotdf[, 3] <- as.numeric(as.matrix(DF1))
    plotdf[, 4] <- as.numeric(as.matrix(DF2))

    colnames(plotdf) <- c("Celltype1", "Celltype2", "Interaction1", "Interaction2")

    if (plot == "heatmap") {
        p <- ggplot(plotdf, aes(x = Celltype1, y = Celltype2)) +
            geom_heat_tri(
                upper = Interaction1, lower = Interaction2,
                upper_name = c(group1), lower_name = c(group2),
                lower_colors = c("#075fd5", "white", "#fd6a78"),
                upper_colors = c("#075fd5", "white", "#fd6a78")
            ) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }
    if (plot == "circle") {
        p <- ggplot(plotdf, aes(x = Celltype1, y = Celltype2)) +
            geom_heat_circle(
                outside = Interaction2,
                inside = Interaction1,
                outside_colors = c("#075fd5", "white", "#fd6a78"),
                inside_colors = c("#075fd5", "white", "#fd6a78")
            ) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }
    if (plot == "groupheatmap") {
        plotdf2 <- cbind(rep(plotdf[, 1], times = 2), rep(plotdf[, 2], times = 2), c(plotdf[, 3], plotdf[, 4]), rep(c(group1, group2), each = nrow(plotdf)))
        plotdf2 <- as.data.frame(plotdf2)
        plotdf2[, 3] <- as.numeric(plotdf2[, 3])
        colnames(plotdf2) <- c("Celltype1", "Celltype2", "Interaction", "Group")

        p <- ggplot(plotdf2, aes(x = Celltype1, y = Celltype2)) +
            geom_tile(aes(fill = Interaction)) +
            scale_fill_gradient2(low = "#0000ff", high = "#ff0000", mid = "#ffffff") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            facet_grid(. ~ Group)

        pdf(savePath, height = 6, width = 12)
        print(p)
        dev.off()

        return(NULL)
    }

    pdf(savePath, height = 6, width = 8)
    print(p)
    dev.off()

    return(NULL)
}

# Double heatmap
TribleHeat <- function(MergeDF1, labelDF1, group1, MergeDF2, labelDF2, group2, MergeDF3, labelDF3, group3, savePath) {
    DF1 <- MergeDF1 * labelDF1
    DF2 <- MergeDF2 * labelDF2
    DF3 <- MergeDF3 * labelDF3

    plotdf <- matrix(data = 0, nrow = nrow(DF1) * ncol(DF1), ncol = 5)
    plotdf <- as.data.frame(plotdf)

    plotdf[, 1] <- rep(rownames(DF1), times = ncol(DF1))
    plotdf[, 2] <- rep(colnames(DF1), each = nrow(DF1))
    plotdf[, 3] <- as.numeric(as.matrix(DF1))
    plotdf[, 4] <- as.numeric(as.matrix(DF2))
    plotdf[, 5] <- as.numeric(as.matrix(DF3))

    colnames(plotdf) <- c("Celltype1", "Celltype2", "High-Interaction", "Low-Interaction", "None-Interaction")

    plotdf2 <- cbind(rep(plotdf[, 1], times = 3), rep(plotdf[, 2], times = 3), c(plotdf[, 3], plotdf[, 4], plotdf[, 5]), rep(c(group1, group2, group3), each = nrow(plotdf)))
    plotdf2 <- as.data.frame(plotdf2)
    plotdf2[, 3] <- as.numeric(plotdf2[, 3])
    colnames(plotdf2) <- c("Celltype1", "Celltype2", "Interaction", "Group")

    p <- ggplot(plotdf2, aes(x = Celltype1, y = Celltype2)) +
        geom_tile(aes(fill = Interaction)) +
        scale_fill_gradient2(low = "#0000ff", high = "#ff0000", mid = "#ffffff") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        facet_grid(. ~ Group)

    pdf(savePath, height = 6, width = 18)
    print(p)
    dev.off()



    return(NULL)
}

## Swarm plot for high- and low- marker groups certain celltype
AbundanceSwarmPlot <- function(AbundanceDF1, AbundanceDF2, AbundanceDF3, groupsName, celltype, marker, savePath, numGroup = 3) {
    if (numGroup == 3) {
        if (!(celltype %in% colnames(AbundanceDF1)) | !(celltype %in% colnames(AbundanceDF2)) | !(celltype %in% colnames(AbundanceDF3))) {
            cat("Celltype ", celltype, " not in matrix!", "\n")
            return(NULL)
        }
        Counts1 <- AbundanceDF1[, celltype]
        Counts2 <- AbundanceDF2[, celltype]
        Counts3 <- AbundanceDF3[, celltype]

        ROIsVec <- c(rownames(AbundanceDF1), rownames(AbundanceDF2), rownames(AbundanceDF3))
        CountsVec <- c(Counts1, Counts2, Counts3)
        GroupsVec <- c(rep(groupsName[1], length(Counts1)), rep(groupsName[2], length(Counts2)), rep(groupsName[3], length(Counts3)))

        plotdf <- cbind(ROIsVec, CountsVec, GroupsVec)
        plotdf <- as.data.frame(plotdf)
        colnames(plotdf) <- c("ROI", "Count", "Group")

        plotdf$Count <- as.numeric(plotdf$Count)

        ## Swarm plot
        p <- ggplot(data = plotdf, aes(x = Group, y = Count, color = Group)) +
            geom_violin(show.legend = FALSE) +
            geom_jitter(aes(fill = Group)) +
            theme_bw() +
            labs(x = NULL) +
            scale_color_manual(values = brewer.pal(3, "Paired")) +
            stat_compare_means(
                comparisons = list(c("high", "low"), c("low", "none"), c("high", "none")),
                method = "t.test"
            )
    }
    if (numGroup == 2) {
        if (!(celltype %in% colnames(AbundanceDF1)) | !(celltype %in% colnames(AbundanceDF2))) {
            cat("Celltype ", celltype, " not in matrix!", "\n")
            return(NULL)
        }
        Counts1 <- AbundanceDF1[, celltype]
        Counts2 <- AbundanceDF2[, celltype]

        ROIsVec <- c(rownames(AbundanceDF1), rownames(AbundanceDF2))
        CountsVec <- c(Counts1, Counts2)
        GroupsVec <- c(rep(groupsName[1], length(Counts1)), rep(groupsName[2], length(Counts2)))

        plotdf <- cbind(ROIsVec, CountsVec, GroupsVec)
        plotdf <- as.data.frame(plotdf)
        colnames(plotdf) <- c("ROI", "Count", "Group")

        plotdf$Count <- as.numeric(plotdf$Count)

        ## Swarm plot
        p <- ggplot(data = plotdf, aes(x = Group, y = Count, color = Group)) +
            geom_violin(show.legend = FALSE) +
            geom_jitter(aes(fill = Group)) +
            theme_bw() +
            labs(x = NULL) +
            scale_color_manual(values = brewer.pal(2, "Paired")) +
            stat_compare_means(
                comparisons = list(c("high", "low")),
                method = "t.test"
            )
    }

    pdf(paste0(savePath, "Swarmplot of ", celltype, " in ", marker, " group.pdf"), height = 4, width = 3)
    print(p)
    dev.off()

    return(NULL)
}

## Swarm plot for high- and low- marker groups certain celltype
MultipleAbundanceSwarmPlot <- function(AbundanceDF1, AbundanceDF2, AbundanceDF3, groupsName, PlotTypes, marker, numGroup = 3, style = "violin") {
    if (numGroup == 2) {
        Counts1 <- AbundanceDF1[, PlotTypes]
        Counts2 <- AbundanceDF2[, PlotTypes]

        ROIsVec <- c(rep(rownames(Counts1), times = ncol(Counts1)), rep(rownames(Counts2), times = ncol(Counts2)))
        CountsVec <- c(as.numeric(as.matrix(Counts1)), as.numeric(as.matrix(Counts2)))
        CelltypeVec <- c(rep(colnames(Counts1), each = nrow(Counts1)), rep(colnames(Counts2), each = nrow(Counts2)))
        GroupsVec <- c(rep(groupsName[1], length(as.matrix(Counts1))), rep(groupsName[2], length(as.matrix(Counts2))))

        plotdf <- cbind(ROIsVec, CountsVec, CelltypeVec, GroupsVec)
        plotdf <- as.data.frame(plotdf)
        colnames(plotdf) <- c("ROI", "Count", "Celltype", "Group")

        plotdf$Count <- as.numeric(plotdf$Count)

        ## Swarm plot
        color <- pal_aaas("default")(10)[c(5, 7)]
        p <- ggplot(data = plotdf, aes(x = Group, y = Count, color = Group, fill = Group))
        if (style == "violin") {
            p <- p + geom_violin(show.legend = FALSE, outlier.shape = NA, alpha = 0.8)
        }
        if (style == "box") {
            p <- p + geom_boxplot(show.legend = FALSE, outlier.shape = NA, alpha = 0.8)
        }
        p <- p +
            geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
            theme_bw() +
            labs(x = NULL, y = "Count", title = "Box Plot of Counts by Group and Cell Type") +
            scale_color_manual(values = color) +
            scale_fill_manual(values = color) +
            stat_compare_means(
                aes(group = Group),
                comparisons = list(c(groupsName[1], groupsName[2])),
                method = "t.test",
                label = "p.signif",
                label.y = max(plotdf$Count) * 0.95,
                size = 4
            ) +
            facet_wrap(~Celltype, nrow = 1) +
            theme(
                plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                axis.title = element_text(size = 12, face = "bold"),
                axis.text = element_text(size = 10),
                strip.text = element_text(size = 12, face = "bold"),
                strip.background = element_blank(),
                panel.border = element_rect(color = "black", fill = NA, size = 1),
                panel.grid.major = element_line(color = "gray", size = 0.5),
                panel.grid.minor = element_blank()
            )
    }
    return(p)
}

## Survival analysis for Phenotype associated label
SurvivalForPhenoAssoLabel <- function(AbundanceDF, GroupCol, time, status, marker, cutoffType = "best", manual = NULL, savePath) {
    AbundanceDF <- AbundanceDF[, c(GroupCol, time, status)]

    if (cutoffType == "best") {
        cutpoint <- surv_cutpoint(data = AbundanceDF, time = time, event = status, variables = GroupCol)
        cutpoint <- summary(cutpoint)$cutpoint
    }
    if (cutoffType == "mean") {
        cutpoint <- mean(AbundanceDF[, GroupCol])
    }
    if (cutoffType == "median") {
        cutpoint <- median(AbundanceDF[, GroupCol])
    }
    if (cutoffType == "manual") {
        cutpoint <- manual
    }

    GroupVec <- ifelse(AbundanceDF[, GroupCol] >= cutpoint, "High Abundance", "Low Abundance")

    plotdf <- cbind(AbundanceDF[, -1], GroupVec)
    plotdf <- as.data.frame(plotdf)
    colnames(plotdf) <- c(time, status, "Label")

    plotdf[, 1] <- as.numeric(plotdf[, 1])
    plotdf[, 2] <- as.numeric(plotdf[, 2])
    plotdf[, 3] <- as.factor(plotdf[, 3])

    ## KM
    fit <- survfit(Surv(RFS_time, RFS_status) ~ Label, data = plotdf)
    p <- ggsurvplot(fit,
        data = plotdf,
        linetype = c("solid", "solid"),
        surv.median.line = "hv", surv.scale = "percent",
        pval = T, risk.table = T,
        conf.int = T, conf.int.alpha = 0.1, conf.int.style = "ribbon",
        risk.table.y.text = T,
        palette = c("#CC3300", "#3300CC"),
        xlab = "Recurrence time"
    )

    pdf(paste0(savePath, "Suvival analysis for phenotype-associated PID of ", marker, ".pdf"), height = 6, width = 8)
    print(p)
    dev.off()

    return(NULL)
}

## take mean
TakeMeans <- function(df, valuecol, groupcol, pidcol) {
    PIDs <- names(table(df[, pidcol]))

    PIDVec <- c()
    ValueVec <- c()
    GroupVec <- c()

    for (i in 1:length(PIDs)) {
        PIDVec <- c(PIDVec, PIDs[i])
        DFtemp <- df[which(df[, pidcol] == PIDs[i]), ]

        ValueVec <- c(ValueVec, mean(DFtemp[, valuecol]))
        GroupVec <- c(GroupVec, DFtemp[1, groupcol])
    }

    returnDF <- cbind(ValueVec, PIDVec, GroupVec)
    returnDF <- as.data.frame(returnDF)

    colnames(returnDF) <- c("Abundance", "PID", "Group")
    returnDF[, 1] <- as.numeric(returnDF[, 1])
    returnDF[, 2] <- as.character(returnDF[, 2])
    returnDF[, 3] <- as.character(returnDF[, 3])

    return(returnDF)
}

## gene co-expression analysis
CoexpAnalysis <- function(sce_, reclusMarkers, ReclusterName, interstType, savePath) {
    sce__ <- sce_[, colData(sce_)[, ReclusterName] %in% interstType]

    exp <- assay(sce__)
    exp <- t(exp[reclusMarkers, ])

    ## calculate the correlation
    corMat <- rcorr(exp, type = "spearman")
    rMat <- corMat[[1]]
    pMat <- corMat[[3]]

    rd <- c()
    pd <- c()
    marker1 <- c()
    marker2 <- c()

    for (i in 1:nrow(rMat)) {
        for (j in i:ncol(rMat)) {
            rd <- c(rd, rMat[i, j])
            pd <- c(pd, pMat[i, j])
            marker1 <- c(marker1, rownames(pMat)[i])
            marker2 <- c(marker2, colnames(pMat)[j])
        }
    }

    pd <- ifelse(is.na(pd), "0", pd)

    rd <- as.character(cut(as.numeric(rd), breaks = c(-Inf, 0.2, 0.4, Inf), labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")))
    pd <- as.character(cut(as.numeric(pd), breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

    plotdf <- as.data.frame(cbind("Marker1" = marker1, "Marker2" = marker2, "r" = as.factor(rd), "p" = as.factor(pd)))

    ## heatmap
    p <- quickcor(exp, method = "spearman", type = "upper", cor.test = T, cluster.type = "all") +
        geom_square() +
        geom_mark(r = NA, sig.thres = 0.05, size = 3.5, colour = "black") +
        scale_fill_gradient2(high = "red", mid = "white", low = "blue") +
        anno_link(plotdf, aes(color = p, size = r)) +
        scale_size_manual(values = c(0.5, 1, 1.5, 2)) +
        scale_color_manual(values = c("green", "blue", "orange")) +
        guides(
            fill = guide_colorbar(title = "correlation", order = 1),
            size = guide_legend(title = "Spearman's r", order = 2),
            color = guide_legend(title = "Spearman's p", order = 3),
            linetype = "none"
        )

    pdf("test.pdf")
    print(p)
    dev.off()

    return(NULL)
}

## Interaction difference barplot
BarplotForInteraction <- function(ResultPath, celltype1, celltype2, ROIs1, ROIs2, savePath) {
    return(NULL)
}

## Calculate celltypes difference in different clinical groups
HeatmapForMarkersInGroups <- function(sce, markers, groupCol, typesCol = "SubType", adjust.p = T, FCthreshold = 1.5, savePath, return.mat = FALSE, sample.size = 2000) {
    celltypes <- names(table(colData(sce)[, typesCol]))

    MarkerVec <- c()
    TypeVec <- c()
    PVec <- c()
    FCVec <- c()

    i <- 1
    for (interstSubType in celltypes) {
        sce_ <- sce[, colData(sce)[, typesCol] %in% interstSubType]
        mat <- as.data.frame(t(assay(sce_)[markers, ]))
        mat$classLabel <- colData(sce_)[, groupCol]
        FCDF <- FCandPvalueCal(mat, xCol = c(1, length(markers)), yCol = (length(markers) + 1), need.sample = TRUE, sample.size = sample.size)

        FCVec <- c(FCVec, as.numeric(FCDF[, 2]))
        PVec <- c(PVec, as.numeric(FCDF[, 3]))
        i <- i + 1
    }

    MarkerVec <- rep(markers, times = length(celltypes))
    TypeVec <- rep(celltypes, each = length(markers))

    plotdf <- matrix(data = NA, nrow = length(MarkerVec), ncol = 4)
    plotdf <- as.data.frame(plotdf)
    plotdf[, 1] <- MarkerVec
    plotdf[, 2] <- TypeVec
    plotdf[, 3] <- PVec
    plotdf[, 4] <- ifelse(FCVec != 0, log2(FCVec), 0)
    plotdf[, 4] <- ifelse(abs(plotdf[, 4]) >= log2(FCthreshold), plotdf[, 4], 0)

    colnames(plotdf) <- c("Marker", "Subtype", "Pvalue", "log2FoldChange")

    if (adjust.p) {
        plotdf$Qvalue <- p.adjust(PVec, method = "BH")
        plotdf$Qlabel <- cut(plotdf$Qvalue, breaks = c(0, 0.001, 0.01, 0.05, 1), labels = c(8, 6, 4, 2))
        plotdf$Qlabel <- as.factor(as.numeric(as.character(plotdf$Qlabel)))

        p <- ggplot(data = plotdf, aes(x = Subtype, y = Marker)) +
            geom_point(aes(size = Qlabel, fill = log2FoldChange), shape = 22, color = "grey80")
    } else {
        plotdf$Plabel <- cut(plotdf$Pvalue, breaks = c(0, 0.001, 0.01, 0.05, 1), labels = c(8, 6, 4, 2))
        plotdf$Plabel <- as.factor(plotdf$Plabel)
        plotdf$Plabel <- as.factor(as.numeric(as.character(plotdf$Plabel)))

        p <- ggplot(data = plotdf, aes(x = Subtype, y = Marker)) +
            geom_point(aes(size = Plabel, fill = log2FoldChange), shape = 22, color = "grey80")
    }

    if (return.mat) {
        return(plotdf)
    }

    p <- p +
        scale_fill_gradient2(low = "#445393", high = "#EE2627", mid = "white") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 12),
            axis.text.y = element_text(size = 12, color = "black"),
            axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.1, color = "black", size = 10),
            axis.ticks = element_blank(),
            legend.key.size = unit(0.15, "inches")
        ) +
        labs(x = "", y = NULL) +
        scale_size_discrete(range = c(2, 8))

    pdf(paste0(savePath, "alltypes differential markers in relapse.pdf"), width = 8, height = 6)
    print(p)
    dev.off()

    return(NULL)
}

## Get Differential expression genes in certain groups
GetDEGinGroup <- function(sce_, groupCol, groupLevel, cellType, markers) {
    groupName <- names(table(colData(sce_)[, groupCol]))
    cat("The clinical group ", groupCol, " is ", groupName, "\n")
    idx <- names(table(colData(sce_)[, groupLevel]))

    df <- matrix(data = 0, nrow = (length(idx) * length(cellType)), ncol = (length(groupCol) + length(markers) + 2))
    df <- as.data.frame(df)

    z <- 1

    for (i in 1:length(idx)) {
        idTemp <- idx[i]
        sce_Temp <- sce_[, colData(sce_)[, groupLevel] == idTemp]
        for (j in 1:length(cellType)) {
            cellTemp <- cellType[j]
            sce_Temp2 <- sce_Temp[, colData(sce_Temp)$SubType == cellTemp]
            expTemp <- assay(sce_Temp2)
            expTemp <- apply(expTemp, MARGIN = 1, FUN = "mean")
            expTemp <- expTemp[markers]
            df[z, ] <- c(expTemp, colData(sce_Temp)[1, groupCol], idTemp, cellTemp)
            z <- z + 1
        }
    }
    colnames(df) <- c(markers, groupCol, "ID", "Subtype")

    dfLong <- pivot_longer(df, cols = -c(RFS_status, ID, Subtype), names_to = "Marker", values_to = "Intensity")
    dfLong$Intensity <- as.numeric(dfLong$Intensity)

    FCList <- list()
    for (celltypeTemp in cellType) {
        dfLongTemp <- subset(dfLong, Subtype == celltypeTemp)
        FCList[[celltypeTemp]] <- CalFC(dfLongTemp, groupCol = groupCol, groupName, nameCol = "Marker", valueCol = "Intensity")
        FCList[[celltypeTemp]]$Subtype <- rep(celltypeTemp, nrow(FCList[[celltypeTemp]]))
    }

    plotdf <- as.data.frame(FCList[[1]])
    for (i in 2:length(FCList)) {
        plotdf <- rbind(plotdf, FCList[[i]])
    }

    return(plotdf)
}

CalFC <- function(df_, groupCol, groupName, nameCol, valueCol, eps = 1e-8) {
    calName <- names(table(df_[, nameCol]))

    resultDF <- as.data.frame(matrix(data = 0, nrow = length(calName), ncol = 2))
    rownames(resultDF) <- calName

    for (i in 1:length(calName)) {
        markerTemp <- calName[i]
        v1 <- as.numeric(as.matrix(df_[subset(df_[, nameCol] == markerTemp & df_[, groupCol] == groupName[1]), valueCol]))
        v2 <- as.numeric(as.matrix(df_[subset(df_[, nameCol] == markerTemp & df_[, groupCol] == groupName[2]), valueCol]))

        FC <- (mean(v2) + eps) / (mean(v1) + eps)
        P <- t.test(v2, v1)$p.value
        resultDF[i, ] <- c(FC, P)
    }
    colnames(resultDF) <- c("FC", "Pvalue")
    resultDF$log2FC <- log2(resultDF$FC)
    resultDF$lgPvalue <- -log10(resultDF$Pvalue + eps)
    resultDF$Qvalue <- p.adjust(resultDF$Pvalue, "BH")
    resultDF$lgQvalue <- -log10(resultDF$Qvalue + eps)
    resultDF$Marker <- rownames(resultDF)

    return(resultDF)
}

## Calcualte the x and y from a position vector
GetXandY <- function(posVec) {
    xVec <- c()
    yVec <- c()

    spliteRes <- sapply(posVec, function(x) {
        return(strsplit(x, split = ",")[[1]])
    })

    for (i in seq(1, length(spliteRes), 2)) {
        xTemp <- strsplit(spliteRes[i], "\\(")[[1]][2]
        yTemp <- strsplit(spliteRes[i + 1], ")")[[1]][1]

        xVec <- c(xVec, as.numeric(xTemp))
        yVec <- c(yVec, as.numeric(yTemp))
    }

    cor <- as.data.frame(matrix(data = NA, ncol = 2, nrow = length(xVec)))
    cor[, 1] <- xVec
    cor[, 2] <- yVec

    colnames(cor) <- c("cor_x", "cor_y")
    return(cor)
}

## Dimention Reduction for building graph
DimenRec <- function(sce_, GroupFeature, markers, do.pca = FALSE, pca.dimen = 10, explaine.pca = FALSE, embeeding.method = "tsne", num.threads = 1) {
    m <- t(assay(sce_))
    m <- m[, match(markers, colnames(m))]
    label <- colData(sce_)[, GroupFeature]

    if (do.pca) {
        # Apply PCA
        pca_result <- prcomp(m, scale. = TRUE)

        # Print summary of the PCA result
        print(summary(pca_result))

        # If you want to use the transformed data (for example, using the first two principal components), you can access it like this:
        m <- pca_result$x[, 1:pca.dimen]

        # If you want to see the proportion of variance explained by each principal component, you can view it like this:
        if (explaine.pca) {
            explained_variance <- pca_result$sdev^2
            explained_variance_ratio <- explained_variance / sum(explained_variance)

            print(explained_variance_ratio)
        }
    }

    if (embeeding.method == "tsne") {
        # Apply t-SNE
        tsne_result <- Rtsne(m, dims = 2, perplexity = 30, check_duplicates = FALSE, num_threads = num.threads)

        # If you want to use the transformed data (for example, using the first two t-SNE dimensions), you can access it like this:
        m <- tsne_result$Y
    }

    m <- as.data.frame(m)
    colnames(m) <- c("cor_x", "cor_y")
    m$label <- label

    return(m)
}

## Build KNN graph and calculate the neighbors entropy and fraction
KGandFracEntro <- function(m, cordim = 2, k = 10, central.weight = NULL, multi.cores = 8) {
    classes <- m$label
    m <- m[, 1:cordim]

    LabelGroup <- sort(names(table(classes)), decreasing = T)
    cat("The labels are ", LabelGroup, "\n")
    cat("Build KNN Graph", "\n")

    # Find the k-nearest neighbors
    knn_result <- get.knn(m, k = k)

    cat("Calculate Entropy and Fraction in Neighbors", "\n")
    # Calculate the entropy and fractions for each vertex, using parallel computation
    results <- mclapply(1:nrow(m), function(i) {
        # Get the classes of the k-nearest neighbors and the vertex itself
        neighbor_classes <- c(rep(classes[i], central.weight), classes[knn_result$nn.index[i, ]])
        df <- as.data.frame(table(factor(neighbor_classes, levels = LabelGroup)))
        # Calculate the entropy, handle the case of single-class event
        entropy_val <- CalEntropy(df)

        # Calculate the fractions of neighbors that belong to each class
        fractions <- CalFraction(df)

        # Return the entropy and fraction of neighbors
        c(entropy_val, fractions)
    }, mc.cores = multi.cores)

    # Convert the results to a matrix
    results_matrix <- do.call(rbind, results)
    colnames(results_matrix) <- c("Entropy", "Fraction_class_1")
    return(results_matrix)
}

## Calculate entropy
CalEntropy <- function(df) {
    if (0 %in% df[, 2]) {
        value <- 0
    } else {
        allcount <- sum(df[, 2])
        p_x <- df[, 2] / allcount
        logp_x <- log2(p_x)
        value <- p_x[1] * logp_x[1] + p_x[2] * logp_x[2]
        value <- (-1) * value
    }
    return(value)
}

## Calculate Fraction
CalFraction <- function(df) {
    allcount <- sum(df[, 2])
    value <- round(df[1, 2] / allcount, 6)

    return(value)
}

## Generate null distribution
GeneNullDistri <- function(m, cordim, sample.size = 20, perm.times = 1000, seed = 619) {
    classes <- m$label
    m <- m[, 1:cordim]
    LabelGroup <- sort(names(table(classes)), decreasing = T)

    cat("Generate Null Distribution in size", sample.size, "for times", perm.times, "\n")

    set.seed(seed)
    len_class <- length(classes)

    fraction_NullDistribution <- c()

    for (i in 1:perm.times) {
        randidx <- sample.int(len_class, size = sample.size, replace = FALSE)
        randneighbor_classes <- classes[randidx]

        df <- as.data.frame(table(factor(randneighbor_classes, levels = LabelGroup)))
        # Calculate the fractions of neighbors that belong to each class
        fractions <- CalFraction(df)
        fraction_NullDistribution <- c(fraction_NullDistribution, fractions)
    }

    return(sort(fraction_NullDistribution))
}

## Plot the points label and entropy
PlotClassLabel <- function(sce_, Classlabels, Entropylabels, markers, sample.size = 1e4, num.threads = 1, SavePath1 = NULL, SavePath2 = NULL, SavePath3 = NULL, seed = 619) {
    ## Data prepare
    exp <- t(assay(sce_))
    exp <- exp[, match(markers, colnames(exp))]

    set.seed(seed)
    ## sample
    if (!is.null(sample.size)) {
        idx <- sample.int(length(Classlabels), size = sample.size)
    } else {
        idx <- 1:(length(Classlabels))
    }

    exp <- exp[idx, ]
    Classlabels <- Classlabels[idx]
    Entropylabels <- Entropylabels[idx]
    PIDlabels <- sce_$PID[idx]

    ## T-SNE
    tsne_result <- Rtsne(exp, dims = 2, perplexity = 30, check_duplicates = FALSE, num_threads = num.threads)
    coor <- tsne_result$Y

    Classlabels <- factor(Classlabels, levels = c("Background", "Pheno_pos", "Pheno_neg"))

    ## Class
    plotdf <- as.data.frame(matrix(data = NA, nrow = length(Classlabels), ncol = 3))
    colnames(plotdf) <- c("x", "y", "Identity")
    plotdf["x"] <- coor[, 1]
    plotdf["y"] <- coor[, 2]
    plotdf["Identity"] <- Classlabels

    myPalette <- c("grey", pal_nejm("default")(2))
    myPalette <- alpha(myPalette, ifelse(myPalette == "grey", 0.5, 1))

    p <- ggplot(plotdf, aes(x = x, y = y, color = Identity)) +
        geom_point(size = 1) +
        scale_colour_manual(values = myPalette) +
        theme_test()

    pdf(SavePath1, height = 6, width = 8)
    print(p)
    dev.off()

    ## Entropy
    plotdf <- as.data.frame(matrix(data = NA, nrow = length(Entropylabels), ncol = 3))
    colnames(plotdf) <- c("x", "y", "Identity")
    plotdf["x"] <- coor[, 1]
    plotdf["y"] <- coor[, 2]
    plotdf["Identity"] <- Entropylabels

    p <- ggplot(plotdf, aes(x = x, y = y, color = Identity)) +
        geom_point(size = 1) +
        scale_colour_viridis(option = "magma") +
        theme_test()

    pdf(SavePath2, height = 6, width = 8)
    print(p)
    dev.off()

    ## Markers
    if (!dir.exists(paste0(SavePath3, "marker TSNE of Treg/"))) {
        dir.create(paste0(SavePath3, "marker TSNE of Treg/"))
    }

    for (marker in markers) {
        plotdfTemp <- cbind(coor, exp[, marker])
        colnames(plotdfTemp) <- c("tSNE1", "tSNE2", "Intensity")
        plotdfTemp[, "Intensity"] <- as.numeric(plotdfTemp[, "Intensity"])
        plotdfTemp <- as.data.frame(plotdfTemp)

        p <- ggplot(plotdfTemp, aes(tSNE1, tSNE2)) +
            geom_point(aes(color = Intensity), size = 0.5) +
            scale_colour_gradient(low = "grey", high = "#EE0000") +
            theme_classic()

        pdf(paste0(SavePath3, "marker TSNE of Treg/", marker, " expression on tSNE reclustering.pdf"), height = 3, width = 4)
        print(p)
        dev.off()
    }

    ## Patient level Tsne
    plotdf <- as.data.frame(matrix(data = NA, nrow = length(PIDlabels), ncol = 3))
    colnames(plotdf) <- c("x", "y", "Identity")
    plotdf["x"] <- coor[, 1]
    plotdf["y"] <- coor[, 2]
    plotdf["Identity"] <- as.factor(PIDlabels)

    p <- ggplot(plotdf, aes(x = x, y = y, color = Identity)) +
        geom_point(size = 1) +
        scale_colour_manual(values = pal_igv("default")(51)) +
        theme_test()

    pdf(paste0(SavePath3, "T-SNE of phenotype associated label on PID.pdf"), height = 6, width = 8)
    print(p)
    dev.off()

    return(NULL)
}

## Calcualte the K-nn of certain cell subpopulation via spatstat
CalKNNBySpatstat <- function(sce, ROIs, CenType, NeighType = "SubType", k = 10, xlim = c(0, 1000), ylim = c(0, 1000), return.cellID = FALSE) {
    # Record the 10-nn cells of TC
    k_NeighborsList <- list()

    for (ROI in ROIs) {
        sce_ <- sce[, sce$ID == ROI]

        coor <- sce_$Position
        coor_df <- GetXandY(coor)
        coor_df$SubType <- colData(sce_)[, NeighType]

        # head(coor_df)
        TargetTypeidx <- (colData(sce_)[, NeighType] == CenType)
        if (length(TargetTypeidx) == 0) {
            next
        }

        # Convert your data frame to a ppp object
        points_ppp <- ppp(coor_df$cor_x, coor_df$cor_y, xrange = xlim, yrange = ylim, marks = coor_df$SubType)

        # Find nearest neighbors
        nn <- nncross(points_ppp, points_ppp, k = 2:(k + 1))

        k_neighbors <- nn[, (k + 1):ncol(nn)]
        k_neighbors <- apply(k_neighbors, MARGIN = 2, function(x) {
            x <- as.numeric(x)
            return(sce_$CellID[x])
        })

        k_neighbors_ <- k_neighbors[TargetTypeidx, ]

        k_NeighborsList[[ROI]] <- as.numeric(k_neighbors_)
    }

    for (i in length(k_NeighborsList):1) {
        if (length(k_NeighborsList[[i]]) == 0) {
            k_NeighborsList[[i]] <- NULL
        }
    }

    if (return.cellID) {
        return(k_NeighborsList)
    }

    Subtypes <- names(table(colData(sce)[, NeighType]))
    InterDF <- as.data.frame(matrix(data = NA, nrow = length(k_NeighborsList), ncol = length(Subtypes)))
    colnames(InterDF) <- Subtypes

    for (i in 1:length(k_NeighborsList)) {
        CellIDTemp <- k_NeighborsList[[i]]
        CelltypeTemp <- colData(sce)[, NeighType][match(CellIDTemp, sce$CellID)]
        CelltypeTemp <- factor(CelltypeTemp, levels = Subtypes)
        CelltypeTemp <- as.data.frame(table(CelltypeTemp))
        InterDF[i, ] <- CelltypeTemp[, 2] / (length(k_NeighborsList[[i]]) / k)
    }

    rownames(InterDF) <- names(k_NeighborsList)
    InterDF$PID <- sapply(rownames(InterDF), function(x) {
        return(strsplit(x, "_")[[1]][1])
    })
    InterDF$RFS_status <- sce$RFS_status[match(InterDF$PID, sce$PID)]
    InterDF$RFS_time <- sce$RFS_time[match(InterDF$PID, sce$PID)]

    return(InterDF)
}

## Calcualte the neighors within d of certain cell subpopulation via spatstat
CalDNNBySpatstat <- function(sce, ROIs, CenType, NeighType = "SubType", d = 10, xlim = c(0, 1000), ylim = c(0, 1000)) {
    celltypes <- names(table(colData(sce)[, NeighType]))

    k_NeighborsList <- list()

    for (ROI in ROIs) {
        sce_ <- sce[, sce$ID == ROI]

        coor <- sce_$Position
        coor_df <- GetXandY(coor)
        coor_df$SubType <- colData(sce_)[, NeighType]

        # head(coor_df)
        TargetTypeidx <- (colData(sce_)[, NeighType] == CenType)
        if (nrow(coor_df[TargetTypeidx, ]) == 0) {
            next
        }

        # Convert your data frame to a ppp object
        points_ppp <- ppp(coor_df$cor_x, coor_df$cor_y, xrange = xlim, yrange = ylim, marks = coor_df$SubType)

        # Find nearest neighbors
        nn <- nncross(points_ppp, points_ppp, k = 2:21)

        k_neighborsDist <- nn[, 1:20]
        k_neighborsID <- nn[, 21:ncol(nn)]

        k_neighborsDist <- k_neighborsDist[TargetTypeidx, ]
        k_neighborsID <- k_neighborsID[TargetTypeidx, ]

        k_NeighborsDF <- as.data.frame(matrix(data = 0, nrow = nrow(k_neighborsID), ncol = length(celltypes)))
        colnames(k_NeighborsDF) <- celltypes

        for (i in 1:nrow(k_neighborsDist)) {
            idxTemp <- unname((k_neighborsDist[i, ] <= d))
            x <- unname(k_neighborsID[i, idxTemp])
            if (length(x) == 0) {
                next
            } else {
                dfTemp <- as.data.frame(table(colData(sce_)[as.numeric(x), NeighType]))
                k_NeighborsDF[i, match(dfTemp[, 1], colnames(k_NeighborsDF))] <- as.numeric(dfTemp[, 2])
            }
        }
        k_NeighborsDF$ID <- rep(ROI, nrow(k_NeighborsDF))
        rownames(k_NeighborsDF) <- sce_$CellID[TargetTypeidx]

        k_NeighborsList[[ROI]] <- k_NeighborsDF
    }

    return(k_NeighborsList)
}

## Specific cell neighbors' subpopulation DEGs analysis
SpeciNeiDEGAnalysis <- function(sce, IDcol, analysisType, metaMarkers, sample.size = 2000) {
    ReturnDF <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 4))

    for (cellType in analysisType) {
        sceType <- sce[, colData(sce)[, "SubType"] == cellType]

        IDx <- colData(sceType)[, IDcol]

        Assosce <- sceType[, IDx]
        NoAssosce <- sceType[, !IDx]

        Assoexp <- t(assay(Assosce))
        NoAssoexp <- t(assay(NoAssosce))

        set.seed(1)

        if (!is.null(sample.size)) {
            sample.idx <- sample.int(nrow(NoAssoexp), size = sample.size)
            NoAssoexp <- NoAssoexp[sample.idx, ]
        }

        Exp <- rbind(Assoexp, NoAssoexp)
        Exp <- as.data.frame(Exp)

        Exp <- Exp[, match(metaMarkers, colnames(Exp))]

        Exp$label <- c(rep(1, nrow(Assoexp)), rep(0, nrow(NoAssoexp)))

        mat_foldchangeMat <- FCandPvalueCal(mat = Exp, xCol = c(1, length(metaMarkers)), yCol = (length(metaMarkers) + 1), need.sample = FALSE)
        mat_foldchangeMat$Q.value <- p.adjust(mat_foldchangeMat$P.value, method = "BH")
        mat_foldchangeMat <- cbind("NeigborType" = rep(cellType, nrow(mat_foldchangeMat)), mat_foldchangeMat)

        ReturnDF <- rbind(ReturnDF, mat_foldchangeMat)
    }

    return(ReturnDF)
}

## 0-1 rescale
rescale <- function(x, new_range = c(0, 1)) {
    old_range <- range(x, na.rm = TRUE)
    x <- (x - old_range[1]) / (old_range[2] - old_range[1]) * diff(new_range) + new_range[1]
    return(x)
}

## plot delaunay triangulation network
plot_network_graph <- function(graph, coords, save_path, sample_fraction = 0.1) {
    # Get the largest connected component
    comp <- components(graph)
    largest_comp <- induced_subgraph(graph, which(comp$membership == which.max(comp$csize)))

    # Sample a subset of vertices from the largest connected component
    V_sub <- sample(V(largest_comp), size = floor(sample_fraction * vcount(largest_comp)))

    # Extract subgraph
    sub_graph <- induced_subgraph(largest_comp, vids = V_sub)
    sub_coords <- coords[as.numeric(names(V_sub)), ]

    # Set the layout for the subgraph
    layout <- layout_with_fr(sub_graph)

    # Define colors for different SubTypes
    unique_subtypes <- unique(sub_coords$SubType)
    color_map <- colorRampPalette(brewer.pal(min(9, length(unique_subtypes)), "Set1"))(length(unique_subtypes))
    names(color_map) <- unique_subtypes
    vertex_colors <- color_map[sub_coords$SubType]

    # Plot the subgraph
    pdf(save_path, height = 12, width = 16)

    plot(sub_graph,
        layout = layout,
        vertex.size = 2, # adjust as needed
        vertex.color = vertex_colors,
        vertex.label = NA, # no labels
        edge.color = "grey50",
        edge.width = 1, # adjust as needed
        main = "Network Visualization (Subset)"
    )

    dev.off()
}
