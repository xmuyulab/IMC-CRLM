# Functions for structural analysis
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(Rtsne)
library(RColorBrewer)
library(ggbeeswarm)

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


HeatmapForCelltypeInNeighbor <- function(sce, colname1, colname2, savePath) {
    ## transfer into heatmap matrix
    plotdf <- GetCelltype2NeighborMat(colData(sce), colname1, colname2)

    ## heatmap
    color <- colorRampPalette(c("#436eee", "white", "#EE0000"))(100)
    p <- pheatmap(plotdf,
        color = color, scale = "column",
        cluster_rows = F, cluster_cols = T,
        legend_labels = c("Abundance high", "Abundance low"), legend = T,
        show_rownames = T, show_colnames = T
    )

    pdf(paste0(savePath, "Cellular Neighbors celltype fraction heatmap.pdf"), width = 8, height = 6)
    print(p)
    dev.off()
    return(NULL)
}

## Comapre the Celluar pattern into groups
CompareCellularPattern <- function(sce, sep, countcol, n_cluster, savePath) {
    groups <- names(table(colData(sce)[, sep]))
    cat("The category is: ", groups, " in ", sep, "\n")

    sce1 <- sce[, colData(sce)[, sep] == groups[1]]
    sce2 <- sce[, colData(sce)[, sep] == groups[2]]

    ## ROI-level Boxplot
    abundance1 <- GetAbundance(sceobj = sce1, countcol = countcol, is.reuturnMeans = F)
    abundance2 <- GetAbundance(sceobj = sce2, countcol = countcol, is.reuturnMeans = F)

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
GetAbundance <- function(sceobj, countcol, is.fraction = TRUE, is.reuturnMeans = FALSE) {
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

    CellCountMat$Tissue <- cellMeta[match(rownames(CellCountMat), cellMeta$ID), ]$Tissue
    CellCountMat$RFS_status <- cellMeta[match(CellCountMat$PID, cellMeta$PID), ]$RFS_status
    CellCountMat$RFS_time <- cellMeta[match(CellCountMat$PID, cellMeta$PID), ]$RFS_time
    CellCountMat$KRAS_mutation <- cellMeta[match(CellCountMat$PID, cellMeta$PID), ]$KRAS_mutation

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
            expTemp <- Temp[, 1:(ncol(Temp) - 5)]
            MeanexpTemp <- apply(expTemp, MARGIN = 2, FUN = "mean")
            CellCountMat2[i, 1:(ncol(Temp) - 5)] <- MeanexpTemp
            CellCountMat2[i, (ncol(Temp) - 4):ncol(Temp)] <- Temp[1, (ncol(Temp) - 4):ncol(Temp)]
        }
        return(CellCountMat2)
    }
}

## Transform count matrix into ggplot plot matrix
TransformIntoPlotMat <- function(mat, valueCol) {
    exp <- mat[, valueCol]
    meta <- mat[, (max(valueCol) + 1):ncol(mat)]

    plotdf <- matrix(data = 0, nrow = nrow(exp) * ncol(exp), ncol = 5)
    colnames(plotdf) <- c("Abundance", "Pattern", "ROI", "Relapse", "RelapseTime")
    plotdf <- as.data.frame(plotdf)

    AbundanceVec <- as.numeric(as.matrix(exp))
    PatternVec <- rep(colnames(exp), each = nrow(exp))
    ROIVec <- rep(rownames(exp), times = ncol(exp))
    RelapsVec <- rep(meta$RFS_status, times = ncol(exp))
    RelapsTimeVec <- rep(meta$RFS_time, times = ncol(exp))

    plotdf$Abundance <- as.numeric(AbundanceVec)
    plotdf$Pattern <- as.factor(PatternVec)
    plotdf$ROI <- as.factor(ROIVec)
    plotdf$Relapse <- as.factor(RelapsVec)
    plotdf$RelapseTime <- as.numeric(RelapsTimeVec)

    return(plotdf)
}

## Boxplot For cellular pattern
BoxPlotForCellular <- function(mat1, mat2, sep, groups, valueCol, savePath) {
    plotdf1 <- TransformIntoPlotMat(mat1, valueCol)
    plotdf2 <- TransformIntoPlotMat(mat2, valueCol)

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
        stat_compare_means(aes(group = Group), label.y = 0.8, method = "t.test", label = "p.signif")

    pdf(paste0(savePath, "Cellular Neighborhood pattern difference of ", sep, ".pdf"), height = 6, width = 8)
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
Reclustering <- function(sce, markers, ReMajorType, ReclusterName, ncluster = 10, savePath) {
    ## extract major types
    sce_ <- sce[, colData(sce)$MajorType %in% ReMajorType]

    exp <- assay(sce_)
    exp <- exp[markers, ]

    ## K-means clustering
    exp <- t(exp) ## row should be sample
    set.seed(619)
    fit <- kmeans(exp, centers = ncluster, nstart = 25, iter.max = 500)
    table(fit$cluster)

    colData(sce_)[, ReclusterName] <- fit$cluster

    ## T-sne visualization
    sampleidx <- sample(1:nrow(exp), size = 15000, replace = F) ### sample 15k cells to visualize

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

    pdf(paste0(savePath, "tSNE reclustering of ", ReclusterName, ".pdf"), height = 6, width = 10)
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
    SubtypeInReclustering(sce_, reclusteringCol = ReclusterName, OrigintypeCol = "SubType", PatternCol = "kmeans_knn_20", savePath)

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
            panel.border = element_rect(colour = "white", fill = NA)
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

        pdf(paste0(savePath, marker, " expression on tSNE reclustering.pdf"), height = 5, width = 8)
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
PlotCelltypes <- function(sce, ROI, TypeCol, SavePath) {
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
FCandPvalueCal <- function(mat, xCol, yCol, need.sample = FALSE) {
    if (need.sample) {
        idx <- sample(1:nrow(mat), size = 2000, replace = FALSE)
        mat <- mat[idx, ]
    }
    groups <- names(table(mat[, yCol]))
    groups <- as.character(sort(as.numeric(groups), decreasing = F))
    if (length(groups) < 2) {
        return(0)
    }

    returnMat <- matrix(data = NA, nrow = length(xCol[1]:xCol[2]), ncol = 3)
    returnMat <- as.data.frame(returnMat)
    colnames(returnMat) <- c("Celltype", "Foldchange", "P.value")

    group1mat <- mat[which(mat[, yCol] == groups[1]), ]
    group2mat <- mat[which(mat[, yCol] == groups[2]), ]

    for (i in xCol[1]:xCol[2]) {
        typeTemp <- colnames(mat)[i]

        v1 <- group1mat[, i]
        v2 <- group2mat[, i]

        ## relaps versus no relaps
        foldchange <- mean(v2) / mean(v1)
        pvalue <- t.test(v2, v1)$p.value

        returnMat[i, ] <- c(typeTemp, foldchange, pvalue)
    }

    return(returnMat)
}

## Plot volcano plot
VolcanoPlot <- function(df, pthreshold = 0.05, fcthreshold = 1.4, feature, filename = NULL) {
    df$Foldchange <- as.numeric(df$Foldchange)
    df$P.value <- as.numeric(df$P.value)

    df$change <- as.factor(ifelse(df$P.value < pthreshold & abs(log2(df$Foldchange)) > log2(fcthreshold),
        ifelse(log2(df$Foldchange) > log2(fcthreshold), "Up-regulate", "Down-regulate"), "Non-significant"
    ))

    # 样本标签
    df$label <- ifelse(df[, 3] < pthreshold & abs(log2(df$Foldchange)) > log2(fcthreshold), as.character(df[, 1]), "")

    # 绘制火山图
    p.vol <- ggplot(
        data = df,
        aes(x = log2(Foldchange), y = -log10(P.value), colour = change, fill = change)
    ) +
        scale_color_manual(values = c("Down-regulate" = "blue", "Non-significant" = "grey", "Up-regulate" = "red")) +
        geom_point(alpha = 0.4, size = 3.5) +
        # 标签
        geom_text_repel(aes(x = log2(Foldchange), y = -log10(P.value), label = label),
            size = 3,
            box.padding = unit(0.6, "lines"), point.padding = unit(0.7, "lines"),
            segment.color = "black", show.legend = FALSE
        ) +
        # 辅助线
        geom_vline(xintercept = c(-(log2(fcthreshold)), (log2(fcthreshold))), lty = 4, col = "black", lwd = 0.8) +
        geom_hline(yintercept = -log10(pthreshold), lty = 4, col = "black", lwd = 0.8) +
        theme_bw() +
        labs(x = "log2(Fold Change)", y = "-log10(P value)", title = paste0("Volcano Plot of Different Expression Markers in ", feature)) +
        # 坐标轴标题、标签和图例相关设置
        theme(
            axis.text = element_text(size = 11), axis.title = element_text(size = 13), # 坐标轴标签和标题
            plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), # 标题
            legend.text = element_text(size = 11), legend.title = element_text(size = 13), # 图例标签和标题
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
        ) # 图边距

    ggsave(p.vol, filename = filename)

    return(NULL)
}

## Assign new label
AssignNewlabel <- function(sce_, phenoLabel, ReclusterName, interstType, cutoffType) {
    colData(sce_)[, phenoLabel] <- "Pheno_neg"
    colData(sce_)[, phenoLabel][which(colData(sce_)[, ReclusterName] %in% interstType)] <- "Pheno_pos"

    CountMat <- GetAbundance(sce_, countcol = phenoLabel, is.fraction = FALSE, is.reuturnMeans = FALSE)
    if (cutoffType == "median") {
        Cutoff <- median(CountMat[, "Pheno_pos"])
    }
    if (cutoffType == "mean") {
        Cutoff <- mean(CountMat[, "Pheno_pos"])
    }

    RList <- list()
    RList[["sce"]] <- sce_
    RList[["high_ROI"]] <- rownames(CountMat[which(CountMat[, "Pheno_pos"] >= Cutoff), ])
    RList[["low_ROI"]] <- rownames(CountMat[which(CountMat[, "Pheno_pos"] < Cutoff), ])

    return(RList)
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
        rowname_ <- rownames(csvTemp)

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
        rowname_ <- rownames(csvTemp)

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

## Swarm plot for high- and low- marker groups certain celltype
AbundanceSwarmPlot <- function(AbundanceDF1, AbundanceDF2, groupsName, celltype, marker, savePath) {
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
        geom_boxplot(show.legend = FALSE) +
        geom_beeswarm(dodge.width = 0.8, shape = 21) +
        theme_bw() +
        labs(x = NULL) +
        scale_color_manual(values = c("#648fff", "#d36b1c")) +
        stat_compare_means(label.y = 0.7, method = "t.test")

    pdf(paste0(savePath, "Swarmplot of ", celltype, " in ", marker, " group.pdf"), height = 3, width = 4)
    print(p)
    dev.off()

    return(NULL)
}

## Survival analysis for Phenotype associated label
SurvivalForPhenoAssoLabel <- function(AbundanceDF1, AbundanceDF2, time, status, marker, savePath) {
    AbundanceDF1 <- AbundanceDF1[, c(time, status)]
    AbundanceDF2 <- AbundanceDF2[, c(time, status)]
    GroupVec <- c(rep("High", nrow(AbundanceDF1)), rep("Low", nrow(AbundanceDF2)))

    plotdf <- rbind(AbundanceDF1, AbundanceDF2)
    plotdf <- cbind(plotdf, GroupVec)
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

    pdf(paste0(savePath, "Suvival analysis for phenotype-associated ROI of ", marker, ".pdf"), height = 6, width = 8)
    print(p)
    dev.off()

    return(NULL)
}
