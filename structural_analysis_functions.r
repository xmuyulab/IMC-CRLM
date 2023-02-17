# Functions for structural analysis
library(pheatmap)
library(ggplot2)

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
MergeAbundanceResult <- function(sce) {
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
