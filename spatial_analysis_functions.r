## functions for permutation test in R
library(pheatmap)
library(ComplexHeatmap)
library(ggDoubleHeat)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(ggprism)
library(grid)
library(patchwork)

library(SingleCellExperiment)
library(survival)
library(forcats)

library(rtiff)
library(tiff)
library(magick)

source("./structural_analysis_functions.r")

### Make groupInfo file
GetGroupInfo <- function(sce, clinical) {
    ROIs <- names(table(sce$ID))

    returnDF <- matrix(data = NA, nrow = length(ROIs), ncol = 2)
    colnames(returnDF) <- c("RFS_time", "RFS_status")
    rownames(returnDF) <- ROIs

    for (i in 1:nrow(returnDF)) {
        PIDTemp <- strsplit(rownames(returnDF)[i], "_")[[1]][1]
        returnDF[i, 1] <- clinical[clinical$PID == PIDTemp, "RFS_time"]
        returnDF[i, 2] <- clinical[clinical$PID == PIDTemp, "RFS_status"]
    }

    return(as.data.frame(returnDF))
}

### get the results of permutation test
getResult <- function(ResultPath, GroupInfo, clinicalGroup, celltypes, p_threshold = 0.05) {
    ## Interaction
    numtypes <- length(celltypes)

    interDF <- matrix(data = 0, nrow = numtypes, ncol = numtypes)
    interDF <- as.data.frame(interDF)
    rownames(interDF) <- celltypes
    colnames(interDF) <- celltypes

    IDsTemp <- rownames(GroupInfo[GroupInfo$RFS_status %in% clinicalGroup, ])
    numID <- length(IDsTemp)

    for (ID in IDsTemp) {
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

    IDsTemp <- rownames(GroupInfo[GroupInfo$RFS_status == clinicalGroup, ])

    for (ID in IDsTemp) {
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

### Transfrom matrix into plot dataframe
BubblePlot <- function(MergeDF, LabelDF, savePath) {
    if (class(MergeDF) == "list") {
        MergeDF <- MergeDF[[1]]
        LabelDF <- LabelDF[[1]]
    }

    nrow_ <- nrow(MergeDF)
    ncol_ <- ncol(MergeDF)

    plotdf <- matrix(data = NA, nrow = nrow_ * ncol_, ncol = 4)
    plotdf <- as.data.frame(plotdf)

    c1 <- rep(rownames(MergeDF), each = nrow_)
    c2 <- rep(colnames(MergeDF), time = ncol_)
    num <- c()
    label <- c()

    for (i in 1:nrow_) {
        num <- c(num, as.numeric(MergeDF[i, ]))
        label <- c(label, as.numeric(LabelDF[i, ]))
    }
    label <- ifelse(label == 1, "Interaction", "Avoidence")

    colnames(plotdf) <- c("Celltype1", "Celltype2", "ROI Fraction", "Spatial Status")

    plotdf["Celltype1"] <- c1
    plotdf["Celltype2"] <- c2
    plotdf["ROI Fraction"] <- num
    plotdf["Spatial Status"] <- as.factor(label)

    ### bubble plot
    p <- ggplot(plotdf, aes(x = Celltype1, y = Celltype2, size = `ROI Fraction`, color = `Spatial Status`)) +
        geom_point() +
        theme(
            panel.background = element_blank(),
            panel.grid.major = element_line(colour = "white"),
            panel.border = element_rect(colour = "white", fill = NA),
            axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.25)
        )

    pdf(savePath, width = 8, height = 8)
    print(p)
    dev.off()

    return(NULL)
}

# DoubleHeat function to create spatial aggregation cell interaction data
SingleHeatmap <- function(data1, label1, savePath) {
    # Multiply the data frames element-wise
    df1 <- data1 * label1

    # Create a symmetric sequence of breaks with 0 in the middle
    breaks <- seq(from = -1, to = 1, length.out = 101)

    # Create custom color palette with blue for negative values, white for zero, and red for positive values
    color_palette <- colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1)

    # Check the specified plot type and create the corresponding plot
    p <- pheatmap(
        df1,
        scale = "none",
        cellwidth = 16, cellheight = 12,
        cluster_row = F, cluster_col = F,
        angle_col = "90",
        breaks = breaks,
        color = color_palette,
        legend = TRUE,
        legend_breaks = c(min(breaks), max(breaks))
    )

    # Save the plot as a PDF
    pdf(savePath, height = 6, width = 8)
    print(p)
    dev.off()

    return(NULL)
}

# Heatmap to visualize of cell interaction data
DoubleHeat <- function(data1, label1, group1, data2, label2, group2, plot = "circle", savePath) {
    # Multiply the data frames element-wise
    df1 <- data1 * label1
    df2 <- data2 * label2

    # Initialize an empty data frame
    plotdf <- matrix(data = 0, nrow = nrow(df1) * ncol(df1), ncol = 4)
    plotdf <- as.data.frame(plotdf)

    # Fill the data frame with the calculated values
    plotdf[, 1] <- rep(rownames(df1), times = ncol(df1))
    plotdf[, 2] <- rep(colnames(df1), each = nrow(df1))
    plotdf[, 3] <- as.numeric(as.matrix(df1))
    plotdf[, 4] <- as.numeric(as.matrix(df2))

    # Set the column names of the data frame
    colnames(plotdf) <- c("Celltype1", "Celltype2", "Interaction1", "Interaction2")

    # Check the specified plot type and create the corresponding plot
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

    # Save the plot as a PDF
    pdf(savePath, height = 6, width = 8)
    print(p)
    dev.off()

    return(NULL)
}

## t-test for interaction number between groups
LoadAnalysisResult <- function(IDs1, celltypes, sce) {
    numtypes <- length(celltypes)

    array1 <- array(data = NA, dim = c(numtypes, numtypes, length(IDs1)))
    rownames(array1) <- celltypes
    colnames(array1) <- celltypes

    i <- 1
    for (ID in IDs1) {
        sceTemp <- sce[, sce$ID == ID]

        ## load result dataframe
        filepathTemp <- paste0(ResultPath, ID, "/TureInteraction/TureInteraction.csv")
        csvTemp <- read.csv(filepathTemp)
        rownames(csvTemp) <- csvTemp[, 1]
        csvTemp <- csvTemp[, -1]
        colnames(csvTemp) <- rownames(csvTemp)
        csvTemp <- as.matrix(csvTemp)

        ## load mask dataframe
        maskfilepathTemp <- paste0(ResultPath, ID, "/ClosePvalue/ClosePvalue.csv")
        maskTemp <- read.csv(maskfilepathTemp)
        maskTemp <- maskTemp[, -1]
        maskTemp <- as.matrix(maskTemp)

        maskTemp <- ifelse(maskTemp <= 0.05, 1, 0)
        csvTemp <- csvTemp * maskTemp

        ## sort the interaction number
        if (nrow(csvTemp) != numtypes) {
            TempMat <- matrix(data = 0, nrow = numtypes, ncol = numtypes)
            colnames(TempMat) <- celltypes
            rownames(TempMat) <- celltypes

            ## match the names of csvTemp and TempMat
            for (j in 1:ncol(csvTemp)) {
                colName <- colnames(csvTemp)[j]
                TempMat[match(rownames(csvTemp), rownames(TempMat)), colName] <- csvTemp[, j]
            }
        } else {
            csvTemp <- csvTemp[match(celltypes, rownames(csvTemp)), ]
            TempMat <- csvTemp[, match(celltypes, colnames(csvTemp))]
        }

        ## calculate the cell number of certain cell types
        cellnumMat1 <- matrix(data = 0, nrow = numtypes, ncol = numtypes)
        rownames(cellnumMat1) <- celltypes
        colnames(cellnumMat1) <- celltypes

        cellnumMat2 <- cellnumMat1

        for (j in 1:ncol(cellnumMat1)) {
            c1 <- colnames(cellnumMat1)[j]
            c1num <- ncol(sceTemp[, sceTemp$Cell_Subtype == c1])
            cellnumMat1[, j] <- cellnumMat1[, j] + c1num
        }
        for (j in 1:nrow(cellnumMat2)) {
            c2 <- rownames(cellnumMat2)[j]
            c2num <- ncol(sceTemp[, sceTemp$Cell_Subtype == c2])
            cellnumMat2[j, ] <- cellnumMat2[j, ] + c2num
        }
        cellnumMat <- sqrt(cellnumMat1 ** 2 + cellnumMat2 ** 2) 

        TempMat <- ifelse(cellnumMat == 0, 0, round(TempMat / cellnumMat, 6))

        array1[, , i] <- TempMat
        i <- i + 1
    }

    return(array1)
}

TestDiff <- function(array1, array2, celltypes, savepath, do.fdr = TRUE, FC_threshold = 1.5, Sig = 0.05, returnMat = F) {
    numtypes <- length(celltypes)
    PvalueMat <- matrix(data = NA, nrow = numtypes, ncol = numtypes)
    rownames(PvalueMat) <- celltypes
    colnames(PvalueMat) <- celltypes

    MaskMat <- PvalueMat
    FoldChangeMat <- PvalueMat

    for (i in 1:nrow(PvalueMat)) {
        for (j in 1:ncol(PvalueMat)) {
            value1 <- as.numeric(na.omit(array1[i, j, ]))
            value2 <- as.numeric(na.omit(array2[i, j, ]))

            ## Mask
            if (mean(value1) < mean(value2)) {
                MaskMat[i, j] <- 1
            } else {
                MaskMat[i, j] <- -1
            }

            ## Fold Change
            if (mean(value2) == 0) {
                FoldChangeMat[i, j] <- 0
            } else {
                FoldChangeMat[i, j] <- mean(value1) / mean(value2)
            }

            ## P value
            if (mean(value2) == 0) {
                PvalueMat[i, j] <- 1
            } else {
                PvalueMat[i, j] <- t.test(value1, value2)$p.value
            }
        }
    }
    if (do.fdr) {
        qvalueMat <- matrix(p.adjust(as.numeric(PvalueMat), method = "BH"), nrow = nrow(PvalueMat))
        colnames(qvalueMat) <- colnames(PvalueMat)
        rownames(qvalueMat) <- rownames(PvalueMat)
        qvalueMat <- -log10(qvalueMat)
        FoldChangeMat <- ifelse(FoldChangeMat == 0, 0, log2(FoldChangeMat))
        FoldChangeMat <- ifelse(abs(FoldChangeMat) < log2(FC_threshold), 0, FoldChangeMat)

        if (returnMat) {
            list_ <- list()
            list_[["qvalueMat"]] <- qvalueMat
            list_[["FoldChangeMat"]] <- FoldChangeMat
            return(list_)
        }
        HeatmapForDiff(qvalueMat, MaskMat, FoldChangeMat, Sig = Sig, savepath)
    } else {
        PvalueMat <- -log10(PvalueMat)
        FoldChangeMat <- ifelse(FoldChangeMat == 0, 0, log2(FoldChangeMat))
        FoldChangeMat <- ifelse(abs(FoldChangeMat) < log2(FC_threshold), 0, FoldChangeMat)

        if (returnMat) {
            list_ <- list()
            list_[["PvalueMat"]] <- PvalueMat
            list_[["FoldChangeMat"]] <- FoldChangeMat
            return(list_)
        }
        HeatmapForDiff(PvalueMat, MaskMat, FoldChangeMat, Sig = Sig, savepath)
    }

    return(NULL)
}

getSig.1 <- function(dc) {
    sc <- ""

    if (dc >= 2) {
        sc <- "***"
    } else if (dc >= 1.3) {
        sc <- "**"
    } else if (dc >= 1) {
        sc <- "*"
    }
    return(sc)
}
getSig.05 <- function(dc) {
    sc <- ""

    if (dc >= 3) {
        sc <- "***"
    } else if (dc >= 2) {
        sc <- "**"
    } else if (dc >= 1.3) {
        sc <- "*"
    }
    return(sc)
}
getSig.01 <- function(dc) {
    sc <- ""
    if (dc >= 4) {
        sc <- "***"
    } else if (dc >= 3) {
        sc <- "**"
    } else if (dc >= 2) {
        sc <- "*"
    }
    return(sc)
}
## Modify
HeatmapForDiff <- function(PvalueMat, MaskMat, FoldChangeMat, Sig = 0.05, savepath) {
    ## Set the FC less than cutoff as n.s.
    TempDF <- ifelse(FoldChangeMat == 0, 0, 1)
    PvalueMat <- PvalueMat * TempDF

    if (Sig == 0.05) {
        sig_mat <- matrix(sapply(PvalueMat, getSig.05), nrow = nrow(PvalueMat))
    }
    if (Sig == 0.01) {
        sig_mat <- matrix(sapply(PvalueMat, getSig.01), nrow = nrow(PvalueMat))
    }
    if (Sig == 0.1) {
        sig_mat <- matrix(sapply(PvalueMat, getSig.1), nrow = nrow(PvalueMat))
    }
    plotdf <- FoldChangeMat #* MaskMat

    # Determine the range of fold changes in your data
    min_value <- min(plotdf)
    max_value <- max(plotdf)

    # Determine the absolute maximum value, either positive or negative
    abs_max <- max(abs(min_value), abs(max_value))

    # Create a symmetric sequence of breaks with 0 in the middle
    breaks <- seq(from = -abs_max, to = abs_max, length.out = 101)

    # Create custom color palette with blue for negative values, white for zero, and red for positive values
    color_palette <- colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1)

    # Create custom labels for the color legend
    legend_labels <- c("Non-Relapse", "Relapse")

    # Create the heatmap with custom breaks, color palette, and legend labels
    p <- pheatmap(
        plotdf,
        cellwidth = 16, cellheight = 12,
        cluster_row = F, cluster_col = F,
        angle_col = "90", display_numbers = sig_mat, fontsize_number = 15,
        breaks = breaks,
        color = color_palette,
        legend = TRUE,
        legend_labels = legend_labels,
        legend_breaks = c(min(breaks), max(breaks))
    )

    pdf(savepath, width = 10, height = 8)
    print(p)
    dev.off()

    return(NULL)
}

getInteracDiff <- function(ResultPath, sce, GroupInfo = NULL, groups = NULL, celltypes, savepath, do.fdr = TRUE, FC_threshold = 1.5, IDs1 = NULL, IDs2 = NULL, Sig = 0.05) {
    if (is.null(IDs1)) {
        IDs1 <- rownames(GroupInfo[GroupInfo$RFS_status == groups[1], ])
        IDs2 <- rownames(GroupInfo[GroupInfo$RFS_status == groups[2], ])
    }

    array1 <- LoadAnalysisResult(IDs1, celltypes, sce)
    array2 <- LoadAnalysisResult(IDs2, celltypes, sce)

    TestDiff(array1, array2, celltypes, savepath, do.fdr = do.fdr, FC_threshold = FC_threshold, Sig = Sig)

    return(NULL)
}

## Invertion Barplot to visualize interaction difference
InterDiffInvertBarplot <- function(ResultPath, sce, GroupInfo = NULL, groupCol = NULL, groups = NULL, PlotTypes, savepath, IDs1 = NULL, IDs2 = NULL) {
    if (is.null(IDs1)) {
        IDs1 <- rownames(GroupInfo[GroupInfo[, groupCol] == groups[1], ])
        IDs2 <- rownames(GroupInfo[GroupInfo[, groupCol] == groups[2], ])
    }
    celltypes <- names(table(sce$SubType))
    array1 <- LoadAnalysisResult(IDs1, celltypes, sce)
    array2 <- LoadAnalysisResult(IDs2, celltypes, sce)

    list_ <- TestDiff(array1, array2, celltypes, savepath = savepath, returnMat = T)
    qvalueMat <- list_[["qvalueMat"]]
    FoldChangeMat <- list_[["FoldChangeMat"]]

    if (!is.null(PlotTypes)) {
        idx <- match(PlotTypes, rownames(qvalueMat))
    } else {
        idx <- match(celltypes, rownames(qvalueMat))
    }

    qvalueMat <- qvalueMat[idx, idx]
    FoldChangeMat <- FoldChangeMat[idx, idx]

    pairs <- rep(sapply(rownames(qvalueMat), function(x) {
        paste0(x, "-", colnames(qvalueMat))
    }))
    qvalues <- as.numeric(qvalueMat)
    foldchangevalues <- as.numeric(FoldChangeMat)
    df <- as.data.frame(cbind(pairs, qvalues, foldchangevalues))
    colnames(df) <- c("CelltypePair", "LgQValue", "Log2FoldChange")
    df$LgQValue <- as.numeric(df$LgQValue)
    df$Log2FoldChange <- as.numeric(df$Log2FoldChange)

    # Rearrange the data frame based on Log2FoldChange
    df <- df[order(df$Log2FoldChange), ]

    # Reorder factor levels of CelltypePair
    df$CelltypePair <- factor(df$CelltypePair, levels = df$CelltypePair)

    ## omit the same sample
    df <- df[seq.int(1, nrow(df), by = 2), ]

    # Create a horizontal bar plot with LgQValue as the transparency of each bar
    p <- ggplot(df, aes(x = CelltypePair, y = Log2FoldChange, fill = Log2FoldChange > 0, alpha = LgQValue)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = c("steelblue", "darkorange"), labels = c("Negative", "Positive"), name = "Direction") +
        scale_alpha_continuous(range = c(0.4, 1)) +
        labs(x = "Cell Type Pair", y = "log2(Fold Change)", title = "Horizontal Bar Plot of Fold Changes") +
        theme_minimal() +
        theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 12, face = "bold"),
            axis.text = element_text(size = 10),
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10),
            legend.position = "bottom"
        )

    return(p)
}

## Invertion Barplot to visualize interaction difference
InterAbundanceBarplot <- function(ResultPath, sce, GroupInfo = NULL, groupCol = NULL, groups = NULL, groupsName, PlotTypes, savepath = NULL, IDs1 = NULL, IDs2 = NULL) {
    if (is.null(IDs1)) {
        IDs1 <- rownames(GroupInfo[GroupInfo[, groupCol] == groups[1], ])
        IDs2 <- rownames(GroupInfo[GroupInfo[, groupCol] == groups[2], ])
    }
    celltypes <- names(table(sce$SubType))
    ## Get interaction
    array1 <- LoadAnalysisResult(IDs1, celltypes, sce)
    array2 <- LoadAnalysisResult(IDs2, celltypes, sce)

    list_ <- TestDiff(array1, array2, celltypes, savepath = savepath, returnMat = T)
    qvalueMat <- list_[["qvalueMat"]]
    FoldChangeMat <- list_[["FoldChangeMat"]]

    ## Get abundance
    sce1 <- sce[, sce$ID %in% IDs1]
    sce2 <- sce[, sce$ID %in% IDs2]

    abundanceMat1 <- GetAbundance(sce1, countcol = "SubType", is.fraction = T)
    abundanceMat2 <- GetAbundance(sce2, countcol = "SubType", is.fraction = T)

    if (is.null(PlotTypes)) {
        PlotTypes <- celltypes
    }

    idx_ <- match(PlotTypes, colnames(abundanceMat1))
    abundanceMat1 <- abundanceMat1[, idx_]
    idx_ <- match(PlotTypes, colnames(abundanceMat2))
    abundanceMat2 <- abundanceMat2[, idx_]

    idx <- match(PlotTypes, rownames(qvalueMat))
    qvalueMat <- qvalueMat[idx, idx]
    FoldChangeMat <- FoldChangeMat[idx, idx]

    pairs <- c()
    cor_coef <- c()
    cor_test <- c()

    for (i in 1:nrow(qvalueMat)) {
        for (j in 1:ncol(qvalueMat)) {
            c1 <- rownames(qvalueMat)[i]
            c2 <- colnames(qvalueMat)[j]
            pairs <- c(pairs, paste0(c1, "-", c2))
            sign_ <- ifelse(FoldChangeMat[i, j] > 0, T, F)
            if (sign_) {
                vTemp <- cbind(abundanceMat1[, c1], abundanceMat1[, c2])
                coefTemp <- cor(vTemp[, 1], vTemp[, 2], method = "spearman")
                testTemp <- cor.test(vTemp[, 1], vTemp[, 2], method = "spearman")$p.value
                cor_coef <- c(cor_coef, coefTemp)
                cor_test <- c(cor_test, testTemp)
            } else {
                vTemp <- cbind(abundanceMat2[, c1], abundanceMat2[, c2])
                coefTemp <- cor(vTemp[, 1], vTemp[, 2], method = "spearman")
                testTemp <- cor.test(vTemp[, 1], vTemp[, 2], method = "spearman")$p.value
                cor_coef <- c(cor_coef, coefTemp)
                cor_test <- c(cor_test, testTemp)
            }
        }
    }

    qvalues <- as.numeric(qvalueMat)
    foldchangevalues <- as.numeric(FoldChangeMat)

    df <- as.data.frame(cbind(pairs, qvalues, foldchangevalues, cor_coef, cor_test))
    colnames(df) <- c("CelltypePair", "LgQValue", "Log2FoldChange", "Cor_coef", "Cor_test")
    df$LgQValue <- as.numeric(df$LgQValue)
    df$Log2FoldChange <- as.numeric(df$Log2FoldChange)
    df$Cor_coef <- as.numeric(df$Cor_coef)
    df$Cor_test <- as.numeric(df$Cor_test)
    df$LgCor_test <- -log10(df$Cor_test + 0.0001)

    # Rearrange the data frame based on Log2FoldChange
    df <- df[order(df$Log2FoldChange), ]

    # Reorder factor levels of CelltypePair
    df$CelltypePair <- factor(df$CelltypePair, levels = df$CelltypePair)

    ## omit the same sample
    df <- df[seq.int(1, nrow(df), by = 2), ]

    # Create a horizontal bar plot with LgQValue as the transparency of each bar
    p1 <- ggplot(df, aes(x = CelltypePair, y = -abs(Log2FoldChange), fill = Log2FoldChange > 0, alpha = LgQValue)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = c("steelblue", "darkorange"), labels = c(groupsName[2], groupsName[1]), name = "Fold Change") +
        scale_alpha_continuous(range = c(0.4, 1), name = "Q Value") +
        labs(x = "Cell Type Pair", y = NULL) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(size = 10),
            axis.text.y = element_blank(),
            legend.title = element_text(size = 10, face = "bold"),
            legend.text = element_text(size = 8),
            legend.position = "left"
        )

    # Create p2 with y-axis text
    p2 <- ggplot(df, aes(x = CelltypePair, y = abs(Cor_coef), fill = Cor_coef > 0, alpha = LgCor_test)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = c("steelblue", "darkorange"), labels = c("Negative", "Positive"), name = "Correlation") +
        scale_alpha_continuous(range = c(0.4, 1), name = "Cor P-Value") +
        labs(x = NULL, y = NULL) +
        theme_minimal() +
        theme(
            axis.text = element_text(size = 10),
            legend.title = element_text(size = 10, face = "bold"),
            legend.text = element_text(size = 8),
            legend.position = "right"
        )

    # Combine the plots using plot_grid
    p <- plot_grid(p1, p2, ncol = 2, align = "hv", axis = "lrtb")

    return(p)
}

## Stack plot to visualize certain interaction difference
InterDiffStackplot <- function(ResultPath, sce, GroupInfo = NULL, groupCol, groups, groupsName, celltypespair, savepath, IDs1 = NULL, IDs2 = NULL) {
    if (is.null(IDs1)) {
        IDs1 <- rownames(GroupInfo[GroupInfo[, groupCol] == groups[1], ])
        IDs2 <- rownames(GroupInfo[GroupInfo[, groupCol] == groups[2], ])
    }
    celltypes <- names(table(sce$SubType))
    array1 <- LoadAnalysisResult(IDs1, celltypes, sce)
    array2 <- LoadAnalysisResult(IDs2, celltypes, sce)

    numPair <- nrow(celltypespair)

    # Prepare data for the plot
    df <- as.data.frame(matrix(data = NA, ncol = 3, nrow = 0))
    len <- length(c(IDs1, IDs2))

    for (i in 1:numPair) {
        c1 <- celltypespair[i, 1]
        c2 <- celltypespair[i, 2]
        row_ <- match(c1, rownames(array1[, , 1]))
        col_ <- match(c2, colnames(array1[, , 1]))

        v1 <- array1[row_, col_, ]
        v2 <- array2[row_, col_, ]

        valueVec <- c(v1, v2)
        groupVec <- c(rep(groups[1], length(v1)), rep(groups[2], length(v2)))
        interVec <- rep(paste0(c1, ":", c2), len)

        dfTemp <- cbind(valueVec, groupVec, interVec)
        df <- rbind(df, dfTemp)
    }
    colnames(df) <- c("Strength", "ClinicalGroup", "Interaction")
    df$ClinicalGroup <- ifelse(df$ClinicalGroup == 1, groupsName[1], groupsName[2])
    df[, 1] <- as.numeric(df[, 1])
    df[, 2] <- as.character(df[, 2])
    df[, 3] <- as.character(df[, 3])

    # Create the plot
    p <- ggbarplot(df,
        x = "ClinicalGroup",
        y = "Strength",
        fill = "Interaction",
        facet.by = "Interaction",
        color = "black",
        palette = ggsci::pal_jco("default")(10),
        add = "mean_sd",
        xlab = "Clinical Group",
        ylab = "Strength",
        legend = "none",
        ggtheme = theme_minimal()
    ) +
        theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 12, face = "bold"),
            axis.text = element_text(size = 10),
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10),
            legend.key.size = unit(1, "cm"),
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        )
    # Add the p-value label above the middle of the line
    ymax <- 0.4
    p <- p +
        stat_compare_means(aes(group = ClinicalGroup),
            method = "t.test",
            label.y = ymax + 0.05 * (max(df$Strength, na.rm = TRUE) - min(df$Strength, na.rm = TRUE))
        )
    pdf(savepath, height = 6, width = 8)
    print(p)
    dev.off()

    return(NULL)
}

## Interaction circle plot
InterCircleplot <- function(ResultPath, sce, GroupInfo = NULL, groupCol, groups, groupsName, PlotTypes, savepath, IDs1 = NULL, IDs2 = NULL) {
    if (is.null(IDs1)) {
        IDs1 <- rownames(GroupInfo[GroupInfo[, groupCol] == groups[1], ])
        IDs2 <- rownames(GroupInfo[GroupInfo[, groupCol] == groups[2], ])
    }
    celltypes <- names(table(sce$SubType))
    array1 <- LoadAnalysisResult(IDs1, celltypes, sce)
    array2 <- LoadAnalysisResult(IDs2, celltypes, sce)

    idx <- match(PlotTypes, rownames(array1))

    array1_ <- matrix(data = NA, nrow = nrow(array1), ncol = ncol(array1))
    array2_ <- array1_
    for (i in 1:nrow(array1_)) {
        for (j in 1:ncol(array1_)) {
            array1_[i, j] <- mean(array1[i, j, ])
        }
    }
    for (i in 1:nrow(array2_)) {
        for (j in 1:ncol(array2_)) {
            array2_[i, j] <- mean(array2[i, j, ])
        }
    }
    rownames(array1_) <- rownames(array1)
    colnames(array1_) <- colnames(array1)
    rownames(array2_) <- rownames(array2)
    colnames(array2_) <- colnames(array2)

    array1_ <- array1_[idx, idx]
    array2_ <- array2_[idx, idx]

    cell_group <- c(rep(rep(colnames(array1_), each = nrow(array1_)), times = 2))
    interacting_partner <- c(rep(rep(rownames(array1_), times = ncol(array1_)), times = 2))
    interaction_strength <- c(as.numeric(array1_), as.numeric(array2_))
    group <- c(rep(groupsName[1], length(array1_)), rep(groupsName[2], length(array2_)))

    cell_data <- cbind(cell_group, interacting_partner, interaction_strength, group)
    cell_data <- as.data.frame(cell_data)

    # Filter out interactions with strength = 0
    cell_data <- cell_data[cell_data$interaction_strength > 0, ]

    # Create an igraph object
    cell_network <- graph_from_data_frame(cell_data, directed = FALSE)

    # Create a circular plot
    p <- ggraph(cell_network, layout = "circle") +
        geom_edge_link(aes(width = interaction_strength, alpha = interaction_strength), color = "steelblue") +
        geom_node_point(size = 5, color = "darkorange") +
        geom_node_text(aes(label = name), repel = TRUE, color = "black", size = 4) +
        theme_void() +
        theme(legend.position = "none") +
        labs(title = "Cell-cell Interaction Circular Plot") +
        facet_wrap(~group)

    return(p)
}

## cox test
cox_test <- function(ResultPath, GroupInfo, groups, celltypes, savepath) {
    IDs1 <- rownames(GroupInfo[GroupInfo$RFS_status == groups[1], ])
    IDs2 <- rownames(GroupInfo[GroupInfo$RFS_status == groups[2], ])

    array1 <- LoadAnalysisResult(IDs1, celltypes)
    array2 <- LoadAnalysisResult(IDs2, celltypes)

    GroupInfo <- GroupInfo[c(IDs1, IDs2), ]

    time <- GroupInfo$RFS_time
    status <- GroupInfo$RFS_status
    value <- c(array1[1, 1, ], array2[1, 1, ])

    data <- matrix(data = NA, nrow = length(time), ncol = 3)
    data <- as.data.frame(data)
    colnames(data) <- c("Number", "RFS_time", "RFS_status")
    data$Number <- value
    data$RFS_time <- as.numeric(time)
    data$RFS_status <- ifelse(status == "Recurrence", 1, 0)
    data <- na.omit(data)

    data$Number[1:61] <- sample(1:100, size = 61)
    data$Number[62:115] <- sample(100:200, size = 54)

    fit <- coxph(Surv(RFS_time, RFS_status) ~ Number, data = data)
    summary(fit)

    return(NULL)
}

## visualize functions for sce object

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

## plot Marker expression value on cell-level
PlotMarker <- function(sce, ROI, Marker, SavePath) {
    if (!Marker %in% rownames(sce)) {
        cat(Marker, " is no in sce object!", "\n")
        return(NULL)
    }

    sce_ <- sce[, sce$ID == ROI]

    coorList <- getCoordinate(sce_)

    value <- sce_[Marker, ]
    value <- as.vector(assay(value))

    plotdf <- as.data.frame(matrix(data = NA, nrow = length(value), ncol = 3))
    colnames(plotdf) <- c("x", "y", "Expression")
    plotdf["x"] <- coorList[[1]]
    plotdf["y"] <- coorList[[2]]
    plotdf["Expression"] <- value

    custom_color_scale <- colorRampPalette(c("lightgrey", "red"), alpha = TRUE)(10)
    custom_color_scale <- alpha(custom_color_scale, ifelse(seq(0.0, 1, length.out = 10) < 0.4, 0.3, 1))

    p <- ggplot(plotdf, aes(x = x, y = y, color = Expression)) +
        geom_point(size = 1) +
        scale_colour_gradientn(colours = custom_color_scale, limits = c(0, 1)) +
        labs(title = paste0(Marker, " expression level in ", ROI)) +
        theme_test()

    pdf(SavePath, height = 6, width = 8)
    print(p)
    dev.off()

    return(p)
}

## plot celltypes on cell-level
PlotCelltypes <- function(sce, ROI, selectCelltypes, whichType = "SubType", SavePath = NULL) {
    # colnames(colData(sce))
    sce_ <- sce[, sce$ID == ROI]

    coorList <- getCoordinate(sce_)

    celltypes <- colData(sce_)[, whichType]
    celltypes <- ifelse(celltypes %in% selectCelltypes, celltypes, "Background")

    celltypes <- fct_relevel(celltypes, "Background") # move 'Background' to first

    myPalette <- brewer.pal(length(levels(celltypes)) - 1, "Dark2")
    myPalette <- c("grey", myPalette)
    myPalette <- alpha(myPalette, ifelse(myPalette == "grey", 0.5, 1))

    plotdf <- as.data.frame(matrix(data = NA, nrow = length(celltypes), ncol = 4))
    colnames(plotdf) <- c("x", "y", "Identity", "value")
    plotdf["x"] <- coorList[[1]]
    plotdf["y"] <- coorList[[2]]
    plotdf["Identity"] <- celltypes

    p <- ggplot(plotdf, aes(x = x, y = y, color = Identity)) +
        geom_point(size = 1) +
        scale_colour_manual(values = myPalette) +
        labs(title = paste0(ROI)) +
        theme_test()

    pdf(SavePath, height = 6, width = 8)
    # pdf("test.pdf", height = 6, width = 8)
    print(p)
    dev.off()

    return(p)
}


## calculate the cell distance
getDistance <- function(coorList1, coorList2) {
    x1 <- round(coorList1[[1]], 6)
    y1 <- round(coorList1[[2]], 6)
    x2 <- round(coorList2[[1]], 6)
    y2 <- round(coorList2[[2]], 6)

    df <- matrix(data = NA, nrow = length(x2), ncol = length(x1))
    rownames(df) <- paste0(x2, ",", y2)
    colnames(df) <- paste0(x1, ",", y1)

    for (i in 1:ncol(df)) {
        x1Temp <- rep(x1[i], times = nrow(df))
        y1Temp <- rep(y1[i], times = nrow(df))

        dist <- sqrt((x1Temp - x2)^2 + (y1Temp - y2)^2)
        dist <- ifelse(dist == 0, 999, dist)
        df[, i] <- dist
    }

    return(df)
}

## plot celltype interaction
PlotCelltypeInteraction <- function(sce, ROI, celltypes, radius = 22, savePath) {
    sce_ <- sce[, sce$ID == ROI]
    c1 <- celltypes[1]
    c2 <- celltypes[2]

    allcoorList <- getCoordinate(sce_)
    x_all <- round(allcoorList[[1]], 6)
    y_all <- round(allcoorList[[2]], 6)

    sce_1 <- sce_[, sce_$SubType == c1]
    sce_2 <- sce_[, sce_$SubType == c2]

    if (dim(sce_1)[2] < 5) {
        cat(c1, " in ", ROI, " were less than 5!", "\n")
    }
    if (dim(sce_2)[2] < 5) {
        cat(c2, " in ", ROI, " were less than 5!", "\n")
    }

    coorList1 <- getCoordinate(sce_1)
    coorList2 <- getCoordinate(sce_2)

    df <- getDistance(coorList1, coorList2)
    MASK <- ifelse(df <= radius, 1, 0)

    plotdf <- as.data.frame(matrix(data = NA, nrow = length(x_all), ncol = 3))
    colnames(plotdf) <- c("x", "y", "Subtype")
    plotdf["x"] <- x_all
    plotdf["y"] <- y_all
    plotdf["Subtype"] <- NA

    for (i in 1:nrow(MASK)) {
        for (j in 1:ncol(MASK)) {
            if (MASK[i, j] == 1) {
                c1Temp <- colnames(MASK)[j]
                c2Temp <- rownames(MASK)[i]

                c1x <- as.numeric(strsplit(c1Temp, ",")[[1]][1])
                c1y <- as.numeric(strsplit(c1Temp, ",")[[1]][2])
                c2x <- as.numeric(strsplit(c2Temp, ",")[[1]][1])
                c2y <- as.numeric(strsplit(c2Temp, ",")[[1]][2])

                index1 <- (plotdf$x == c1x) & (plotdf$y == c1y)
                plotdf$Subtype[index1] <- c1

                index2 <- (plotdf$x == c2x) & (plotdf$y == c2y)
                plotdf$Subtype[index2] <- c2
            }
        }
    }
    plotdf$Subtype[is.na(plotdf$Subtype)] <- "Background"

    p <- ggplot(plotdf, aes(x = x, y = y, color = Subtype)) +
        geom_point(size = 1) +
        labs(title = paste0("Interaction of ", c1, " and ", c2, " in ", ROI)) +
        theme_test()

    pdf(savePath, height = 6, width = 8)
    print(p)
    dev.off()

    return(p)
}

## load cell mask files and rename
LoadCellMask <- function(cellMaskPath, ROI) {
    files <- list.files(cellMaskPath)
    filesROI <- sapply(files, function(x) {
        temp <- strsplit(x, split = "_")[[1]]
        paste0(temp[2], "_", temp[3])
    })

    MaskPath <- paste0(cellMaskPath, names(filesROI[filesROI == ROI]))

    ## read image
    CellMaskMat <- ReadTiff(MaskPath)

    return(CellMaskMat)
}

## load channel files
LoadChannelImage <- function(channelIamgePath, ROI, channel) {
    pathTemp1 <- paste0(channelIamgePath, ROI, "/")
    channels <- list.files(pathTemp1)

    channelImages <- sapply(channels, function(x) {
        temp <- strsplit(x, split = "_")[[1]]
        return(temp[3])
    })

    channelsPath <- c()
    for (i in channel) {
        pathTemp2 <- paste0(pathTemp1, names(channelImages[channelImages == i]))
        channelsPath <- c(channelsPath, pathTemp2)
    }

    ## Load channel images
    ImageList <- list()
    for (i in 1:length(channelsPath)) {
        imgTemp <- png::readPNG(channelsPath[i])
        imgTemp <- TransfromtoGray(imgTemp)
        # imgTemp <- PercentileFilter(imgTemp,cutoff=0.95)
        imgTemp <- MedianFilter(imgTemp)
        imgTemp <- Transfromto0255(imgTemp)

        ImageList[[i]] <- imgTemp
    }

    if (length(ImageList) == 2) {
        ImageList[[3]] <- matrix(data = 0, nrow = nrow(ImageList[[2]]), ncol = ncol(ImageList[[2]]))
    }
    if (length(ImageList) == 1) {
        ImageList[[2]] <- matrix(data = 0, nrow = nrow(ImageList[[1]]), ncol = ncol(ImageList[[1]]))
        ImageList[[3]] <- matrix(data = 0, nrow = nrow(ImageList[[1]]), ncol = ncol(ImageList[[1]]))
    }

    ## combine different channel
    imagearray <- array(data = NA, dim = c(nrow(ImageList[[1]]), ncol(ImageList[[1]]), length(ImageList)))

    for (i in 1:length(channelsPath)) {
        imagearray[, , i] <- ImageList[[i]]
    }

    # png::writePNG(imgTemp,"median.png")

    return(imagearray)
}

## read tiff
ReadTiff <- function(filePath, filter = F) {
    tif <- readTIFF(filePath)
    tif <- Transfromto0255(tif)
    if (filter) {
        tif <- ImageFilter(tif)
    }

    return(tif)
}

## 0 - 255 tranformation
Transfromto0255 <- function(mat) {
    max_ <- max(mat)
    min_ <- min(mat)

    mat <- (mat - min_) / (max_ - min_)
    return(mat * 255)
}

## Transform RGB image into gray-scale image
TransfromtoGray <- function(array) {
    mat <- array[, , 1] + array[, , 2] + array[, , 3] / 3
    return(mat)
}

## Image filter, remove the value more than percentile
PercentileFilter <- function(mat, cutoff = 0.99) {
    cutoff <- quantile(mat, probs = cutoff)
    ifelse(mat >= cutoff, 0, mat)
    return(mat)
}

## Image fileter, median filter
MedianFilter <- function(mat, filterSize = 3) {
    rowNum <- nrow(mat)
    colNum <- ncol(mat)

    margin <- (filterSize - 1) %/% 2

    for (i in (1 + margin):(rowNum - margin)) {
        for (j in (1 + margin):(colNum - margin)) {
            seq_ <- as.numeric(mat[(i - margin):(i + margin), (j - margin):(j + margin)])
            mat[i, j] <- median(seq_)
        }
    }
    return(mat)
}

## Visualize the cellsubtype, cell mask and channel of ROI
VisTypeMaskChannel <- function(sce, ROI, celltypes, whichType = "Subtype", channel, maskPath, channelPath, SavePath) {
    if (!dir.exists(SavePath)) {
        dir.create(SavePath)
    }
    SavePath <- paste0(SavePath, ROI, "/")
    if (!dir.exists(SavePath)) {
        dir.create(SavePath)
    }

    MaskMat <- LoadCellMask(maskPath, ROI)
    # ChannelArray <- LoadChannelImage(channelPath, ROI, channel)
    PlotCelltypes(sce, ROI, celltypes, whichType = whichType, paste0(SavePath, celltypes, " in ", whichType, " on cell level.pdf"))

    png::writePNG(MaskMat, paste0(SavePath, "CellMask.png"), dpi = 100)
    # png::writePNG(ChannelArray, paste0(SavePath, channel[1], "-", channel[2], "_channel.png"), dpi = 100)
    cat(ROI, ": CellMask, celltypes and channel image were done!", "\n")
    return(NULL)
}

## Visualize the cellsubtype on mask and channel of ROI on mask (undone)
VisTypeonMask <- function(sce_, ROI, celltypes2plot, channels2plot, CellMaskPath, ChannelPath, SavePath) {
    if (!dir.exists(SavePath)) {
        dir.create(SavePath)
    }
    SavePath <- paste0(SavePath, ROI, "/")
    if (!dir.exists(SavePath)) {
        dir.create(SavePath)
    }

    MaskMat <- LoadCellMask(CellMaskPath, ROI)
    ChannelArray <- LoadChannelImage(ChannelPath, ROI, channels2plot)
    p <- PlotCelltypes(sce_, ROI, celltypes2plot, returenFigure = T)

    png::writePNG(MaskMat, paste0(SavePath, "CellMask.png"), dpi = 100)
    png::writePNG(ChannelArray, paste0(SavePath, channels2plot[1], "-", channels2plot[2], "_channel.png"), dpi = 100)

    ## mask lay on celltype
    mask <- image_read(paste0(SavePath, "CellMask.png"))

    pdf("test.pdf")
    print(p)
    grid.raster(mask)
    dev.off()

    cat(ROI, ": CellMask, celltypes and channel image were done!", "\n")
    return(NULL)
}

## Visualize the CA-9 expression value and Tumor cells on ROI
VisualizeMarkerAndType <- function(sce, ROI, celltype, whichType = "Subtype", marker, SavePath) {
    SavePath <- paste0(SavePath, ROI, "/")
    if (!dir.exists(SavePath)) {
        dir.create(SavePath)
    }

    celltypesname <- paste(celltype, collapse = "-")
    ## Plot type
    PlotCelltypes(sce = sce, ROI = ROI, selectCelltypes = celltype, whichType = whichType, SavePath = paste0(SavePath, celltypesname, " in ", whichType, " on cell level.pdf"))

    ## Plot marker
    if (!is.null(marker)) {
        PlotMarker(sce = sce, ROI = ROI, Marker = marker, SavePath = paste0(SavePath, marker, " of ", celltypesname, " on cell level.pdf"))
    }
    cat(paste0(ROI, " is plot done!"), "\n")
    return(NULL)
}
