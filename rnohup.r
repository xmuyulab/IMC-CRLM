# Identify tumor microenviroment pattern
library(SingleCellExperiment)
library(survival)
library(survminer)

library(pheatmap)
library(ggrepel)
library(ggalluvial)

library(dplyr)
library(tidyverse)

source("./spatial_analysis_functions.r")
source("./structural_analysis_functions.r")

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
names(table(sce$Tissue))

savePath <- "/mnt/data/lyx/IMC/analysis/structure/"

tissues <- c("IM", "CT", "TAT")

NeigoborDomain <- c("CNP5", "CNP10", "CNP15", "CNP20", "CNP25", "CNP30")

k <- 10

## All cells
if (T) {
    ## Bind CNPs
    scimapResult <- read.csv(paste0("/mnt/data/lyx/IMC/analysis/structure/cellular_neighbor_", "All", ".csv"))
    scimapResult <- scimapResult[, -1]

    sce <- BindResult(sce, scimapResult, NeigoborDomain)

    ## Plot the CNP contribution

    for (structure in NeigoborDomain) {
        savePath2 <- paste0(savePath, structure, "/")
        if (!dir.exists(savePath2)) {
            dir.create(savePath2, recursive = T)
        }

        ## Cell subtype fraction in cellular neighbors pattern
        HeatmapForCelltypeInNeighbor(sce, "SubType", structure, savePath2)

        ## Compare the CNP between tissue
        BarPlotForCelltypeCounts(sce = sce, tissueCol = "Tissue", groupCol = "RFS_status", typeCol = structure, savePath = savePath2)
    }
}


## Seperate Tissue
for (tissue in tissues) {
    ## Select tissue
    sce_ <- sce[, sce$Tissue == tissue]
    celltypes <- names(table(sce_$SubType))
    celltypes <- celltypes[celltypes != "UNKNOWN"]

    savePath2 <- paste0(savePath, tissue, "/")

    ## Cellular neighbors analysis in different domain size
    for (structure in NeigoborDomain) {
        ## savepath
        savePath3 <- paste0(savePath2, structure, "/")
        if (!dir.exists(savePath3)) {
            dir.create(savePath3, recursive = T)
        }

        ## Cell subtype fraction in cellular neighbors pattern
        HeatmapForCelltypeInNeighbor(sce_, "SubType", structure, savePath3)

        ## Cellular pattern difference in whole ROI between Relaps and Non-Relaps
        CompareCellularPattern(sce_, sep = "RFS_status", countcol = structure, n_cluster = 10, clinicalFeatures = NULL, savePath = savePath3)

        ## Celllular neighborhood pattern survival analysis
        CNP_countsDF <- GetAbundance(sce_, countcol = structure, clinicalFeatures = c("RFS_status", "RFS_time"), is.fraction = TRUE, is.reuturnMeans = T)
        CNPs <- names(table(colData(sce_)[, structure]))
        if (!dir.exists(paste0(savePath3, "KM/"))) {
            dir.create(paste0(savePath3, "KM/"), recursive = T)
        }
        for (CNP in CNPs) {
            plotdf <- CNP_countsDF[, c(CNP, "RFS_status", "RFS_time")]
            KMForCNP(plotdf, CNP, savePath = paste0(savePath3, "KM/", "Cellular Neighborhood pattern Suvival analysis of ", CNP, ".pdf"))
        }

        ## Sankey diagram to display the cells flow
        if (T) {
            interstingType <- c("B", "CD4T", "CD8T", "NK", "Treg", "Macro_CD11b", "Macro_CD163", "Macro_CD169", "Macro_HLADR", "Mono_CD11c")

            df <- as.data.frame(colData(sce_))
            df <- df[, c("SubType", structure, "RFS_status")]
            colnames(df) <- c("SubType", "CNP", "RFS_status")

            df$RFS_status <- factor(df$RFS_status, levels = c(1, 0), labels = c("Recurrence", "No Recurrence"))
            df <- df[df[, 1] %in% interstingType, ]

            ## sampling
            set.seed(619)
            sampleidx <- sample.int(nrow(df), size = 5000)
            dfplot <- df[sampleidx, ]

            my_color <- c()
            my_color <- c(my_color, pal_igv("default")(length(table(dfplot[, 1]))))
            my_color <- c(my_color, pal_jco("default")(length(table(dfplot[, 2]))))
            my_color <- c(my_color, pal_nejm("default")(length(table(dfplot[, 3]))))

            p <- ggplot(data = dfplot, aes(axis1 = SubType, axis2 = CNP, axis3 = RFS_status)) +
                geom_alluvium(aes(fill = RFS_status)) +
                geom_stratum(aes(fill = after_stat(stratum))) +
                geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2, color = "black") +
                scale_x_discrete(limits = c("SubType", "CNP", "RFS_status"), expand = c(0.15, 0.05)) +
                scale_fill_manual(values = my_color) +
                theme_void() +
                theme(legend.position = "none")

            pdf(paste0(savePath3, "Sankey of celltype_CNP_RFS in ", structure, ".pdf"), height = 8, width = 6)
            print(p)
            dev.off()
        }

        ## Cellular pattern difference in certain central type between relapse and non relapse
        meta <- colData(sce_)

        meta_ <- meta[, match(c("ID", structure, "SubType", "RFS_status"), colnames(meta))]
        meta_ <- as.data.frame(meta_)
        colnames(meta_) <- c("ID", "CNP", "CentralType", "RFS_status")

        df_fraction <- meta_ %>%
            group_by(ID, CentralType, RFS_status, CNP) %>%
            summarise(count = n()) %>%
            mutate(fraction = count / sum(count)) %>%
            select(-count)

        # Reshape the dataframe
        df_wide <- df_fraction %>%
            pivot_wider(names_from = CNP, values_from = fraction, values_fill = 0)

        visuaDF <- as.data.frame(matrix(data = 0, nrow = 0, ncol = 5))

        for (celltype in celltypes) {
            df_wideTemp <- df_wide[df_wide$CentralType %in% celltype, ]
            mat <- df_wideTemp[, c(4:13, 3)]
            mat <- as.data.frame(mat)
            xCol <- c(1, k)
            yCol <- k + 1
            mat_foldchangeMat <- FCandPvalueCal(mat, xCol = xCol, yCol = yCol)
            mat_foldchangeMat$Q.value <- p.adjust(mat_foldchangeMat$P.value, method = "BH")
            mat_foldchangeMat <- cbind(rep(celltype, nrow(mat_foldchangeMat)), mat_foldchangeMat)

            visuaDF <- rbind(visuaDF, mat_foldchangeMat)
        }

        ### Visualize
        if (T) {
            colnames(visuaDF) <- c("CentralType", "CNP", "FC", "P.value", "Q.value")
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

            visuaDF$label <- ifelse(visuaDF$dir != "", visuaDF$CNP, "")

            mycol <- c(ggsci::pal_npg("nrc", alpha = 0.8)(2), "gray")
            names(mycol) <- c("up-regulated", "down-regulated", "NOT")

            visuaDF$log2FC <- log2(visuaDF$FC)

            p1 <- ggplot(visuaDF, aes(x = CentralType, y = log2FC)) +
                geom_jitter(aes(x = CentralType, y = log2FC, color = dir), size = 0.2, width = 0.3) +
                theme_classic() +
                geom_text_repel(aes(label = label), size = 3) +
                scale_color_manual(values = mycol) +
                ylab("log2FC") +
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

            pdf(paste0(savePath3, "CNP Difference of Central types between Relapse.pdf"), width = 8, height = 6)
            print(p1)
            dev.off()
        }

        ## Plot CNP on image
        if (T) {
            SavePath4 <- paste0(savePath3, "CNP_oncells/")
            if (!dir.exists(SavePath4)) {
                dir.create(SavePath4, recursive = T)
            }
            colData(sce)[, structure] <- as.factor(colData(sce)[, structure])

            ROIs <- names(table(colData(sce)$ID))
            for (ROI in ROIs) {
                PlotCNP(sce, ROI, TypeCol = structure, SavePath = paste0(SavePath4, ROI, "_"))
            }
        }
    }
}
