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

## Load CNP
scimapResult <- read.csv("/mnt/data/lyx/IMC/analysis/structure/cellular_neighbor_All.csv")
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
if (T) {
    sce_ <- sce[, sce$Tissue %in% c("IM", "TAT")]
    sce_ <- sce_[, sce_$SubType %in% c("Treg")]
    mat <- as.data.frame(t(assay(sce_)[metaMarkers, ]))
    mat$phenoLabel <- sce_$Tissue

    mat$phenoLabel <- ifelse(mat$phenoLabel %in% c("IM"), "1", "0")

    FCDF <- FCandPvalueCal(mat, xCol = c(1, length(metaMarkers)), yCol = (length(metaMarkers) + 1), need.sample = TRUE)
    VolcanoPlot(FCDF, pthreshold = 0.05, fcthreshold = 1.4, feature = "Treg in different tissue", filename = paste0("/mnt/data/lyx/IMC/analysis/reclustering/Treg markers diff in IM and TAT.pdf"), Qvalue = TRUE)
}

## Treg reclustering analysis
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

## investigate the tissue specificity of Active Treg
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

    p <- ggplot(plotdf2, aes(x = Tissue, y = Abundance)) +
        geom_boxplot(alpha = 0.7, color = "black", fill = "#56B4E9", outlier.shape = NA) + # move fill color here
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

## Plot the location of Active Treg
if (T) {
    sce <- sce[, sce$Tissue == "IM"]
    sce_ <- readRDS(paste0(savePathTemp1, "Treg_recluster.rds"))

    head(colData(sce_))
    sce__ <- sce_[, sce_$Treg_Type == "Pheno_pos"]
    sce__ <- sce__[, sce__$Tissue == "IM"]

    # sort(table(sce__$ID), decreasing = T)

    source("./spatial_analysis_functions.r")

    PD1TregID <- sce__$CellID
    sce[, match(PD1TregID, sce$CellID)]$SubType <- "Active Treg"

    ROI_num <- sort(table(sce__$ID), decreasing = T)
    ROIs <- names(ROI_num[ROI_num > 20])

    celltype <- c("Active Treg")
    for (ROI in ROIs) {
        VisualizeMarkerAndType(sce = sce, ROI = ROI, celltype = celltype, whichType = "SubType", marker = NULL, SavePath = savePathTemp1)
    }
}

## investigate the spatial interaction of Active Treg
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

    sce[, match(sce_$CellID, sce$CellID)]$SubType <- "Active Treg"
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
    k_NeighborsList <- CalKNNBySpatstat(sce, ROIs = ROIs, CenType = "Active Treg", NeighType = "SubType", k = 10, xlim = c(0, 1000), ylim = c(0, 1000), return.cellID = TRUE)

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


    ## Distance between PD-1+ Treg to tumor boundary
    if (F) {
        sce_ <- readRDS(paste0(savePathTemp1, "Treg_recluster.rds"))
        sce_ <- sce_[, sce_$Tissue == "IM"]
        head(colData(sce_))

        sceDist <- readRDS("/mnt/data/lyx/IMC/analysis/GradientAna/sce_IM_with_distance2Tumor.rds")
        head(colData(sceDist))

        colData(sce_)$Dis2Tumor <- sceDist$Dis2Tumor[match(sce_$CellID, sceDist$CellID)]

        df <- as.data.frame(colData(sce_))

        ## plot
        df_type1 <- df[df$Treg_Type %in% c("Pheno_neg", "Background"), ]
        df_type1$Treg_Type <- "Normal Treg"

        df_type2 <- df[df$Treg_Type %in% "Pheno_pos", ]
        df_type2$Treg_Type <- "Active Treg"

        remain_sample <- names(table(df_type2$ID))[table(df_type2$ID) >= 20]


        df_type <- rbind(df_type1, df_type2)
        df_type <- df_type[df_type$ID %in% remain_sample, ]

        df_type <- as.data.frame(df_type)
        df_type <- df_type[, -c(1, 3)]

        ## Define distance breaks
        maxDis <- 1000
        breaks <- c(0, 50, 100, 200, 300, 500, maxDis)

        ## continuous
        if (T) {
            df_stat <- df_type[, c("ID", "Treg_Type", "Dis2Tumor")]

            ## Statistical Testing
            df_stat <- df_stat %>%
                group_by(ID, Treg_Type) %>%
                summarise(across(c(1:(ncol(df_stat) - 2)), mean, na.rm = TRUE))

            wilcox_test_result <- wilcox.test(Dis2Tumor ~ Treg_Type, data = df_stat)
            p_value <- wilcox_test_result$p.value

            ## Plot figure
            df_type$Treg_Type <- as.factor(df_type$Treg_Type)

            p <- ggplot(df_type, aes(x = Dis2Tumor, fill = Treg_Type, color = Treg_Type)) +
                geom_line(stat = "density", size = 1) +
                scale_color_brewer(palette = "Set1") +
                labs(x = "Distance to Tumor Boundary", y = "Density") +
                theme_bw() +
                theme(
                    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
                    axis.title.x = element_text(size = 10, face = "bold"),
                    axis.title.y = element_text(size = 10, face = "bold"),
                    axis.text = element_text(size = 12),
                    legend.title = element_text(size = 14, face = "bold"),
                    legend.text = element_text(size = 12),
                    legend.position = "right"
                )
            ## Add p-value to the plot
            p <- p + annotate("text", x = Inf, y = Inf, label = sprintf("wilcox test: p-value = %.4e", p_value), hjust = "right", vjust = "top", size = 4, color = "black")
        }

        pdf(paste0(savePathTemp1, "Density to Tumor/", "Density different between Active Treg and other treg to Tumor Boundary.pdf"), width = 8, height = 6)
        print(p)
        dev.off()
    }
}

## Grouping IM ROIs into different abundance group
if (T) {
    sce
    sce <- sce[, sce$Tissue %in% "IM"]

    structure <- "CNP20"

    TargeType <- "Treg"
    savePathTemp1 <- paste0(savePath, TargeType, "/")

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
    # celltypes <- celltypes[-12] ## remove "Active Treg"

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

    getInteracDiff(ResultPath, sce,
        celltypes = celltypes, savepath = paste0(savePathTemp1, "Active Treg Interaction Difference between relapse_FC_2.pdf"),
        do.fdr = FALSE, FC_threshold = 2, IDs1 = high_ROI, IDs2 = low_ROI
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
    PlotTypes <- c("CD8T", "Macro_HLADR", "Macro_CD163", "Macro_CD169", "Treg")
    p <- MultipleAbundanceSwarmPlot(HighCountMat1, lowCountMat1, groupsName = c("high", "low"), PlotTypes = PlotTypes, marker = "Acitve Treg", numGroup = 2, style = "box")
    pdf(paste0(savePathTemp1, "ViolinPlot of multiple celltypes in ", "Acitve Treg", " group.pdf"), height = 4.5, width = 8)
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

## Neighbors diff between Active Treg and Normal Treg
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
    mat_foldchangeMat <- FCandPvalueCal(mat, xCol = xCol, yCol = yCol, groups = c("Treg", "Acitve Treg"), need.sample = FALSE)

    VolcanoPlot(mat_foldchangeMat, pthreshold = 0.05, fcthreshold = 2, feature = NULL, filename = paste0(savePathTemp1, "Volcano of neigobors in Active Treg and Other Treg.pdf"))
}

## Mono_CD11c reclustering analysis
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
