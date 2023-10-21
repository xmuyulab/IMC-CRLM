# functions for cell type abundance analysis
library(SingleCellExperiment)

library(ggplot2)
library(ggpubr)
library(ggcor)
library(ggridges)
library(ggpointdensity)
library(ggrepel)
library(patchwork)
library(ggthemes)
library(ggsci)
library(ggunchained)
library(RColorBrewer)

library(survival)
library(survminer)
library(tableone)
library(forestplot)

library(tidyverse)
library(dplyr)

## Plot the marker distribution from different batch
PlotMarkerExpDistribution <- function(sce, savePath) {
  ## create expression matrix
  exp <- assay(sce)

  genes <- rep(rownames(exp), times = ncol(exp))
  values <- as.numeric(exp)
  batches <- rep(sce$Batch, each = nrow(exp))

  plotdf <- cbind("Marker" = genes, "Value" = values, "Batch" = batches)
  plotdf <- as.data.frame(plotdf)
  plotdf$Value <- as.numeric(plotdf$Value)

  p <- ggplot(plotdf, aes(x = Value, y = Batch, fill = Batch)) +
    geom_density_ridges(alpha = 0.25) +
    theme_ridges() +
    facet_wrap(~Marker, ncol = 6)

  pdf(paste0(savePath, "Marker Expression distribution from batches.pdf"), height = 9, width = 12)
  print(p)
  dev.off()

  return(NULL)
}

## transform annotated data into SingleCellExperiment object
SCE_Transform <- function(scedf, assay_col, cellmeta_col, clinical = NULL) {
  cat("The dim of dataframe is ", dim(scedf), "\n")
  assay_ <- t(as.matrix(scedf[, assay_col[1]:assay_col[2]]))
  cellmeta_ <- scedf[, cellmeta_col[1]:cellmeta_col[2]]
  cellmeta_ <- dplyr::left_join(cellmeta_, clinical, by = "PID")

  sce <- SingleCellExperiment(
    assays = assay_,
    colData = cellmeta_
  )

  return(sce)
}

## Multiple clinical cross boxplot
CrossBoxplotForAbundance <- function(countdf, celltype, TypeName = NULL, clinicalFeatures, savePath) {
  if (length(celltype) > 1) {
    countdf[, TypeName] <- apply(countdf[, celltype], MARGIN = 1, FUN = "sum")
  } else {
    TypeName <- celltype
  }
  df <- countdf[, c(TypeName, clinicalFeatures)]
  colnames(df)[1] <- "Abundance"

  # Create long format dataframe for clinical features
  long_df <- df %>% tidyr::pivot_longer(cols = starts_with(clinicalFeatures[-1]), names_to = "Feature", values_to = "Value")
  long_df$Group <- paste0(long_df$Tissue, "_", long_df$Feature, "_", long_df$Value)

  # Compute percentiles and median for each Group (Tissue), Celltype, and Feature
  summary_stats <- long_df %>%
    group_by(Group) %>%
    summarise(
      Q1 = quantile(Abundance, 0.25),
      Median = median(Abundance),
      Q3 = quantile(Abundance, 0.75)
    )

  # Add the approximate x-axis positions for the labels and dashed lines
  total_features <- length(unique(long_df$Group))
  features_per_class <- total_features / 3
  label_positions <- seq(features_per_class / 2, (total_features - features_per_class / 2), length.out = 3)
  dashed_line_positions <- seq(features_per_class, total_features - features_per_class, length.out = 2)
  class_labels <- data.frame(
    Class = c("CT", "IM", "TAT"),
    LabelPosition = label_positions
  )
  # Create a named vector for new labels
  new_labels <- c(
    "CT_Age_60_0" = "Age<60", "IM_Age_60_0" = "Age<60", "TAT_Age_60_0" = "Age<60",
    "CT_Age_60_1" = "Age>=60", "IM_Age_60_1" = "Age>=60", "TAT_Age_60_1" = "Age>=60",
    "CT_RFS_status_0" = "Non-Relapse", "IM_RFS_status_0" = "Non-Relapse", "TAT_RFS_status_0" = "Non-Relapse",
    "CT_RFS_status_1" = "Relapse", "IM_RFS_status_1" = "Relapse", "TAT_RFS_status_1" = "Relapse",
    "CT_Gender_0" = "Male", "IM_Gender_0" = "Male", "TAT_Gender_0" = "Male",
    "CT_Gender_1" = "Female", "IM_Gender_1" = "Female", "TAT_Gender_1" = "Female",
    # "TAT_fong_score_4_1" = "FongScore>=4", "IM_fong_score_4_1" = "FongScore>=4", "CT_fong_score_4_1" = "FongScore>=4",
    # "CT_fong_score_4_0" = "FongScore<4", "IM_fong_score_4_0" = "FongScore<4", "TAT_fong_score_4_0" = "FongScore<4",
    "TAT_KRAS_mutation_1" = "KRAS Mut", "IM_KRAS_mutation_1" = "KRAS Mut", "CT_KRAS_mutation_1" = "KRAS Mut",
    "CT_KRAS_mutation_0" = "WT", "IM_KRAS_mutation_0" = "WT", "TAT_KRAS_mutation_0" = "WT",
    "TAT_TBS_8_0" = "TBS<8", "CT_TBS_8_0" = "TBS<8", "IM_TBS_8_0" = "TBS<8",
    "IM_TBS_8_1" = "TBS>=8", "CT_TBS_8_1" = "TBS>=8", "TAT_TBS_8_1" = "TBS>=8",
    "TAT_Differential_grade_1" = "Moderately Differntiate", "CT_Differential_grade_1" = "Moderately Differntiate", "IM_Differential_grade_1" = "Moderately Differntiate",
    "IM_Differential_grade_0" = "Lowly Differntiate", "CT_Differential_grade_0" = "Lowly Differntiate", "TAT_Differential_grade_0" = "Lowly Differntiate"
  )

  # Create custom cross-boxplot for each clinical feature group
  p <- ggplot(long_df, aes(x = Group, y = Abundance, group = Feature)) +
    geom_point(aes(color = Group), alpha = 0.5, size = 0.25, position = position_jitter(width = 0.15)) +
    geom_crossbar(data = summary_stats, aes(x = Group, y = Median, ymin = Median, ymax = Median), width = 0.5, color = "black", size = 1.0, inherit.aes = FALSE) +
    geom_errorbar(data = summary_stats, aes(x = Group, ymin = Q1, ymax = Q3), width = 0, color = "black", size = 0.25, inherit.aes = FALSE) +
    geom_text(data = class_labels, aes(x = LabelPosition, y = (max(long_df$Abundance) + 0.25), label = Class), size = 5, inherit.aes = FALSE) +
    geom_vline(xintercept = (dashed_line_positions + 0.5), linetype = "dashed", color = "black", size = 1) +
    scale_x_discrete(labels = new_labels) +
    scale_color_igv() +
    labs(x = "Tissue", y = "Abundance", title = "Cell Subpopulation Differences for Each Clinical Feature Group") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 14, face = "bold"),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(15, 15, 15, 15)
    )

  pdf(paste0(savePath, "Abudance Cross boxplot of ", TypeName, ".pdf"), height = 6, width = 8)
  print(p)
  dev.off()
}

## transfer sce to cell counts matrix
Transform_CellCountMat <- function(sceobj, group = c("IM", "CT"), clinicalFeatures, type = "SubType", is.fraction = FALSE) {
  if (length(group) == 2) {
    sceobj <- sceobj[, sceobj$Tissue == group[1] | sceobj$Tissue == group[2]]
  }
  if (length(group) == 3) {
    sceobj <- sceobj[, sceobj$Tissue == group[1] | sceobj$Tissue == group[2] | sceobj$Tissue == group[3]]
  }

  cellMeta <- colData(sceobj)

  ## ROI, major celltype and cell subtype names and other clinical information
  ROIs <- names(table(cellMeta$ID))

  SubTypes <- names(table(cellMeta[, type]))
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
    SubTem <- as.data.frame(t(table(coldataTemp[, type])))
    if (is.fraction) {
      CellCountMat[match(ROI, rownames(CellCountMat)), match(SubTem$Var2, colnames(CellCountMat))] <- SubTem$Freq / cellnum
    }
    if (!is.fraction) {
      CellCountMat[match(ROI, rownames(CellCountMat)), match(SubTem$Var2, colnames(CellCountMat))] <- SubTem$Freq
    }
  }

  ## fill the NA
  for (i in 1:ncol(CellCountMat)) {
    CellCountMat[, i][is.na(CellCountMat[, i])] <- 0
  }

  ## match the row of clinical and plotdf
  CellCountMat$PID <- as.vector(sapply(rownames(CellCountMat), function(x) {
    strsplit(x, "_")[[1]][1]
  }))

  for (feature in clinicalFeatures) {
    CellCountMat[, feature] <- cellMeta[match(rownames(CellCountMat), cellMeta$ID), feature]
  }

  return(CellCountMat)
}

## Volcano Plot
abundanceVolcanoPlot <- function(countdf, pthreshold, fcthreshold, xCol, clinicalFeatures, savePath) {
  ## grouping data
  cat("The tissues is: ", names(table(countdf$Tissue)), "\n")

  ## savePath
  savePath <- paste0(savePath, "VolcanoPlots/")
  if (!dir.exists(savePath)) {
    dir.create(savePath)
  }

  for (feature in clinicalFeatures) {
    cat("Process feature: ", feature, "\n")
    mat1 <- subset(countdf, Tissue == "IM")
    colnames(mat1)
    mat1_foldchangeMat <- FCandPvalueCal(mat1, xCol = xCol, yCol = feature)
    if (length(mat1_foldchangeMat) != 1) {
      VolcanoPlot(mat1_foldchangeMat, pthreshold = pthreshold, fcthreshold = fcthreshold, feature, filename = paste0(savePath, "IM_CellFraction_Volcano of ", feature, ".pdf"))
    }

    mat2 <- subset(countdf, Tissue == "CT")
    colnames(mat2)
    mat2_foldchangeMat <- FCandPvalueCal(mat2, xCol = xCol, yCol = feature)
    if (length(mat2_foldchangeMat) != 1) {
      VolcanoPlot(mat2_foldchangeMat, pthreshold = pthreshold, fcthreshold = fcthreshold, feature, filename = paste0(savePath, "CT_CellFraction_Volcano of ", feature, ".pdf"))
    }

    mat3 <- subset(countdf, Tissue == "TAT")
    colnames(mat3)
    mat3_foldchangeMat <- FCandPvalueCal(mat3, xCol = xCol, yCol = feature)
    if (length(mat3_foldchangeMat) != 1) {
      VolcanoPlot(mat3_foldchangeMat, pthreshold = pthreshold, fcthreshold = fcthreshold, feature, filename = paste0(savePath, "TAT_CellFraction_Volcano of ", feature, ".pdf"))
    }
  }

  return(NULL)
}

abundanceDotPlotMat <- function(countdf, xCol, clinicalFeatures, tissue) {
  ## grouping data
  cat("Abundance analysis for multiple clinical features.", "\n")

  returnDF <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 5))

  for (feature in clinicalFeatures) {
    cat("Process feature: ", feature, "\n")
    mat1 <- subset(countdf, Tissue == tissue)
    mat1_foldchangeMat <- FCandPvalueCal(mat1, xCol = xCol, yCol = feature)
    if (class(mat1_foldchangeMat) == "numeric") {
      cat("The feature ", feature, " do not contain two groups!", "\n")
      next
    }
    mat1_foldchangeMat$Feature <- rep(feature, times = nrow(mat1_foldchangeMat))
    mat1_foldchangeMat$Tissue <- rep(tissue, times = nrow(mat1_foldchangeMat))
    returnDF <- rbind(returnDF, mat1_foldchangeMat)
  }
  colnames(returnDF) <- c("Celltype", "Foldchange", "P.value", "Feature", "Tissue")
  returnDF <- summaryClinical(returnDF)

  return(returnDF)
}

summaryClinical <- function(DF) {
  DF$FeatureGroup <- NA
  features <- names(table(DF$Feature))
  for (feature in features) {
    if (feature == "Prognosis") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "Worse versus Good"
    }
    if (feature == "RFS_status") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "Relapse versus Non-relapse"
    }
    if (feature == "RFS_time_12") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "RFS_time>=12 versus <12"
    }
    if (feature == "RFS_time_24") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "RFS_time>=24 versus <24"
    }
    if (feature == "Recurrence_site") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "Extrahepatic relapse versus Liver relapse"
    }
    if (feature == "zs_rec_riskmodel") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "High-risk versus Low-risk"
    }
    if (feature == "fong_score_3") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "FongScore>=3 versus <3"
    }
    if (feature == "Gender") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "Female versus Male"
    }
    if (feature == "Age_60") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "Age>=60 versus <60"
    }
    if (feature == "KRAS_mutation") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "KRAS Mut versus WT"
    }
    if (feature == "NAS_mutation") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "NAS Mut versus WT"
    }
    if (feature == "BRAF_mutation") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "BRAF Mut versus WT"
    }
    if (feature == "TBS_3") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "TBS>=3 versus <3"
    }
    if (feature == "TBS_8") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "TBS>=8 versus <8"
    }
    if (feature == "CRLM_number_4") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "CRLM_number>=4 versus <4"
    }
    if (feature == "CRLM_size_3") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "CRLM_size>=3 versus <3"
    }
    if (feature == "CRLM_size_5") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "CRLM_size>=5 versus <5"
    }
    if (feature == "Liver_involvement_num") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "Liver involvement 2 versus involvement 1"
    }
    if (feature == "CEA_5") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "CEA>=5 versus CEA<5"
    }
    if (feature == "CA199_37") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "CA199>=37 versus <37"
    }
    if (feature == "Pathology") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "Mucinous versus Adenocarcinoma"
    }
    if (feature == "Differential_grade") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "Well-differentiated versus Moderately-differentiated"
    }
    if (feature == "Lymph_positive") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "Lymphonode-positive versus negative"
    }
  }
  return(DF)
}

## calculate fold-change and p-value
FCandPvalueCal <- function(mat, xCol, yCol) {
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
    labs(x = "log2(Fold Change)", y = "-log10(P value)", title = paste0("Volcano Plot of Different Expression Celltypes in ", feature)) +
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

## Multiple clinical features Dotplot
MultiCliDotplot <- function(plotdf, tissue, savePath) {
  plotdf$P.label <- ifelse(plotdf$P.value <= 0.001, "***", ifelse(plotdf$P.value <= 0.01, "**", ifelse(plotdf$P.value <= 0.05, "*", "n.s.")))
  plotdf$lgPvalue <- ifelse(plotdf$P.value <= 0.001, 4, ifelse(plotdf$P.value <= 0.01, 3, ifelse(plotdf$P.value <= 0.05, 2, 1)))
  plotdf$log2Foldchange <- log2(as.numeric(plotdf$Foldchange))

  plotdf$log2Foldchange <- ifelse(plotdf$log2Foldchange > 1.5, 1.5, plotdf$log2Foldchange)
  plotdf$log2Foldchange <- ifelse(plotdf$log2Foldchange < (-1.5), (-1.5), plotdf$log2Foldchange)

  p <- ggplot(plotdf, aes(x = Celltype, y = FeatureGroup, size = lgPvalue, colour = log2Foldchange)) +
    geom_point() +
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_line(colour = "white"),
      panel.border = element_rect(colour = "white", fill = NA),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    scale_color_gradient2(low = "royalblue", mid = "white", high = "firebrick3")

  pdf(paste0(savePath, "Multiple clinical features abundance differnt of ", tissue, ".pdf"), height = 6, width = 10)
  print(p)
  dev.off()
  return(NULL)
}

## Boxplot for cell subpopulations under certain clinical groups
AbundanceBoxPlot <- function(sce, countdf, celltypes2Plot, expCol, tissueCol, clinicalGroupCol, ClinicalGroupName) {
  ## To long form
  AbunBoxDF <- pivot_longer(countdf, cols = expCol[1]:expCol[2], values_to = "Abundance", names_to = "Celltype")

  AbunBoxDF[, clinicalGroupCol] <- as.factor(ifelse(AbunBoxDF[, clinicalGroupCol] == 1, ClinicalGroupName[1], ClinicalGroupName[2]))
  colnames(AbunBoxDF)[match(clinicalGroupCol, colnames(AbunBoxDF))] <- "ClinicalGroup"

  AbunBoxDF$MajorType <- colData(sce)[match(AbunBoxDF$Celltype, sce$SubType), "MajorType"]

  p <- ggplot(AbunBoxDF, aes(x = Celltype, y = Abundance, fill = ClinicalGroup)) +
    geom_split_violin(trim = T, draw_quantiles = 0.5) +
    scale_y_continuous(name = "Cell Abundance") +
    scale_x_discrete(name = "Cell Population") +
    scale_fill_lancet() +
    stat_compare_means(aes(group = ClinicalGroup),
      method = "t.test",
      hide.ns = TRUE, # Hide non-significant comparisons
      label = "p.signif",
      label.y.npc = "middle"
    ) +
    facet_grid(Tissue ~ MajorType, scales = "free") +
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

  return(p)
}

## Boxplot for some cell subpopulations abudance comparison
AbundanceBoxPlot2 <- function(countdf, celltypes2Plot, expCol, tissueCol, tissue, clinicalGroupCol, ClinicalGroupName) {
  ## To long form
  AbunBoxDF <- pivot_longer(countdf, cols = expCol[1]:expCol[2], values_to = "Abundance", names_to = "Celltype")

  AbunBoxDF[, clinicalGroupCol] <- as.factor(ifelse(AbunBoxDF[, clinicalGroupCol] == 1, ClinicalGroupName[1], ClinicalGroupName[2]))
  colnames(AbunBoxDF)[match(clinicalGroupCol, colnames(AbunBoxDF))] <- "ClinicalGroup"

  AbunBoxDF <- as.data.frame(AbunBoxDF)

  AbunBoxDF <- AbunBoxDF[AbunBoxDF[, tissueCol] %in% tissue, ]
  AbunBoxDF <- AbunBoxDF[AbunBoxDF[, "Celltype"] %in% celltypes2Plot, ]


  p <- ggplot(AbunBoxDF, aes(x = Tissue, y = Abundance, fill = ClinicalGroup)) +
    geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
    scale_fill_manual(values = ggsci::pal_jco("default")(2))  +
    facet_wrap(~Celltype, nrow = 2, scales = "free_y") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      text = element_text(size = 12),
      axis.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5),
      strip.background = element_blank()
    )+
    stat_compare_means(aes(group = ClinicalGroup),
      method = "t.test",
      hide.ns = FALSE,
      label = "p.signif",
      label.y.npc = "middle"
    )

  return(p)
}

## Convert count matrix into plot dataframe
CountMat2Plotdf <- function(plotdf, MetaCol, expCol) {
  ### matrix transform
  abundanceVec <- as.numeric(as.matrix(plotdf[, expCol[1]:expCol[2]]))
  cellsubtypeVec <- rep(colnames(plotdf[, expCol[1]:expCol[2]]), each = nrow(plotdf))
  metaList <- list()
  for (i in 1:length(MetaCol)) {
    metaList[[MetaCol[i]]] <- rep(as.character(plotdf[, MetaCol[i]]), times = ncol(plotdf))
  }

  plotdf2 <- cbind(abundanceVec, cellsubtypeVec, metaList[[1]])
  plotdf2 <- as.data.frame(plotdf2)
  plotdf2[, 1] <- as.numeric(plotdf2[, 1])

  for (i in 2:length(metaList)) {
    plotdf2 <- cbind(plotdf2, metaList[[i]])
    plotdf2[, i] <- as.factor(plotdf2[, i])
  }

  colnames(plotdf2) <- c("Abundance", "Celltype", MetaCol)

  return(plotdf2)
}

## Correlation analysis of RFS_time and cell abundance
abundanceSurvivalCorrelation <- function(plotdf) {
  plotdf <- CountMat2Plotdf(plotdf)

  IMplotdf <- subset(plotdf, Tissue == "IM")

  p <- ggplot(IMplotdf, aes(x = RecurrenceTime, y = Abundance, color = Recurrence)) +
    geom_point(size = 0.5, shape = 16) +
    geom_smooth(method = "lm", se = FALSE, color = "turquoise4") +
    facet_wrap(~Celltype)
  pdf("test.pdf", height = 8, width = 12)
  print(p)
  dev.off()
}

## meta analysis
MultipleUniCOX <- function(df, do.scale = FALSE, UniCOX = TRUE) {
  result <- matrix(data = NA, nrow = 0, ncol = 6)
  result <- as.data.frame(result)

  colnames(df) <- sapply(colnames(df), function(x) {
    sub(pattern = "\\+", replacement = "", x)
  })
  colnames(df) <- sapply(colnames(df), function(x) {
    sub(pattern = " ", replacement = "_", x)
  })

  features <- colnames(df)[1:(ncol(df) - 2)]
  cat("The features in multi-cox are: ", features, "\n")

  # features <- features[-(match(c("RFS_time", "RFS_status"), features))]

  if (do.scale) {
    continuousVar <- c()
    for (feature in features) {
      a <- df[, match(feature, colnames(df))]
      if (length(table(a)) >= nrow(df) / 2) {
        continuousVar <- c(continuousVar, feature)
      }
    }

    df[, match(continuousVar, colnames(df))] <- apply(df[, match(continuousVar, colnames(df))], MARGIN = 2, function(a) {
      scale(a)
    })
  }

  if (UniCOX) {
    univ_formulas <- sapply(
      features,
      function(x) {
        as.formula(paste0("Surv(RFS_time, RFS_status)~", x))
      }
    )

    univ_models <- lapply(univ_formulas, function(x) {
      coxph(x, data = df)
    })

    univ_results <- lapply(univ_models, function(x) {
      ftest <- cox.zph(x)
      ftest
      mul_cox1_result <- summary(x)
      multi1 <- round(mul_cox1_result$conf.int[, c(1, 3, 4)], 4)
      multi1 <- as.data.frame(t(multi1))

      multi2 <- ShowRegTable(
        x,
        exp = TRUE,
        digits = 2, pDigits = 3,
        printToggle = TRUE, quote = FALSE, ciFun = confint
      )

      result <- cbind(multi1, multi2)
      result <- cbind(Features = rownames(result), result)

      return(result)
    })

    for (x in univ_results) {
      result <- rbind(result, x)
    }

    return(result)
  } else {
    formula_ <- paste0("Surv(RFS_time, RFS_status)~", features[1])
    for (i in 2:length(features)) {
      formula_ <- paste0(formula_, "+", features[i])
    }
    mul_cox <- coxph(as.formula(formula_), data = df)
    mul_cox1 <- summary(mul_cox)

    multi1 <- as.data.frame(round(mul_cox1$conf.int[, c(1, 3, 4)], 2))
    multi2 <- ShowRegTable(mul_cox,
      exp = TRUE,
      digits = 2,
      pDigits = 3,
      printToggle = TRUE,
      quote = FALSE,
      ciFun = confint
    )
    result <- cbind(Features = rownames(multi2), multi1, multi2)
    return(result)
  }
}

abundanceMetaAnalysis <- function(plotdf, celltypes2Plot, clinical, features, tissue, savePath) {
  ## PIDs
  plotdf <- as.data.frame(plotdf)
  plotdf$PID <- rownames(plotdf)

  features_ <- c("PID", features)

  clinical2 <- clinical[match(plotdf$PID, clinical$PID), features_]

  clinical_features <- features
  print("Clinical Features: ")
  print(clinical_features)

  rownames(clinical2) <- c(1:nrow(clinical2))

  ## bind more clinical information
  plotdf <- dplyr::left_join(plotdf, clinical2, by = "PID")
  plotdf <- plotdf[, c(celltypes2Plot, features)]

  result <- MultipleUniCOX(plotdf)

  ## Definate space
  ins <- function(x) {
    c(as.character(x), rep(NA, ncol(result) - 1))
  }

  numcelltype <- length(celltypes2Plot)

  ## result matrix
  if (1) {
    result_df <- rbind(
      c("Features", NA, NA, NA, "HR(95%CI)", "p-value"),
      ins("Cell Population Abundance"),
      result[c(1:numcelltype), ],
      ins("Other Clnical Information"),
      result[c((numcelltype + 1):nrow(result)), ],
      c(NA, NA, NA, NA, NA, NA)
    )
  }

  ## hight-light rows
  is_summary_vector <- c()
  for (i in 1:nrow(result_df)) {
    if (is.na(result_df[i, 2])) {
      is_summary_vector <- c(is_summary_vector, TRUE)
    } else {
      is_summary_vector <- c(is_summary_vector, FALSE)
    }
  }

  ## forest plot
  p <- forestplot(result_df[, c(1, 5, 6)],
    mean = as.numeric(result_df[, 2]),
    lower = as.numeric(result_df[, 3]),
    upper = as.numeric(result_df[, 4]),
    zero = 1,
    boxsize = 0.6,
    graph.pos = "right",
    hrzl_lines = list(
      "1" = gpar(lty = 1, lwd = 2),
      "2" = gpar(lty = 2) # ,
      # "21" = gpar(lwd = 2, lty = 1, columns = c(1:4))
    ),
    graphwidth = unit(.25, "npc"),
    xlab = "HR(exp(coef))",
    xticks = c(0.4, 1, 3, 5, 7, 10),
    is.summary = is_summary_vector,
    txt_gp = fpTxtGp(
      label = gpar(cex = 1),
      ticks = gpar(cex = 1),
      xlab = gpar(cex = 1.5),
      title = gpar(cex = 2)
    ),
    lwd.zero = 1,
    lwd.ci = 1.5,
    lwd.xaxis = 2,
    lty.ci = 1.5,
    ci.vertices = T,
    ci.vertices.height = 0.2,
    clip = c(0.1, 8),
    ineheight = unit(8, "mm"),
    line.margin = unit(8, "mm"),
    colgap = unit(6, "mm"),
    fn.ci_norm = "fpDrawDiamondCI",
    title = "Multi-variables forest plot",
    col = fpColors(
      box = "#021eaa",
      lines = "#021eaa",
      zero = "black"
    )
  )

  pdf(paste0(savePath, "Multi-uniCOX of ", tissue, ".pdf"), width = 12, height = 9)
  print(p)
  dev.off()

  return(NULL)
}

## Kaplan-Meier curve
KMVisualize <- function(df, celltype, cutoff = "best", savePath = NULL) {
  df <- df[, c(celltype, "RFS_status", "RFS_time")]

  if (cutoff == "median") {
    df[, 1] <- ifelse(df[, 1] >= median(df[, 1]), "high", "low")
  }
  if (cutoff == "mean") {
    df[, 1] <- ifelse(df[, 1] >= mean(df[, 1]), "high", "low")
  }
  if (cutoff == "best") {
    cutpoint <- surv_cutpoint(data = df, time = "RFS_time", event = "RFS_status", variables = celltype)
    cutpoint <- summary(cutpoint)$cutpoint
    df[, 1] <- ifelse(df[, 1] >= cutpoint, "high", "low")
  }

  ## km curve
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

## correlation analysis
CorrelationAnalysis <- function(df, savePath) {
  p <- quickcor(df, cluster = TRUE, type = "upper", cor.test = TRUE) +
    geom_colour(data = get_data(type = "upper")) +
    geom_mark(data = get_data(type = "upper"), size = 3, color = "black", fontface = 1) +
    scale_fill_gradientn(colours = c("#77C034", "white", "#C388FE")) +
    geom_panel_grid(colour = "black", size = 0.5)

  pdf(savePath, height = 15, width = 15)
  print(p)
  dev.off()
  return(NULL)
}

## merge the conunts of ROIs from same patients
MergeCountsofROI <- function(countdf, tissue, expCol, scale = F) {
  countdf <- subset(countdf, Tissue == tissue)
  pids <- names(table(countdf$PID))

  PIDDF <- matrix(data = NA, nrow = length(pids), ncol = (expCol[2] - expCol[1] + 1))
  rownames(PIDDF) <- pids
  colnames(PIDDF) <- colnames(countdf)[expCol[1]:expCol[2]]

  for (pid in pids) {
    expTemp <- subset(countdf, PID == pid)
    expTemp <- expTemp[, expCol[1]:expCol[2]]
    expTemp <- apply(expTemp, MARGIN = 2, FUN = "mean")
    PIDDF[pid, ] <- expTemp
  }
  if (scale) {
    PIDDF <- scale(PIDDF, center = F, scale = T)
  }
  PIDDF <- as.data.frame(PIDDF)

  return(PIDDF)
}


## Stack barplot to show the celltype fraction
BarPlotForCelltypeFraction <- function(sce, rowSep, colSep, savePath) {
  meta <- colData(sce)

  Majortypes <- meta[, "MajorType"]
  majortypes <- names(table(Majortypes))
  Subtypes <- meta[, "SubType"]
  IDs <- meta[, "ID"]
  Tissues <- meta[, rowSep]
  Groups <- meta[, colSep]

  plotdf <- cbind(Majortypes, Subtypes, IDs, Tissues, Groups)
  plotdf <- as.data.frame(plotdf)
  colnames(plotdf) <- c("Majortype", "Subtype", "ROI", "Tissue", "Relapse")

  df <- plotdf[, c(2, 3, 4)]
  # Calculate the total cell count for each ROI and Tissue combination
  cell_count <- df %>%
    group_by(ROI, Tissue) %>%
    summarise(Total = n())

  # Calculate the cell count for each Subtype, ROI, and Tissue combination
  subtype_count <- df %>%
    group_by(Subtype, ROI, Tissue) %>%
    summarise(Count = n())

  # Calculate cell fraction for each Subtype, ROI, and Tissue combination
  cell_fraction <- subtype_count %>%
    inner_join(cell_count, by = c("ROI", "Tissue")) %>%
    mutate(Fraction = Count / Total) %>%
    select(Subtype, ROI, Tissue, Fraction)

  # Create a data frame with unique Subtype and Tissue combinations
  roi_tissue <- unique(df[, c("ROI", "Tissue")])

  # Order ROIs based on Tissue
  roi_tissue <- roi_tissue %>%
    arrange(Tissue, ROI)

  # Create a color palette for the subtypes
  subtype_palette <- c(pal_lancet("lanonc")(7), pal_jama("default")(7), pal_jco("default")(10))

  # Create the ggplot stacked bar plot with enhanced visual style
  p <- ggplot(cell_fraction, aes(x = factor(ROI, levels = roi_tissue$ROI), y = Fraction, fill = Subtype)) +
    geom_bar(stat = "identity", color = "black", width = 1) +
    labs(x = "", y = "Cell Fraction", title = "Cell Fractions in Each ROI by Tissue") +
    scale_fill_manual(values = subtype_palette) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      legend.position = "top",
      legend.title = element_blank()
    )

  # Create the tissue color bar below the stacked bar plot
  color_bar <- ggplot(roi_tissue, aes(x = factor(ROI, levels = ROI), y = 0.5, fill = Tissue)) +
    geom_bar(stat = "identity", width = 1) +
    labs(x = "", y = "Tissue") +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      legend.position = "bottom"
    )

  # Reduce the space between the main plot and the color bar
  p <- p + theme(plot.margin = margin(t = 0, r = 0, b = -10, l = 0, unit = "pt"))
  color_bar <- color_bar + theme(plot.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "pt"))

  # Combine the stacked bar plot and the color bar using the cowplot package
  combined_plot <- p / color_bar
  combined_plot <- combined_plot + plot_layout(heights = c(14, 1))

  pdf(paste0(savePath, "Abundance Barplot from tissues.pdf"), width = 15, height = 8)
  print(combined_plot)
  dev.off()

  return(NULL)
}

BarPlotForCelltypeCounts <- function(sce, tissueCol, groupCol, savePath) {
  ## Get meta data
  meta <- colData(sce)

  countdf <- Transform_CellCountMat(sceobj = sce, group = c("IM", "CT", "TAT"), clinicalFeatures = c("ID", tissueCol, groupCol), is.fraction = T)

  plotdf <- pivot_longer(countdf, cols = 1:20, names_to = "Subpopulation", values_to = "Fraction")

  plotdf <- plotdf[, -c(1:2)]
  plotdf2 <- plotdf %>%
    group_by(Tissue, RFS_status, Subpopulation) %>%
    summarise(across(c(1:(ncol(plotdf) - 3)), mean, na.rm = TRUE))

  plotdf2$MajorType <- meta[match(plotdf2$Subpopulation, meta$SubType), "MajorType"]

  # Create a list of colors
  colors <- ggsci::pal_ucscgb("default")(26)

  plotdf2 <- as.data.frame(plotdf2)

  plotdf2[, 1] <- as.character(plotdf2[, 1])
  plotdf2[, 2] <- as.factor(plotdf2[, 2])
  plotdf2[, 3] <- as.character(plotdf2[, 3])
  plotdf2[, 4] <- as.numeric(plotdf2[, 4])
  plotdf2[, 5] <- as.character(plotdf2[, 5])

  # Create the bar plot
  p <- ggplot(plotdf2, aes(x = RFS_status, y = Fraction, fill = Subpopulation)) +
    # Add bars
    geom_bar(stat = "identity") +
    # Use different sets of colors for different MajorTypes
    scale_fill_manual(values = colors) +
    # Rotate the x-axis labels to make them readable
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    # Add title and labels
    labs(
      title = "Fraction of each Subpopulation by RFS_status and Tissue",
      x = "RFS Status", y = "Fraction", fill = "Subpopulation"
    ) +
    # Flip the coordinates
    coord_flip() +
    # Separate the plot by Tissue and MajorType
    facet_grid(MajorType ~ Tissue, scales = "free") +
    # Hide the borders of the facet grid, make the background of the facet to be empty
    theme(
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    )

  pdf(paste0("Abundance Barplot of Relapse.pdf"), height = 7.5, width = 10)
  print(p)
  dev.off()

  return(NULL)
}

## Plot the density dotplot of certain 2 channles
PlotDensityDotplot <- function(sce, marker1, marker2, MajorType = NULL, sampleSize = 100000, savePath) {
  if (!is.null(MajorType)) {
    sce <- sce[, sce$MajorType == MajorType]
  }
  exp <- assay(sce)
  exp <- exp[c(marker1, marker2), ]

  set.seed(619)
  exp_ <- exp[, sample.int(ncol(exp), sampleSize, replace = F)]
  exp_ <- as.data.frame(t(exp_))
  colnames(exp_) <- c("x", "y")

  p <- ggplot(data = exp_, aes(x = x, y = y)) +
    geom_pointdensity(size = 0.3) +
    scale_color_viridis() +
    # geom_smooth(method = lm) +
    stat_cor(method = "spearman", label.x = 0.5, label.y = 0.9) +
    xlab(marker1) +
    ylab(marker2) +
    theme(axis.title.x = element_text(size = 16, face = "bold", vjust = 0.5, hjust = 0.5)) +
    theme(axis.title.y = element_text(size = 16, face = "bold", vjust = 0.5, hjust = 0.5)) +
    theme_bw() +
    theme(
      panel.grid.major = element_line(colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.minor = element_blank(),
      text = element_text(size = 12, family = "serif")
    ) +
    theme(legend.position = "none")

  pdf(paste0(savePath, marker1, "_", marker2, " in ", MajorType, " densityDotplot.pdf"), width = 4, height = 4)
  print(p)
  dev.off()

  return(NULL)
}

## Get mean expression profile
GetMeanExpressionProfile <- function(sce, LevelCol, markers, clinicalGroupCol) {
  ids <- names(table(colData(sce)[, LevelCol]))

  df <- matrix(data = NA, nrow = length(ids), ncol = length(markers))
  clinicalMat <- matrix(data = NA, nrow = 0, ncol = length(clinicalGroupCol))
  for (i in 1:length(ids)) {
    idx <- ids[i]
    ## Subset the sce object
    sceTemp <- sce[, colData(sce)[, LevelCol] %in% idx]
    numCell <- ncol(sceTemp)
    ## Get the expression profile
    expTemp <- assay(sceTemp)
    expTemp <- expTemp[match(markers, rownames(expTemp)), ]
    df[i, ] <- apply(expTemp, MARGIN = 1, FUN = "mean")
    clinicalTemp <- colData(sceTemp)[1, match(clinicalGroupCol, colnames(colData(sceTemp)))]
    clinicalMat <- rbind(clinicalMat, as.matrix(clinicalTemp))
  }

  ## rename
  rownames(df) <- ids
  colnames(df) <- markers
  colnames(clinicalMat) <- clinicalGroupCol
  rownames(clinicalMat) <- ids

  ## return list
  list_ <- list()
  list_[["df"]] <- df
  list_[["clinicalGroup"]] <- clinicalMat
  return(list_)
}
