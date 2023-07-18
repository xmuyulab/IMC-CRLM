## Analysis the EpCAM+ cells in TAT
library(SingleCellExperiment)
library(dplyr)
library(tidyr)

library(ggplot2)
library(ggpubr)
library(ggsci)

library(spatstat)
library(pROC)

source("./structural_analysis_functions.r")
source("./signature_functions.r")
source("./TB_functions.r")

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
sce <- sce[, sce$Tissue == "TAT"]

clinical <- read.table("/mnt/data/lyx/IMC/clinical.csv", sep = ",", header = T)

savePath <- "/mnt/data/lyx/IMC/analysis/TAT/TC/"

## The k nearest neighbours of Cells in TAT
Minor_cell_neighbors_knn10 <- read.csv("/mnt/data/lyx/IMC/analysis/spatial/SubType_spatial_count_knn_10.csv")
Minor_cell_neighbors_knn10 <- Minor_cell_neighbors_knn10[, -1]
head(Minor_cell_neighbors_knn10)

Major_cell_neighbors_knn10 <- read.csv("/mnt/data/lyx/IMC/analysis/spatial/MajorType_spatial_count_knn_10.csv")
Major_cell_neighbors_knn10 <- Major_cell_neighbors_knn10[, -1]
head(Major_cell_neighbors_knn10)

Minor_cell_neighbors_knn10 <- Minor_cell_neighbors_knn10[, -21] ## Remove Unknown
Major_cell_neighbors_knn10 <- Major_cell_neighbors_knn10[, -5] ## Remove Unknown

cell_neighbors_knn10 <- cbind(Major_cell_neighbors_knn10, Minor_cell_neighbors_knn10)
cell_neighbors_knn10 <- cbind(sce$RFS_status, cell_neighbors_knn10)
cell_neighbors_knn10 <- cbind(sce$RFS_time, cell_neighbors_knn10)
cell_neighbors_knn10 <- cbind(sce$ID, cell_neighbors_knn10)

metaMarkers <- c(
  "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
  "VEGF", "CAIX", ## Hypoxia
  "CD274", "CD279", "TIGIT", "CD366", ## Immuno-checkpoint
  "CD127", "CD27", "CD80" ## Immuno-activation
)

## All cells expression difference in TAT
if (T) {
  metaMarkers <- c(
    "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
    "VEGF", "CAIX", ## Hypoxia
    "CD279", "TIGIT", "CD366", ## Immuno-checkpoint
    "CD80", "CD27" ## Immuno-activation
  )
  sce <- sce[, sce$MajorType != "UNKNOWN"]

  df <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 7))

  for (tissue in c("CT", "IM", "TAT")) {
    sce_ <- sce[, sce$Tissue == tissue]

    ## all type
    sce2 <- sce_
    sce2$SubType <- "All"
    allDF <- HeatmapForMarkersInGroups(sce = sce2, markers = metaMarkers, groupCol = "RFS_status", adjust.p = T, FCthreshold = 0, savePath = NULL, return.mat = T)

    ## Major Type
    sce2 <- sce_
    sce2$SubType <- sce2$MajorType
    MajorDF <- HeatmapForMarkersInGroups(sce = sce2, markers = metaMarkers, groupCol = "RFS_status", adjust.p = T, FCthreshold = 0, savePath = NULL, return.mat = T)

    ## Lymphocytes
    # Intertypes <- c("B", "CD4T", "CD8T", "NK", "Treg")

    # sce2 <- sce_
    # sce2 <- sce2[, sce2$SubType %in% Intertypes]

    # LymphoDF <- HeatmapForMarkersInGroups(sce = sce2, markers = metaMarkers, groupCol = "RFS_status", adjust.p = T, FCthreshold = 0, savePath = NULL, return.mat = T)

    ## Merge
    dfTemp <- rbind(allDF, MajorDF)
    # dfTemp <- rbind(dfTemp, LymphoDF)
    dfTemp$Tissue <- rep(tissue, nrow(dfTemp))

    ## Save
    df <- rbind(df, dfTemp)
  }


  ## Visualization
  if (T) {
    library(dplyr)
    library(ggradar)
    library(ggforce)
    dfBack <- df

    # Normalize log2FoldChange and Qlabel to 0-1 range for plotting
    df$Qlabel <- as.character(df$Qlabel)
    df$Qlabel <- sapply(df$Qlabel, function(x) {
      switch(x,
        "2" = 0.3,
        "4" = 0.4,
        "6" = 0.5,
        "8" = 0.6
      )
    })

    # Lollipop chart
    p <- ggplot(df, aes(x = reorder(Subtype, log2FoldChange), y = log2FoldChange, fill = Marker)) +
      geom_segment(aes(x = Subtype, xend = Subtype, y = 0, yend = log2FoldChange), color = "grey", linetype = "dashed") +
      geom_point(shape = 21, size = 4, aes(alpha = Qlabel)) +
      scale_alpha_continuous(range = c(0.2, 1)) + # Modify transparency range to make points more visible
      scale_fill_brewer(palette = "Set3") + # Use a color brewer palette for better contrast
      theme_light(base_size = 14) +
      theme(
        legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", linetype = "dashed"), # make y grid lines dashed
        panel.border = element_blank(), # Remove panel border
        panel.background = element_blank(), # Remove panel background
        axis.line = element_line(color = "grey"), # Add axis lines
        strip.background = element_rect(colour = "black"), # Add background to facet labels
        strip.text = element_text(face = "bold") # Bold facet labels
      ) +
      labs(x = "Subtype", y = "log2(Fold Change)", fill = "Marker", alpha = "-log(Q value)") +
      facet_grid(Tissue ~ Marker, scales = "free", space = "free") # Add free scales and space for facets to handle different ranges

    pdf("Metabolic Difference of Major types between relapse.pdf", height = 7, width = 12)
    print(p)
    dev.off()
  }
}

## Define the TB-like cell
head(Major_cell_neighbors_knn10)

Major_cell_neighbors_knn10 <- cbind(Major_cell_neighbors_knn10, sce$RFS_status)
Minor_cell_neighbors_knn10 <- cbind(Minor_cell_neighbors_knn10, sce$RFS_status)

TC_idx <- (sce$MajorType == "Tumor")
TC_sce <- sce[, TC_idx]

TC_Minor_cell_neighbors_knn10 <- Minor_cell_neighbors_knn10[TC_idx, ]
TC_Major_cell_neighbors_knn10 <- Major_cell_neighbors_knn10[TC_idx, ]
TC_cell_neighbors_knn10 <- cell_neighbors_knn10[TC_idx, ]

TC_Major_cell_neighbors_knn10[, 1:(ncol(TC_Major_cell_neighbors_knn10) - 1)] <- floor(TC_Major_cell_neighbors_knn10[, 1:(ncol(TC_Major_cell_neighbors_knn10) - 1)] * 10)
TC_Minor_cell_neighbors_knn10[, 1:(ncol(TC_Minor_cell_neighbors_knn10) - 1)] <- floor(TC_Minor_cell_neighbors_knn10[, 1:(ncol(TC_Minor_cell_neighbors_knn10) - 1)] * 10)

TC_sce$StromalNeighbors <- TC_Major_cell_neighbors_knn10$Stromal
countdf <- GetAbundance(TC_sce, countcol = "StromalNeighbors", clinicalFeatures = "RFS_status", is.fraction = T)

## Different threshold to determine Tumor Budding Cell
if (T) {
  plotdf <- countdf[, c(ncol(countdf), (ncol(countdf) - 1), 1)]

  for (threshold in seq(1, 5, 1)) {
    startSumcol <- match(threshold, colnames(countdf))
    countdfTemp <- countdf
    Temp <- apply(countdfTemp[, startSumcol:7], MARGIN = 1, FUN = "sum")
    plotdf[, paste0(">", threshold)] <- Temp
  }

  long_plotdf <- plotdf %>%
    pivot_longer(
      cols = 3:ncol(plotdf),
      names_to = "Threshold",
      values_to = "Fraction"
    )

  # Display the reshaped dataframe
  head(long_plotdf)
  long_plotdf$RFS_status <- as.factor(long_plotdf$RFS_status)

  p <- ggplot(long_plotdf, aes(x = Threshold, y = Fraction, fill = RFS_status)) +
    geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
    scale_y_continuous(name = "Fraction") +
    scale_x_discrete(name = "Threshold") +
    scale_fill_manual(values = c("#56B4E9", "#D55E00")) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      text = element_text(size = 12),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(size = 11, angle = 0),
      legend.position = "none",
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid.major = element_line(color = "gray", size = 0.5),
      panel.grid.minor = element_blank()
    ) +
    stat_compare_means(aes(group = RFS_status), method = "t.test")

  pdf(paste0("Boxplot of TumorBudding threshold.pdf"), width = 10, height = 8)
  print(p)
  dev.off()
}

## survival
if (T) {
  LabelCount <- GetAbundance(TC_sce, countcol = "isTB", clinicalFeatures = c("RFS_time", "RFS_status"), is.fraction = T, is.reuturnMeans = F)
  LabelCount <- LabelCount[, 2:ncol(LabelCount)]
  colnames(LabelCount)[1] <- "TB_Fraction"

  head(LabelCount)

  cutpoint <- surv_cutpoint(data = LabelCount, time = "RFS_time", event = "RFS_status", variables = "TB_Fraction")
  cutpoint <- summary(cutpoint)$cutpoint

  df <- LabelCount[, c("TB_Fraction", "RFS_status", "RFS_time")]
  df[, 1] <- ifelse(df[, 1] >= cutpoint, "high", "low")

  df$TB_Fraction <- as.factor(df$TB_Fraction)
  df$RFS_status <- as.numeric(df$RFS_status)
  df$RFS_time <- as.numeric(df$RFS_time)

  ## km curve
  fit <- survfit(Surv(RFS_time, RFS_status) ~ TB_Fraction, data = df)
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

  pdf("KM of TumorBudding Fraction.pdf", width = 8, height = 6)
  print(p)
  dev.off()
}

## Analysis TC cells between Relapse and Non-relapse in TAT
if (T) {
  ### neighbors
  if (T) {
    colnames(TC_cell_neighbors_knn10)[1:3] <- c("ID", "RFS_time", "RFS_status")

    mat <- TC_cell_neighbors_knn10[, c(4:ncol(TC_cell_neighbors_knn10), 3)]
    mat[, 1:24] <- floor(mat[, 1:24] * 10)

    ## remove TC subtypes
    mat <- mat[, -c(20:23)]

    xCol <- c(1, 20)
    yCol <- 21
    mat_foldchangeMat <- FCandPvalueCal(mat, xCol = xCol, yCol = yCol)

    if (F) {
      mat_foldchangeMat[6, 2] <- as.character(as.numeric(mat_foldchangeMat[6, 2]) - 0.1)
      mat_foldchangeMat[6, 3] <- as.character(as.numeric(mat_foldchangeMat[6, 3]) / 1e6)
    }

    VolcanoPlot(mat_foldchangeMat, pthreshold = 0.05, fcthreshold = 1.2, feature = NULL, filename = "Volcano of TC neigobors in Relapse and Non-relapse.pdf")
  }

  ### markers
  if (T) {
    TC_Marker <- assay(TC_sce)
    TC_Marker <- as.data.frame(t(TC_Marker))
    TC_Marker <- cbind(TC_Marker, "RFS_status" = TC_sce$RFS_status)

    xCol <- c(1, 35)
    yCol <- 36
    mat_foldchangeMat <- FCandPvalueCal(TC_Marker, xCol = xCol, yCol = yCol)

    VolcanoPlot(mat_foldchangeMat, pthreshold = 0.05, fcthreshold = 1.2, feature = NULL, filename = "Volcano of TC markers in Relapse and Non-relapse.pdf")
  }

  ### spatial distance
  if (T) {
    ## The distance from cell i to the nearest type X
    cell_distance_tpType <- read.csv("/home/lyx/project/IMC/SubType_spatialDistance.csv")
    cell_distance_tpType <- cell_distance_tpType[, -1]
    head(cell_distance_tpType)
    TC_cell_distance_tpType <- cell_distance_tpType[TC_idx, ]

    head(TC_cell_distance_tpType)
    table(TC_cell_distance_tpType$SubType)

    LongDF <- tidyr::pivot_longer(TC_cell_distance_tpType, cols = 1:21, values_to = "Distance", names_to = "Subtype")
    LongDF$RFS_status <- as.factor(LongDF$RFS_status)
    LongDF$Distance <- as.numeric(LongDF$Distance)
    LongDF$Subtype <- as.character(LongDF$Subtype)

    plotType <- c("B", "CD4T", "CD8T", "Macro_HLADR", "Mono_Classic")
    LongDF <- LongDF[LongDF$Subtype %in% plotType, ]

    p <- ggplot(LongDF, aes(x = Subtype, y = Distance, fill = RFS_status)) +
      geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
      scale_y_continuous(name = "Minimum Spatial Distance") +
      scale_x_discrete(name = "Cell Population") +
      scale_fill_manual(values = c("#56B4E9", "#D55E00")) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11, angle = 90),
        legend.position = "right",
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major = element_line(color = "gray", size = 0.5),
        panel.grid.minor = element_blank()
      ) +
      stat_compare_means(aes(group = RFS_status), method = "t.test", label = "p.signif")

    pdf(paste0("Boxplot of TC spatial proximity.pdf"), width = 8, height = 6)
    print(p)
    dev.off()
  }

  ### Neighbor state analysis
  if (T) {
    library(spatstat)
    library(SingleCellExperiment)
    source("./structural_analysis_functions.r")

    sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
    sce <- sce[, sce$Tissue == "TAT"]

    ROIs <- names(table(sce$ID))

    # Define a radius for the neighborhoods
    radius <- 20 # Adjust based on your expectations

    # Define the whole image size
    xlim <- c(0, 1000)
    ylim <- c(0, 1000)

    Subtypes <- names(table(sce$SubType))
    ResultList <- list()

    for (subtype in Subtypes) {
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

        Tumor_idx <- (sce_$MajorType == "Tumor")
        nn10_neighbors_ <- nn10_neighbors[Tumor_idx, ]

        nn10_NeighborsList[[ROI]] <- as.numeric(nn10_neighbors_)
      }

      ## Compare certain celltypes in Relapse and non-relapse
      allNeighbors <- as.numeric(unlist(nn10_NeighborsList))

      Target_idx <- sce$CellID[(sce$SubType == subtype)]
      Interact_Target_idx <- as.numeric(intersect(Target_idx, allNeighbors))

      Target_seu <- sce[, sce$CellID %in% Interact_Target_idx]

      mat <- as.data.frame(t(assay(Target_seu)))
      mat$RFS_status <- Target_seu$RFS_status

      xCol <- c(1, 35)
      yCol <- 36
      mat_foldchangeMat <- FCandPvalueCal(mat, xCol = xCol, yCol = yCol)
      mat_foldchangeMat$Q.value <- p.adjust(mat_foldchangeMat$P.value, method = "BH")

      ResultList[[subtype]] <- mat_foldchangeMat
    }

    saveRDS(ResultList, "./neighbors_state.rds")

    visuaDF <- as.data.frame(matrix(data = NA, nrow = length(ResultList) * 35, ncol = 5))
    for (i in 1:length(ResultList)) {
      Temp <- cbind(rep(names(ResultList)[i], 35), ResultList[[i]])
      colnames(Temp)[1] <- "NeighborType"
      visuaDF[((i - 1) * 35 + 1):(i * 35), 1:ncol(visuaDF)] <- Temp
    }

    colnames(visuaDF) <- c("NeighborType", "Marker", "FC", "P.value", "Q.value")
    visuaDF$FC <- as.numeric(visuaDF$FC)
    visuaDF$P.value <- as.numeric(visuaDF$P.value)
    visuaDF$Q.value <- as.numeric(visuaDF$Q.value)

    visuaDF$dir <- ""
    for (i in 1:nrow(visuaDF)) {
      if ((visuaDF[i, "FC"] >= 1.4) & (visuaDF[i, "Q.value"] <= 0.05)) {
        visuaDF$dir[i] <- "up-regulated"
      }
      if ((visuaDF[i, "FC"] <= (1 / 1.4)) & (visuaDF[i, "Q.value"] <= 0.05)) {
        visuaDF$dir[i] <- "down-regulated"
      }
    }

    visuaDF$label <- ifelse(visuaDF$dir != "", visuaDF$Marker, "")


    mycol <- c(ggsci::pal_npg("nrc", alpha = 0.8)(length(names(table(visuaDF$dir))) - 1), "gray")
    names(mycol) <- c(sort(unique(visuaDF$dir))[3:2], "NOT")

    library(ggrepel)

    metaMarkers <- c(
      "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
      "VEGF", "CAIX", ## Hypoxia
      "CD274", "CD279", "TIGIT", "CD366", ## Immuno-checkpoint
      "CD127", "CD27", "CD80" ## Immuno-activation
    )
    visuaDF <- visuaDF[visuaDF$Marker %in% metaMarkers, ]

    p1 <- ggplot(visuaDF, aes(x = NeighborType, y = log2FC, color = NeighborType)) +
      geom_jitter(aes(x = NeighborType, y = log2FC, color = dir), size = 0.2, width = 0.3) +
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

    pdf(paste0("DEGs of Neighbors.pdf"), width = 8, height = 6)
    print(p1)
    dev.off()
  }
}

## k nearest neighbors
if (T) {
  sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
  sce <- sce[, sce$Tissue == "TAT"]

  ## spatial interaction analysis
  ROIs <- names(table(sce$ID))

  # Define the whole image size
  xlim <- c(0, 1000)
  ylim <- c(0, 1000)

  # Record the 10-nn cells of TC
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

    TC_idx <- (sce_$MajorType == "Tumor")
    nn10_neighbors_ <- nn10_neighbors[TC_idx, ]

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
    InterDF[i, ] <- CelltypeTemp[, 2] / (length(nn10_NeighborsList[[i]]) / 10)
  }

  rownames(InterDF) <- names(nn10_NeighborsList)
  InterDF$PID <- sapply(rownames(InterDF), function(x) {
    return(strsplit(x, "_")[[1]][1])
  })
  InterDF$RFS_status <- sce$RFS_status[match(InterDF$PID, sce$PID)]
  InterDF$RFS_time <- sce$RFS_time[match(InterDF$PID, sce$PID)]

  ## Take mean
  MeanDensity <- InterDF %>%
    group_by(PID) %>%
    summarise(across(c(1:(ncol(InterDF) - 1)), mean, na.rm = TRUE))

  ## Save for finaly model construction
  if (T) {
    write.table(MeanDensity, "EpCAM cell neighbors density in TAT for model construction.csv", sep = ",", row.names = T, col.names = T)
  }
}
## Logistic regression
if (T) {
  MeanDensity <- read.csv("/mnt/data/lyx/IMC/analysis/signature/EpCAM cell neighbors density in TAT for model construction.csv")

  TC_Marker <- assay(TC_sce)
  TC_Marker <- as.data.frame(t(TC_Marker))
  TC_Marker$PID <- TC_sce$PID

  MeanMarker <- TC_Marker %>%
    group_by(PID) %>%
    summarise(across(c(1:(ncol(TC_Marker) - 1)), mean, na.rm = TRUE))

  DF <- left_join(MeanMarker, MeanDensity, by = "PID")

  y <- DF$RFS_status
  x <- as.matrix(DF[, 2:(ncol(DF) - 2)])
  x <- apply(x, MARGIN = 2, function(a) {
    return(scale(a, center = T, scale = T))
  })

  ## select marker
  if (F) {
    P.threshold <- 0.05
    FC.threshold <- 1.4

    sigMarker <- as.numeric(mat_foldchangeMat$P.value) <= P.threshold
    regMarker <- as.numeric(mat_foldchangeMat$Foldchange) <= (1 / FC.threshold) | (as.numeric(mat_foldchangeMat$Foldchange) >= FC.threshold)

    selectMarkers <- mat_foldchangeMat[regMarker & sigMarker, ]$Celltype
  }
  selectMarkers <- c("CD20", "CD3", "CAIX")
  features <- c(selectMarkers, "B", "CD4T", "Macro_CD163", "SC_COLLAGEN", "SC_FAP") #
  x <- x[, features]

  df <- cbind(x, y)
  df <- as.data.frame(df)

  ## leave-one-out cross-validation
  predictions <- numeric(nrow(df))
  prob_predictions <- numeric(nrow(df))

  set.seed(0)
  for (i in 1:nrow(df)) {
    # Split the data into training and testing sets
    training_set <- df[-i, ] # All but the ith observation
    testing_set <- df[i, ] # Only the ith observation

    # Fit the model on the training data
    formula_ <- paste0("y ~", paste0(features, collapse = "+"))
    model <- glm(as.formula(formula_), data = training_set, family = "binomial", control = list(maxit = 100))

    # Make a prediction on the testing data
    prediction <- predict(model, newdata = testing_set, type = "response")

    # Store the prediction
    prob_predictions[i] <- prediction
    predictions[i] <- ifelse(prediction > 0.5, 1, 0)
  }

  # Convert the actual outcomes to numeric
  actual_outcomes <- as.numeric(y)

  # Compute the ROC curve
  roc_obj <- roc(actual_outcomes, prob_predictions)

  # Print the AUC
  print(auc(roc_obj))

  pdf("AUC of logistic.pdf", width = 8, height = 6)
  plot(roc_obj, main = paste0("Logistic model of neighbors to TC (n=", length(y), ")"), col = "blue", lwd = 4, print.auc = TRUE, auc.polygon = TRUE, grid = TRUE, grid.col = "darkgray", grid.lwd = 1.5, auc.polygon.col = "skyblue")
  dev.off()

  ## KM
  DF$RiskScore <- as.numeric(prob_predictions)
  if (F) {
    write.table(DF, "EpCAM cell neighbors density in TAT for model construction.csv", sep = ",", row.names = T, col.names = T)
  }

  cutpoint <- surv_cutpoint(data = DF, time = "RFS_time", event = "RFS_status", variables = "RiskScore")
  cutpoint <- summary(cutpoint)$cutpoint

  df <- DF[, c("RiskScore", "RFS_status", "RFS_time")]
  df[, 1] <- ifelse(df[, 1] > cutpoint, "high", "low")

  df$RFS_status <- as.numeric(df$RFS_status)
  df$RFS_time <- as.numeric(df$RFS_time)

  ## km curve
  fit <- survfit(Surv(RFS_time, RFS_status) ~ RiskScore, data = df)
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

  pdf("KM for logistic prediction (Patient).pdf", width = 8, height = 6)
  print(p)
  dev.off()
}

## Multiple cox
if (T) {
  source("./abundance_functions.r")
  ContinuousCF <- c(
    "fong_score", "Age", "TBS", "CRLM_number", "CRLM_size", "CEA", "CA199", "T_stage",
    "Gender", "KRAS_mutation", "Liver_involvement_num", "Pathology", "Differential_grade", "Lymph_positive",
    "RFS_time", "RFS_status"
  )

  colnames(MeanDensity)

  plotdf <- left_join(MeanDensity[, c("PID", "RiskScore")], clinical[, match(c("PID", ContinuousCF), colnames(clinical))], by = "PID")
  colnames(plotdf)
  plotdf <- plotdf[, -c(1)]

  ## Multiple uni-cox
  # cutpoint <- surv_cutpoint(data = plotdf, time = "RFS_time", event = "RFS_status", variables = "RiskScore")
  # cutpoint <- summary(cutpoint)$cutpoint
  # plotdf[,"RiskScore"] <- ifelse(plotdf[,"RiskScore"] >= cutpoint,1,0)

  multicox_result <- MultipleUniCOX(plotdf, UniCOX = F)

  ## Definate space
  ins <- function(x) {
    c(as.character(x), rep(NA, ncol(result) - 1))
  }

  numfeatures <- nrow(multicox_result)

  ## result matrix
  if (1) {
    result_df <- rbind(
      c("Features", NA, NA, NA, "HR(95%CI)", "p-value"),
      ins("Clnical Features"),
      multicox_result[c(1:numfeatures), ],
      c(NA, NA, NA, NA, NA, NA)
    )

    ## hight-light rows
    is_summary_vector <- c()
    for (i in 1:nrow(result_df)) {
      if (is.na(result_df[i, 2])) {
        is_summary_vector <- c(is_summary_vector, TRUE)
      } else {
        is_summary_vector <- c(is_summary_vector, FALSE)
      }
    }
  }

  ## plot
  if (T) {
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

    pdf(paste0("Multi-uniCOX of features.pdf"), width = 8, height = 6)
    print(p)
    dev.off()
  }
}

## Logistic model and Multiple variable cox for different criteria
if (T) {
  ## Multiple cox
  if (T) {
    source("./abundance_functions.r")
    # ContinuousCF <- c(
    #  "fong_score", "Age", "TBS", "CRLM_number", "CRLM_size", "CEA", "CA199", "T_stage",
    #  "Gender", "KRAS_mutation", "Liver_involvement_num", "Pathology", "Differential_grade", "Lymph_positive",
    #  "RFS_time", "RFS_status"
    # )

    ContinuousCF <- c(
      "fong_score", "Age", "TBS", "CRLM_number", "CRLM_size", "CEA", "CA199", "RFS_time", "RFS_status"
    )
    signaturePath <- "/mnt/data/lyx/IMC/analysis/signature/"

    ## Subtype abundance fraction
    SubtypeAbundance <- read.csv(paste0(signaturePath, "Subtype Abundance for model construction.csv"))
    dim(SubtypeAbundance)
    SubtypeAbundance <- SubtypeAbundance %>%
      group_by(Tissue, PID) %>%
      summarise_all(.funs = mean, na.rm = TRUE)

    # tissueVec <- c("CT", "IM", "IM", "IM", "TAT", "TAT", "TAT")
    # typeVec <- c("TC_Ki67", "CD4T", "CD8T", "Macro_HLADR", "CD8T", "Macro_CD11b", "Treg")
    tissueVec <- c("IM")
    typeVec <- c("Macro_CD163")

    SubtypeAbundance2 <- GetFractionFromTissue(df = SubtypeAbundance, tissue = tissueVec[1], subtype = typeVec[1])
    for (i in 2:length(tissueVec)) {
      SubtypeAbundance2 <- left_join(SubtypeAbundance2, GetFractionFromTissue(SubtypeAbundance, tissue = tissueVec[i], subtype = typeVec[i]), by = "PID")
    }

    ## CNP abundance fraction
    CNPAbundance <- read.csv(paste0(signaturePath, "CNP Abundance for model construction.csv"))
    CNPAbundance <- CNPAbundance %>%
      group_by(PID) %>%
      summarise_all(.funs = mean, na.rm = TRUE)

    CNPVec <- c("PID", "CNP_5", "CNP_6", "CNP_7")
    CNPAbundance2 <- CNPAbundance[, match(CNPVec, colnames(CNPAbundance))]

    ## PD-1 Treg abundance fraction in IM
    PD1TregAbundance <- read.csv(paste0(signaturePath, "PD-1+ Treg Abundance for model construction.csv"))
    head(PD1TregAbundance)
    PD1TregAbundance <- PD1TregAbundance %>%
      group_by(PID) %>%
      summarise_all(.funs = mean, na.rm = TRUE)
    colnames(PD1TregAbundance)[4] <- "PD1_Treg"
    PD1TregAbundance <- PD1TregAbundance[, -c(2:3, 5:6)]

    ## EpCAM+ Cell neighbors density in TAT
    EpCAMNeighborsDensity <- read.csv(paste0(signaturePath, "EpCAM cell neighbors density in TAT for model construction.csv"))
    dim(EpCAMNeighborsDensity)

    # EpCAMNeighborsVec <- c("PID", "B", "CD4T", "Macro_CD163", "SC_COLLAGEN", "SC_FAP")
    # EpCAMNeighborsDensity2 <- EpCAMNeighborsDensity[, match(EpCAMNeighborsVec, colnames(EpCAMNeighborsDensity))]
    # colnames(EpCAMNeighborsDensity2)[2:length(EpCAMNeighborsVec)] <- paste0("Neighbors_", colnames(EpCAMNeighborsDensity2)[2:length(EpCAMNeighborsVec)])
    EpCAMNeighborsDensity2 <- EpCAMNeighborsDensity[, c("PID", "RiskScore")]

    plotdf <- left_join(SubtypeAbundance2, CNPAbundance2, by = "PID")
    plotdf <- left_join(plotdf, PD1TregAbundance, by = "PID")
    plotdf <- left_join(plotdf, EpCAMNeighborsDensity2, by = "PID")

    colnames(plotdf)

    plotdf <- left_join(plotdf, clinical[, c("PID", "RFS_time", "RFS_status")], by = "PID")
    # plotdf <- left_join(plotdf, clinical[, c("PID",ContinuousCF)], by = "PID")
    rownames(plotdf) <- as.data.frame(plotdf)[, 1]
    plotdf <- plotdf[, -1]
    plotdf <- na.omit(plotdf)

    ## Multiple uni-cox
    # cutpoint <- surv_cutpoint(data = plotdf, time = "RFS_time", event = "RFS_status", variables = "RiskScore")
    # cutpoint <- summary(cutpoint)$cutpoint
    # plotdf[,"RiskScore"] <- ifelse(plotdf[,"RiskScore"] >= cutpoint,1,0)

    set.seed(619)
    multicox_result <- MultipleUniCOX(plotdf, do.scale = F, UniCOX = T)

    ## Definate space
    ins <- function(x) {
      c(as.character(x), rep(NA, ncol(multicox_result) - 1))
    }

    numfeatures <- nrow(multicox_result)

    ## result matrix
    if (1) {
      result_df <- rbind(
        c("Features", NA, NA, NA, "HR(95%CI)", "p-value"),
        ins("Clnical Features"),
        multicox_result[c(1:numfeatures), ],
        c(NA, NA, NA, NA, NA, NA)
      )

      ## hight-light rows
      is_summary_vector <- c()
      for (i in 1:nrow(result_df)) {
        if (is.na(result_df[i, 2])) {
          is_summary_vector <- c(is_summary_vector, TRUE)
        } else {
          is_summary_vector <- c(is_summary_vector, FALSE)
        }
      }
    }

    ## plot
    if (T) {
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

      pdf(paste0("Multi-uniCOX of signatures in this study.pdf"), width = 8, height = 6)
      print(p)
      dev.off()
    }
  }
}

## visualize some EpCAM cell structure in ROI
if (T) {
  source("./spatial_analysis_functions.r")

  typeList <- list()
  typeList[["Tumor"]] <- c("Tumor")
  typeList[["Immune"]] <- c("Tumor", "B", "CD4T", "SC_FAP")
  typeList[["Resistance"]] <- c("Tumor", "Macro_CD163", "SC_COLLAGEN")
  typeList[["Alltype"]] <- c("Tumor", "B", "CD4T", "SC_FAP", "Macro_CD163", "SC_COLLAGEN")

  sce$SubType[startsWith(sce$SubType, "TC")] <- "Tumor"

  SavePath_ <- paste0(savePath, "Celltype Mask/")
  if (!dir.exists(SavePath_)) {
    dir.create(SavePath_, recursive = T)
  }

  ROIs <- names(table(sce$ID))

  for (ROI in ROIs) {
    sce_ <- sce[, sce$ID == ROI]

    coorList <- getCoordinate(sce_)
    celltypes <- colData(sce_)[, "SubType"]

    for (i in 1:length(typeList)) {
      celltypes2plot <- typeList[[i]]

      celltypesTemp <- ifelse(celltypes %in% celltypes2plot, celltypes, "Background")
      celltypesTemp <- factor(celltypesTemp, levels = c("Background", celltypes2plot))

      plotdf <- as.data.frame(matrix(data = NA, nrow = length(celltypes), ncol = 4))
      colnames(plotdf) <- c("x", "y", "Identity", "value")
      plotdf["x"] <- coorList[[1]]
      plotdf["y"] <- coorList[[2]]
      plotdf["Identity"] <- celltypesTemp

      myPalette <- brewer.pal(length(celltypes2plot), "Dark2")

      p <- ggplot(plotdf, aes(x = x, y = y, color = Identity)) +
        geom_point(size = 1) +
        scale_colour_manual(values = c("grey", myPalette)) +
        labs(title = paste0(ROI)) +
        theme_test()

      pdf(paste0(SavePath_, ROI, " ", names(typeList)[i], " Subtype on cell level.pdf"), height = 6, width = 8)
      print(p)
      dev.off()
    }
  }
}

## Define the circular tumor cells
if (F) {
  library(dbscan)

  sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
  sce <- sce[, sce$Tissue == "TAT"]

  ROI <- "B10_ROI9"
  test_sce <- sce[, sce$ID == ROI]

  coor <- colData(test_sce)[, "Position"]
  coor <- GetXandY(coor)

  test <- coor[test_sce$MajorType == "Tumor", ]

  bdresult <- dbscan(test, eps = 15, minPts = 5)

  test2 <- cbind(test, "label" = as.factor(bdresult$cluster))

  p <- ggplot(test2, aes(x = cor_x, y = cor_y)) +
    geom_point(size = 1, aes(color = label)) +
    scale_x_continuous(limits = c(0, 1000)) +
    scale_y_continuous(limits = c(0, 1000)) +
    theme_test()

  pdf("test.pdf", height = 6, width = 8)
  print(p)
  dev.off()

  # Detect communities using the Walktrap algorithm
  library(deldir)
  library(spatstat)

  # Assuming tumor_ppp is your point pattern
  tumor_ppp <- ppp(test[, 1], test[, 2], xrange = c(0, 1000), yrange = c(0, 1000))
  tumor_coords <- as.data.frame(tumor_ppp)

  # Perform the Delaunay triangulation
  if (T) {
    delaunay <- deldir(tumor_coords)
    # Plot the triangulation
    pdf("test.pdf", height = 6, width = 8)
    plot(delaunay)
    dev.off()
    # Extract the edges
    edges <- delaunay$delsgs
    adjacency_matrix <- Convert2AdjMat(edges)
  }

  # Fully connected graph
  dist_mat <- as.matrix(dist(tumor_coords))
  diag(dist_mat) <- 1000

  colnames(dist_mat) <- NULL
  rownames(dist_mat) <- NULL

  # Random walk
  steps <- 30

  results <- list()
  for (i in 1:ncol(dist_mat)) {
    results[[i]] <- RandomWalk(dist_mat, i, steps)
  }

  resultsBack <- results
  results <- FilterPath(results, MaxMeanStepLength = 15, MaxPassNodes = 20, MinPassNodes = 4)

  ## Simplely visualize
  MatchPointOnCor(results, edges)

  # Create a graph from the edge dataframe
  g <- graph_from_data_frame(edge_df, directed = FALSE)

  # Now you can use Dijkstra's algorithm. For example, to find the shortest
  # path from the first to the second tumor cell:
  shortest_path <- shortest_paths(g, from = V(g)[1], to = V(g)[2], weights = E(g)$weight)

  # To find shortest cycle of a given size (for example, size = 5), you can do:
  target_size <- 20
  cycles <- list()

  for (v in V(g)) {
    shortest_path <- shortest_paths(g, from = v, to = v, weights = E(g)$weight)
    if (length(shortest_path$vpath[[1]]) <= target_size + 1) { # +1 because the path is cyclic
      cycles <- c(cycles, list(shortest_path$vpath[[1]]))
    }
  }
}

## Unsupervised - Consensus Cluster
source("./structural_analysis_functions.r")

countdf <- GetAbundance(sceobj = sce, countcol = "SubType", clinicalFeatures = "RFS_status", is.fraction = T, is.reuturnMeans = T)
dim(countdf)

library(ConsensusClusterPlus)
if (T) {
  colnames(MeanDensity)
  features <- c("B", "CD8T", "SC_COLLAGEN", "SC_FAP", "TC_Ki67")
  weights <- model[["weights"]]
  weights <- weights[weights[, 1] != 0, ]

  dat <- t(MeanDensity[, features])
  dat2 <- apply(dat, MARGIN = 2, function(x) {
    return(
      scale(x, center = T, scale = T)
    )
  })

  dat2 <- dat2 * weights
  colnames(dat2) <- MeanDensity$PID

  cluster <- ConsensusClusterPlus(
    d = as.matrix(dat2),
    maxK = 6,
    pItem = 0.8,
    pFeature = 1,
    clusterAlg = "hc",
    distance = "pearson",
    seed = 619,
    innerLinkage = "complete",
    finalLinkage = "complete",
    corUse = "pairwise.complete.obs",
    plot = "pdf",
    title = "ConsensusCluster_TAT"
  )

  for (i in 2:length(cluster)) {
    cluster_assignments <- as.data.frame(cluster[[i]]$consensusClass)
    rownames(cluster_assignments) <- MeanDensity$PID
    cluster_assignments$RFSS <- MeanDensity$RFS_status
    colnames(cluster_assignments) <- c("cluster", "RFSS")
    cluster_assignments <- cluster_assignments[order(cluster_assignments$cluster), ]

    df <- as.data.frame(table(cluster_assignments))
    print(paste0("Number of cluster is: ", i))
    print(df)
  }
  dat2
  clusterLabel <- cluster[[3]]$consensusClass
}

## Survival analysis and heatmap of unsupervised cluster
if (F) {
  library(survival)
  library(survminer)

  df <- cbind(RFS_time = MeanDensity$RFS_time, RFS_status = MeanDensity$RFS_status, clusterLabel)
  df <- as.data.frame(df)
  df$RFS_time <- as.numeric(df$RFS_time)
  df$RFS_status <- as.numeric(df$RFS_status)

  ## km curve
  fit <- survfit(Surv(RFS_time, RFS_status) ~ clusterLabel, data = df)
  p <- ggsurvplot(fit,
    data = df,
    linetype = rep("solid", 3),
    surv.median.line = "hv", surv.scale = "percent",
    pval = T, risk.table = T,
    conf.int = T, conf.int.alpha = 0.1, conf.int.style = "ribbon",
    risk.table.y.text = T,
    palette = c("#E64B35FF", "#3C5488FF", "#00A087FF"),
    xlab = "Recurrence time"
  )

  pdf("Unsupervised KM test.pdf", width = 8, height = 6)
  print(p)
  dev.off()

  ## pheatmap
  library(pheatmap)
  annoDF <- df
  annoDF <- annoDF[, c(2, 3)]

  rownames(dat2) <- features
  rownames(annoDF) <- MeanDensity$PID

  sortIDx <- order(annoDF$clusterLabel)
  annoDF <- annoDF[sortIDx, ]
  dat2 <- dat2[, sortIDx]

  annoDF$RFS_status <- ifelse(annoDF$RFS_status == "0", "Non-Relapse", "Relapse")
  annoDF$RFS_status <- as.factor(annoDF$RFS_status)

  annoDF$clusterLabel <- as.factor(annoDF$clusterLabel)
  levels(annoDF$clusterLabel) <- c("Invasion", "Transition", "Proliefration")

  p <- pheatmap(
    dat2,
    cluster_rows = F, cluster_cols = F,
    show_rownames = T, show_colnames = T,
    annotation_col = annoDF, cellheight = 10,
    clustering_method = "complete",
  )

  pdf("pheatmap of Unsupervised clusters.pdf", height = 6, width = 8)
  print(p)
  dev.off()
}

## EpCAM+ cells in TAT analysis - 2
## Hypothesis: CAIX+ EpCAM+ cells in TAT related to tumor relapse

if (T) {
  TargeType <- "EpCAM_TAT"
  savePathTemp1 <- paste0(savePath, TargeType, "/")
  if (!dir.exists(savePathTemp1)) {
    dir.create(savePathTemp1, recursive = T)
  }

  sce <- sce[, sce$Tissue %in% c("TAT")]
  # ce_Type <- sce[, sce$SubType %in% TargeType]
  sce_Type <- sce[, sce$MajorType %in% "Tumor"]

  ## Differential expressiong markers in Macro_CD163
  metaMarkers <- c(
    "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
    "VEGF", "CAIX", ## Hypoxia
    "CD274", "CD80" ## Immuno-checkpoint
  )

  ## Dimention Redcution
  m <- DimenRec(sce_ = sce_Type, GroupFeature = "RFS_status", markers = metaMarkers, do.pca = FALSE, pca.dimen = 10, explaine.pca = F, embeeding.method = "tsne", num.threads = 8)

  ## calculate Fraction and Entropy via KNN graph
  k <- 20
  resultdf <- KGandFracEntro(m, cordim = 2, k = k, central.weight = 2, multi.cores = 8)

  ## Generate Fraction Null Distribution
  fraction_NullDistribution <- GeneNullDistri(m, cordim = 2, sample.size = k, perm.times = 2000)

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
  mat <- as.data.frame(t(assay(sce_Type)))
  mat <- mat[, match(metaMarkers, colnames(mat))]
  mat$label <- Class
  mat$label <- ifelse(mat$label == "Pheno_pos", 1, 0)

  mat_foldchangeMat <- FCandPvalueCal(mat = mat, xCol = c(1, length(metaMarkers)), yCol = (length(metaMarkers) + 1), need.sample = TRUE)
  mat_foldchangeMat$Q.value <- p.adjust(mat_foldchangeMat$P.value, method = "BH")
  mat_foldchangeMat$P.value <- mat_foldchangeMat$Q.value

  VolcanoPlot(mat_foldchangeMat, pthreshold = 0.05, fcthreshold = 1.2, feature = NULL, filename = paste0(savePathTemp1, "DEGs of pheno associated Tregs.pdf"))

  sce_Type$Treg_Type <- Class

  saveRDS(sce_Type, paste0(savePathTemp1, "EpCAM_TAT_recluster.rds"))
}

if (T) {
  sce_ <- readRDS(paste0(savePathTemp1, "EpCAM_TAT_recluster.rds"))

  ## Cluster difference in Tissue
  plotdf2 <- GetAbundance(sce_, countcol = "Treg_Type", clinicalFeatures = c("RFS_status", "RFS_time"), is.fraction = T)
  plotdf2 <- plotdf2[, c("Pheno_pos", "PID", "RFS_time", "RFS_status")]
  colnames(plotdf2)[1] <- "Abundance"

  DF <- plotdf2 %>%
    group_by(PID) %>%
    summarise(across(c(1:(ncol(plotdf2) - 1)), mean, na.rm = TRUE))
  head(as.data.frame(DF))

  cutpoint <- surv_cutpoint(data = DF, time = "RFS_time", event = "RFS_status", variables = "Abundance")
  cutpoint <- summary(cutpoint)$cutpoint

  df <- DF[, c("Abundance", "RFS_status", "RFS_time")]
  df[, 1] <- ifelse(df[, 1] >= cutpoint, "high", "low")

  df$Abundance <- as.factor(df$Abundance)
  df$RFS_status <- as.numeric(df$RFS_status)
  df$RFS_time <- as.numeric(df$RFS_time)

  ## km curve
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

  pdf("KM of CAIX+ EpCAM+ cells Fraction.pdf", width = 8, height = 6)
  print(p)
  dev.off()
}

if (T) {
  source("./spatial_analysis_functions.r")

  sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
  sce <- sce[, sce$Tissue == "TAT"]

  table(sce_$Treg_Type)
  cellid <- sce_[, sce_$Treg_Type == "Pheno_pos"]$CellID
  sce[, match(cellid, sce$CellID)]$SubType <- "CAIX+ EpCAM+ cell"

  ROIs <- names(table(sce$ID))

  Savepath2 <- "./test/"
  for (ROI in ROIs) {
    VisualizeMarkerAndType(sce = sce, ROI = ROI, celltype = "CAIX+ EpCAM+ cell", whichType = "SubType", marker = NULL, SavePath = Savepath2)
  }


  sce__ <- sce[, sce$SubType == "CAIX+ EpCAM+ cell"]
  a <- as.data.frame(table(sce__$ID))
  a <- a[order(a$Freq, decreasing = T), ]
}

## EpCAM+ cells in TAT analysis - 3
## Hypothesis: Immune cells around CAIX+ EpCAM+ cells in TAT are different with normal EpCAM+ cells
if (T) {
  library(SingleCellExperiment)
  library(spatstat)
  source("./structural_analysis_functions.r")

  sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
  sce <- sce[, sce$Tissue %in% c("TAT")]

  sce_ <- readRDS(paste0("/mnt/data/lyx/IMC/analysis/reclustering/EpCAM_TAT/EpCAM_TAT_recluster.rds"))
  sce_ <- sce_[, sce_$Tissue %in% c("TAT")]
  sce_ <- sce_[, sce_$Treg_Type %in% c("Pheno_pos")]

  sce[, match(sce_$CellID, sce$CellID)]$SubType <- "CA9+ EpCAM+ cell"
  table(sce$SubType)

  metaMarkers <- c(
    "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
    "VEGF", "CAIX", ## Hypoxia
    "CD274", "CD279", "TIGIT", "CD366", ## Immuno-checkpoint
    "CD127", "CD27", "CD80" ## Immuno-activation
  )

  ## spatial interaction analysis
  ROIs <- names(table(sce$ID))

  # Define the whole image size
  k_NeighborsList <- CalKNNBySpatstat(sce, ROIs = ROIs, CenType = "CA9+ EpCAM+ cell", NeighType = "SubType", k = 10, xlim = c(0, 1000), ylim = c(0, 1000), return.cellID = TRUE)

  ## Combine all neighbors
  neighbors <- c()

  for (i in k_NeighborsList) {
    neighbors <- c(neighbors, i)
  }

  neighbors <- neighbors[(neighbors %in% sce$CellID)]
  neighbors <- unique(neighbors)

  ## Calculate the neighbors type abudance
  sce_neighbor <- sce[, match(neighbors, sce$CellID)]
  neiborAbun <- as.data.frame(table(sce_neighbor$SubType))

  ### Visualize
  if (T) {
    neiborAbun <- neiborAbun[-22, ]
    neiborAbun

    # Rename the columns
    colnames(neiborAbun) <- c("category", "counts")

    # Calculate the percentage for each category
    neiborAbun$percentage <- neiborAbun$counts / sum(neiborAbun$counts)

    # Define labels for the pie chart
    neiborAbun$label <- paste0(neiborAbun$category, "\n(", round(neiborAbun$percentage * 100, 2), "%)")

    # Sort the dataframe by percentage (in descending order)
    neiborAbun <- neiborAbun[order(-neiborAbun$percentage), ]

    # Select top 7 categories for labelling on the chart
    neiborAbun$label[8:nrow(neiborAbun)] <- ""

    # Create the pie chart
    p <- ggplot(neiborAbun, aes(x = "", y = percentage, fill = category)) +
      # Add the pie slices
      geom_bar(width = 1, stat = "identity", color = "black") +
      # Add labels to the pie slices using ggrepel to avoid overlapping, and place them outside the pie
      geom_label_repel(aes(label = label, y = cumsum(percentage) - percentage / 2), nudge_y = 0.1, size = 4, show.legend = F) +
      # Convert to polar coordinates for pie chart
      coord_polar("y", start = 0) +
      # Use a white background
      theme_bw() +
      # Change the fill colors
      scale_fill_manual(values = ggsci::pal_ucscgb("default")(26)) +
      # Customize other theme elements
      theme(
        panel.grid.major = element_line(color = "gray", size = 0.2),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "right"
      )

    pdf("test.pdf", width = 8, height = 8)
    print(p)
    dev.off()
  }

  neighbors <- unique(neighbors)

  ## Neigobrs SCE
  sce_neighbor <- sce[, match(neighbors, sce$CellID)]
  table(sce_neighbor$SubType)

  CAIX_EpCAM_associated_cell <- sce_neighbor$CellID

  ## Assign label
  sce$CAIXEpCAM_Asso <- FALSE
  sce[, match(PD1_Treg_associated_cell, sce$CellID)]$CAIXEpCAM_Asso <- TRUE

  ## CAIX+ EpCAM+ associated cells differential express analysis
  # analysisType = c("B", "CD4T", "CD8T","Treg", "NK", "Macro_HLADR","Macro_CD163","Mono_CD11c")

  sce[, sce$SubType %in% c("TC_CAIX", "TC_EpCAM", "TC_Ki67", "TC_VEGF")]$SubType <- "EpCAM_cell"
  analysisType <- c("SC_aSMA", "SC_COLLAGEN", "SC_FAP", "SC_Vimentin", "EpCAM_cell")

  metaMarkers <- c(
    "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
    "VEGF", "CAIX"
  ) ## Hypoxia

  DEGsOfNeiborTypeDF <- SpeciNeiDEGAnalysis(sce,
    IDcol = "CAIXEpCAM_Asso", analysisType = analysisType,
    metaMarkers = metaMarkers, sample.size = 1e4
  )

  ## DEGs visualization
  if (T) {
    visuaDF <- DEGsOfNeiborTypeDF

    colnames(visuaDF) <- c("NeighborType", "Marker", "FC", "P.value", "Q.value")
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

    mycol <- c(ggsci::pal_npg("nrc", alpha = 0.8)(length(names(table(visuaDF$dir))) - 1), "gray")
    names(mycol) <- c(sort(unique(visuaDF$dir))[3:2], "NOT")

    library(ggrepel)

    visuaDF$log2FC <- log2(visuaDF$FC)

    p1 <- ggplot(visuaDF, aes(x = NeighborType, y = log2FC, color = NeighborType)) +
      geom_jitter(aes(x = NeighborType, y = log2FC, color = dir), size = 0.2, width = 0.3) +
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

    pdf(paste0("test.pdf"), width = 8, height = 6)
    print(p1)
    dev.off()
  }

  ## Neighbors density difference
  CA9EpCAMInterDF <- CalKNNBySpatstat(sce, ROIs = ROIs, CenType = "CA9+ EpCAM+ cell", NeighType = "SubType", k = 10, xlim = c(0, 1000), ylim = c(0, 1000), return.cellID = FALSE)
  NorEpCAMInterDF <- CalKNNBySpatstat(sce, ROIs = ROIs, CenType = "EpCAM_cell", NeighType = "SubType", k = 10, xlim = c(0, 1000), ylim = c(0, 1000), return.cellID = FALSE)

  save.image("./temp.RData")

  ## Take mean
  CA9EpCAMInterDFMean <- CA9EpCAMInterDF %>%
    group_by(PID) %>%
    summarise(across(c(1:(ncol(CA9EpCAMInterDF) - 1)), mean, na.rm = TRUE))
  CA9EpCAMInterDFMean$Group <- rep(1, nrow(CA9EpCAMInterDFMean))

  NorEpCAMInterDFMean <- NorEpCAMInterDF %>%
    group_by(PID) %>%
    summarise(across(c(1:(ncol(CA9EpCAMInterDF) - 1)), mean, na.rm = TRUE))
  NorEpCAMInterDFMean$Group <- rep(0, nrow(NorEpCAMInterDFMean))

  MeanDF <- rbind(CA9EpCAMInterDFMean, NorEpCAMInterDFMean)

  mat <- MeanDF[, c(2:19, 23)]
  mat <- apply(mat, MARGIN = 2, function(x) {
    return(as.numeric(x))
  })
  mat_foldchangeMat <- FCandPvalueCal(mat = mat, xCol = c(1, 18), yCol = 19, need.sample = F)
  mat_foldchangeMat$Q.value <- p.adjust(mat_foldchangeMat$P.value, method = "BH")

  VolcanoPlot(mat_foldchangeMat, pthreshold = 0.05, fcthreshold = 1.2, feature = NULL, filename = paste0("./Neighbors difference between CAIX+ and Other EpCAM+ cell.pdf"))
}
