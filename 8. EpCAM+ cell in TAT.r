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
source("./GA_functions.r")

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
sce <- sce[, sce$Tissue == "TAT"]

clinical <- read.table("/mnt/data/lyx/IMC/clinical.csv", sep = ",", header = T)

savePath <- "/mnt/data/lyx/IMC/analysis/TAT/TC/"

## The k nearest neighbours of Cells in TAT
Minor_cell_neighbors_knn10 <- read.csv("/mnt/data/lyx/IMC/analysis/spatial/SubType_spatial_count_knn_10_TAT.csv")
Minor_cell_neighbors_knn10 <- Minor_cell_neighbors_knn10[, -1]
head(Minor_cell_neighbors_knn10)

Major_cell_neighbors_knn10 <- read.csv("/mnt/data/lyx/IMC/analysis/spatial/MajorType_spatial_count_knn_10_TAT.csv")
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

  DF <- MeanDensity

  y <- DF$RFS_status
  x <- as.matrix(DF[, 2:(ncol(DF) - 2)])
  x <- apply(x, MARGIN = 2, function(a) {
    return(scale(a, center = T, scale = T))
  })

  ## select marker
  features <- c("B", "SC_COLLAGEN", "CD4T", "Macro_CD163", "SC_FAP")
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
  plot(roc_obj, main = paste0(savePath, "Logistic model of neighbors to TC (n=", length(y), ")"), col = "blue", lwd = 4, print.auc = TRUE, auc.polygon = TRUE, grid = TRUE, grid.col = "darkgray", grid.lwd = 1.5, auc.polygon.col = "skyblue")
  dev.off()

  ## KM
  DF$RiskScore <- as.numeric(prob_predictions)
  if (T) {
    write.table(DF[, !(colnames(DF) %in% c("RFS_time", "RFS_status"))], "EpCAM cell neighbors density in TAT for model construction.csv", sep = ",", row.names = T, col.names = T)
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

## Heatmap for tow classes bile duct environment
if(T){
  library(pheatmap)

df <- InterDF[, match(c("PID", "RFS_time", "RFS_status", "B", "CD4T", "SC_FAP", "Macro_CD163", "SC_COLLAGEN"), colnames(InterDF))]

y <- df[, 1:3]

x <- df[, 4:ncol(df)]
# x <- apply(x,MARGIN = 1, function(a){return(scale(as.numeric(a)))})
# x <- t(x)

d <- dist(x, method = "manhattan")
cl <- hclust(d, method = 'ward.D')
groups <- cutree(cl, k = 5)
y$Label <- unname(groups)
table(y$RFS_status, y$Label)

RowAnno <- y[, c(3:4)]
RowAnno[,1] <- as.factor(RowAnno[,1])
RowAnno[,2] <- as.factor(RowAnno[,2])

plotdf <- t(x)

plotdf <- plotdf[,match(rownames(RowAnno[order(RowAnno$Label),]),colnames(plotdf))]

p <- pheatmap(plotdf,
  cluster_rows = F, cluster_cols = F, scale = "column",
  clustering_distance_cols = "manhattan",clustering_method = "ward.D",
  annotation_col = RowAnno
)

pdf("test.pdf", width = 15, height = 2)
print(p)
dev.off()
}

