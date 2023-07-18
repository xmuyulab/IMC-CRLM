## Analysis the difference between Relapse and non-relapse

library(SingleCellExperiment)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggforce)

library(survival)
library(survminer)
library(pROC)
library(dbscan)

library(dplyr)
library(tidyr)

source("./structural_analysis_functions.r")

## load TAT sce object
sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
sce <- sce[, sce$Tissue %in% c("TAT", "IM", "CT")]

savePath <- "/mnt/data/lyx/IMC/analysis/Temp/"


metaMarkers <- c(
  "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
  "VEGF", "CAIX" ## Hypoxia
)

immuneMarkers <- c(
  "CD274", "CD279", "TIGIT", "CD366", ## Immuno-checkpoint
  "CD127", "CD27", "CD80" ## Immuno-activation
)

## Cell subpopulation expression difference in tissues
if (T) {
  sce <- sce[, sce$MajorType != "UNKNOWN"]

  df <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 7))

  for (tissue in c("CT", "IM", "TAT")) {
    sce_ <- sce[, sce$Tissue == tissue]

    ## all type
    sce2 <- sce_
    sce2$SubType <- "All"
    allDF <- HeatmapForMarkersInGroups(sce = sce2, markers = metaMarkers, groupCol = "RFS_status", adjust.p = T, FCthreshold = 0, savePath = NULL, return.mat = T, sample.size = 5e4)

    ## Major Type
    sce2 <- sce_
    sce2$SubType <- sce2$MajorType
    MajorDF <- HeatmapForMarkersInGroups(sce = sce2, markers = metaMarkers, groupCol = "RFS_status", adjust.p = T, FCthreshold = 0, savePath = NULL, return.mat = T, sample.size = 1.5e4)

    ## Lymphocytes
    # Intertypes <- c("B", "CD4T", "CD8T", "NK", "Treg","Macro_CD163","Macro_CD169","Mono_CD11c")

    # sce2 <- sce_
    # sce2 <- sce2[, sce2$SubType %in% Intertypes]

    # LymphoDF <- HeatmapForMarkersInGroups(sce = sce2, markers = c(metaMarkers,immuneMarkers), groupCol = "RFS_status", adjust.p = T, FCthreshold = 0, savePath = NULL, return.mat = T, sample.size = 1e4)

    ## Merge
    dfTemp <- rbind(allDF, MajorDF)
    # dfTemp <- rbind(dfTemp, LymphoDF)
    dfTemp$Tissue <- rep(tissue, nrow(dfTemp))

    ## Save
    df <- rbind(df, dfTemp)
  }


  ## Visualization
  if (T) {
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

    FC_threshold <- 1.5

    allSubTypes <- names(table(sce$SubType))

    df$MajorType <- ifelse(df$Subtype %in% allSubTypes,
      colData(sce)[match(df$Subtype, sce$SubType), "MajorType"],
      "MajorType"
    )

    plotdf <- df[!(df$MajorType %in% c("MajorType")), ]

    # Lollipop chart
    p <- ggplot(plotdf, aes(x = reorder(Subtype, log2FoldChange), y = log2FoldChange, fill = Subtype)) +
      geom_segment(aes(x = Subtype, xend = Subtype, y = 0, yend = log2FoldChange), color = "grey", linetype = "dashed") +
      geom_point(shape = 21, size = 4, aes(alpha = Qlabel)) +
      scale_alpha_continuous(range = c(0.2, 1)) + # Modify transparency range to make points more visible
      scale_fill_manual(values = pal_ucscgb("default")(n = length(unique(plotdf$Marker)))) + # Custom color palette
      geom_hline(yintercept = log2(FC_threshold), linetype = "dashed", color = "black") + # Add dashed line
      geom_hline(yintercept = log2(1 / FC_threshold), linetype = "dashed", color = "black") +
      theme_light(base_size = 14) +
      theme(
        legend.position = "right",
        axis.text.x = element_blank(),
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

    pdf(paste0(savePath, "Metabolic Difference of Major types between relapse.pdf"), height = 6, width = 15)
    print(p)
    dev.off()
  }
}

## EpCAM+ cells structure detect
library(topicmodels)
library(deldir)
library(spatstat)

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
sce <- sce[, sce$Tissue == "TAT"]

ROIs <- names(table(sce$ID))

colData(sce)[sce$MajorType %in% "Tumor", "SubType"] <- "EpCAM+ cell"

## LDA to clustering the EpCAM+ cells
if (F) {
  # Get the neighbors within D
  k_NeighborsList <- CalDNNBySpatstat(sce, ROIs = ROIs, CenType = "EpCAM+ cell", NeighType = "SubType", d = 20, xlim = c(0, 1000), ylim = c(0, 1000))

  ## Combine all neighbors
  neighbors <- bind_rows(k_NeighborsList)

  ## omit all zero cells
  neighbors$PID <- sapply(neighbors$ID, function(x) {
    strsplit(x, split = "_")[[1]][1]
  })
  neighbors <- left_join(neighbors, clinical[, match(c("PID", "RFS_status"), colnames(clinical))], by = "PID")

  ## split
  intersType <- c("B", "CD4T", "SC_FAP", "SC_COLLAGEN", "Macro_CD163")
  neighbors_ <- neighbors[, match(c(intersType, "RFS_status"), colnames(neighbors))]

  neighbors_ <- neighbors_[apply(neighbors_[, 1:(ncol(neighbors_) - 1)], 1, function(x) {
    sum(x) != 0
  }), ]

  metaCol <- neighbors_[, ncol(neighbors_)]
  neighborsCol <- neighbors_[, 1:length(intersType)]

  ## LDA
  lda_model <- LDA(neighborsCol, k = 2, control = list(seed = 619))

  as.matrix(terms(lda_model, 6))

  gamma <- lda_model@gamma
  gamma_df <- as.data.frame(gamma)
  doc_labels <- apply(gamma_df, 1, which.max)

  test <- cbind(metaCol, doc_labels)
  colnames(test) <- c("RFS_status", "labels")

  tab <- table(test[, 1], test[, 2])
}

centralTypes <- c("EpCAM+ cell")
NeighborTypes <- c("SC_COLLAGEN")

sce$cluster_size <- NA
sce$InterestypeDensity <- NA

for (ROI in ROIs) {
  savePathTemp <- paste0(savePath, ROI, "/")
  if (!dir.exists(savePathTemp)) {
    dir.create(savePathTemp, recursive = T)
  }

  ## Extract the figure and certain cell subpopulation
  test <- sce[, sce$ID == ROI]

  ## Delaunay triangulation
  filtered_data <- test[, test$SubType %in% c(centralTypes, NeighborTypes)]

  filtered_coor <- GetXandY(filtered_data$Position)

  # Create point pattern
  filtered_ppp <- ppp(filtered_coor[, 1], filtered_coor[, 2], xrange = c(0, 1000), yrange = c(0, 1000))
  filtered_coords <- as.data.frame(filtered_ppp)
  filtered_coords$SubType <- filtered_data$SubType
  filtered_coords$CellID <- filtered_data$CellID

  # Perform the Delaunay triangulation
  delaunay <- deldir(filtered_coords)

  # Plot the triangulation
  pdf(paste0(savePathTemp, "Delaunay triangulation graph.pdf"), height = 6, width = 8)
  plot(delaunay)
  dev.off()

  # Convert the delaunay object to an edge dataframe
  edges <- with(delaunay$delsgs, data.frame(from = ind1, to = ind2, length = sqrt((x1 - x2)^2 + (y1 - y2)^2)))

  # Remove the distance larger than threshold
  d_threshold <- 20
  edges <- edges[edges$length <= d_threshold, ]

  # Filter out vertices that are not in the edge list
  filtered_coords <- filtered_coords[unique(c(edges$from, edges$to)), ]

  # Create a graph object using igraph package
  g <- graph_from_data_frame(edges, directed = FALSE)

  ## igraph clustering
  if (F) {
    # Perform Walktrap community detection
    communities <- walktrap.community(g)

    # Get membership vector (cluster assignments)
    membership <- communities$membership
  }

  ## Louvain clustering
  if (T) {
    # Compute the Louvain communities
    communities <- cluster_louvain(g, resolution = 0.2)

    # Get membership vector (cluster assignments)
    membership <- communities$membership
  }

  ## DBSCAN clustering
  if (F) {
    # Run DBSCAN on the distance matrix
    # eps and MinPts need to be chosen based on your data
    dist_matrix <- get.adjacency(g, sparse = FALSE)
    dbscan_result <- dbscan::dbscan(dist_matrix, eps = 20, minPts = 2)

    # Get membership vector (cluster assignments), note that noise points are assigned to cluster 0
    membership <- dbscan_result$cluster
  }

  # Calculate cluster sizes
  cluster_sizes <- table(membership)

  # Remove small clusters (less than 10 members)
  # small_clusters <- names(cluster_sizes[cluster_sizes < 7])
  # membership[membership %in% small_clusters] <- NA

  # Add cluster results back to the data
  filtered_coords$cluster <- membership

  # Identify clusters without the central type
  clusters_with_central_type <- unique(filtered_coords$cluster[filtered_coords$SubType %in% centralTypes])
  clusters_without_central_type <- setdiff(unique(filtered_coords$cluster), clusters_with_central_type)
  filtered_coords$cluster[filtered_coords$cluster %in% clusters_without_central_type] <- NA


  # Calculate type count within each cluster
  type_count <- table(filtered_coords$cluster, filtered_coords$SubType)

  # Calculate type density (proportion of each type in each cluster)
  type_density <- prop.table(type_count, 1) # normalize rows
  type_density <- data.frame(InterestypeDensity = 1 - type_density[, 1])
  type_density$cluster <- rownames(type_density)

  # Print results
  #print(cluster_sizes)
  #print(type_density)

  filtered_coords$cluster <- as.character(filtered_coords$cluster)
  filtered_coords <- as.data.frame(left_join(filtered_coords, type_density, by = "cluster"))
  filtered_coords$cluster_size <- cluster_sizes[as.numeric(filtered_coords$cluster)]

  # Add cluster column and cluster size to test data frame by matching the CellID
  meta <- as.data.frame(colData(test))
  all_coor <- GetXandY(test$Position)

  meta$x_coordinate <- all_coor[, 1]
  meta$y_coordinate <- all_coor[, 2]

  meta$cluster <- 0
  meta$cluster_size <- 0
  meta$InterestypeDensity <- 0

  meta$cluster[match(filtered_coords$CellID, meta$CellID)] <- filtered_coords$cluster
  meta$cluster_size[match(filtered_coords$CellID, meta$CellID)] <- filtered_coords$cluster_size
  meta$InterestypeDensity[match(filtered_coords$CellID, meta$CellID)] <- filtered_coords$InterestypeDensity
  meta <- replace(meta, is.na(meta), 0)

  # Scatterplot of cells colored by cluster
  colors <- c("grey", rep(rainbow(100), each = 1))

  p <- ggplot(meta, aes(x = x_coordinate, y = y_coordinate, color = as.factor(cluster))) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = colors) +
    theme_minimal() +
    labs(color = "Cluster") +
    ggtitle("Cells colored by cluster")

  pdf(paste0(savePathTemp, "cluster.pdf"), height = 6, width = 8)
  print(p)
  dev.off()

  # Scatterplot of cells colored by cluster size
  p <- ggplot(meta, aes(x = x_coordinate, y = y_coordinate, color = cluster_size)) +
    geom_point(alpha = 0.5) +
    scale_color_gradientn(colors = c("grey", "blue", "yellow"), values = c(0, 0.0001, 1)) +
    theme_minimal() +
    labs(color = "Cluster size") +
    ggtitle("Cells colored by cluster size")

  pdf(paste0(savePathTemp, "cluster size.pdf"), height = 6, width = 8)
  print(p)
  dev.off()

  # Scatterplot of cells colored by type density
  p <- ggplot(meta, aes(x = x_coordinate, y = y_coordinate, color = InterestypeDensity)) +
    geom_point(alpha = 0.5) +
    scale_color_gradientn(colors = c("grey", "blue", "yellow"), values = c(0, 0.0001, 1)) +
    theme_minimal() +
    labs(color = "Type density") +
    ggtitle("Cells colored by type density")

  pdf(paste0(savePathTemp, "Interestype Density.pdf"), height = 6, width = 8)
  print(p)
  dev.off()

  ## match the cluster size and density to original sce
  sce$cluster_size[match(meta$CellID, sce$CellID)] <- meta$cluster_size
  sce$InterestypeDensity[match(meta$CellID,sce$CellID)] <- meta$InterestypeDensity
}
saveRDS(sce,paste0(savePath,"sce.rds"))

## Investigate the marker expression and cluster size
sce <- sce[, sce$SubType %in% centralTypes]
log2clustersize <- log2(sce$cluster_size + 1)

## Get certain marker expression
marker <- "CAIX"
value <- assay(sce[marker, ])

# Define breaks
breaks <- seq(0, 8, by = 0.5)
# Cut the vector into intervals
vec_cut <- cut(vec, breaks = breaks)
vec_group_labels <- as.character(vec_cut)

TC_Marker <- cbind(t(value), vec_group_labels)
TC_Marker <- as.data.frame(TC_Marker)
colnames(TC_Marker) <- c("CAIX","log2ClusterSize")

## downsample
set.seed(1)
TC_Marker <- TC_Marker[sample.int(nrow(TC_Marker),size = 1000),]
TC_Marker[,1] <- as.numeric(TC_Marker[,1])

p <- ggplot(TC_Marker, aes(x = log2ClusterSize, y = CAIX)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Log2 Cluster Size") +
  ylab("CAIX Expression Value")

png("test.png")
print(p)
dev.off()

## clustering size with survival ?
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
