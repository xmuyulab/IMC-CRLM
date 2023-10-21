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
library(deldir)
library(spatstat)

library(dplyr)
library(tidyr)

source("./spatial_analysis_functions.r")
source("./structural_analysis_functions.r")

## load TAT sce object
sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
sce <- sce[, sce$Tissue %in% c("TAT", "IM", "CT")]

savePath <- "/mnt/data/lyx/IMC/analysis/NewStructure/"


metaMarkers <- c(
  "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
  "VEGF", "CAIX" ## Hypoxia
)

immuneMarkers <- c(
  "CD274", "CD279", "TIGIT", "CD366", ## Immuno-checkpoint
  "CD127", "CD27", "CD80" ## Immuno-activation
)

## Cell subpopulation expression difference in tissues
if (F) {
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

    # plotdf <- df[!(df$MajorType %in% c("MajorType")), ]
    plotdf <- df

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
sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
sce <- sce[, sce$Tissue == "TAT"]

savePath <- "/mnt/data/lyx/IMC/analysis/NewStructure/"

ROIs <- names(table(sce$ID))

colData(sce)[sce$MajorType %in% "Tumor", "SubType"] <- "EpCAM+ cell"

centralTypes <- c("EpCAM+ cell")


for (type in c("Inhibit", "Immune")) {
  sceTemp <- sce

  sceTemp$cluster_size <- NA
  sceTemp$InterestypeDensity <- NA

  savePath_ <- paste0(savePath, type, "/")

  if (type == "Inhibit") {
    NeighborTypes <- c("SC_COLLAGEN", "Macro_CD163")
  }
  if (type == "Immune") {
    NeighborTypes <- c("B", "CD4T", "SC_FAP")
  }

  for (ROI in ROIs) {
    savePathTemp <- paste0(savePath_, ROI, "/")
    if (!dir.exists(savePathTemp)) {
      dir.create(savePathTemp, recursive = T)
    }

    ## Extract the figure and certain cell subpopulation
    test <- sceTemp[, sceTemp$ID == ROI]

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
    if (F) {
      # Compute the Louvain communities
      communities <- cluster_louvain(g, resolution = 0.2)

      # Get membership vector (cluster assignments)
      membership <- communities$membership
    }

    ## Just extreact connected subgraph
    if (T) {
      # Identify the connected components
      comp <- components(g)

      # The membership element of the result gives the characteristic vector
      membership <- comp$membership
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

    ## match the cluster size and density to original sceTemp
    sceTemp$cluster_size[match(meta$CellID, sceTemp$CellID)] <- meta$cluster_size
    sceTemp$InterestypeDensity[match(meta$CellID, sceTemp$CellID)] <- meta$InterestypeDensity

    ## Visualize
    VisualizeMarkerAndType(sce = sceTemp, ROI = ROI, celltype = c(centralTypes, NeighborTypes), whichType = "SubType", marker = NULL, SavePath = savePathTemp)
  }
  saveRDS(sceTemp, paste0(savePath, "TAT_sce_with_", type, "_stucture.rds"))
}

# Inhibit
if (F) {
  sce <- readRDS(paste0(savePath, "TAT_sce_with_Inhibit_stucture.rds"))

  ## Investigate the marker expression and cluster size
  sce_ <- sce[, sce$SubType %in% centralTypes]
  savePath_ <- paste0(savePath, "Inhibit/")

  # sce_ <- sce[,sce$cluster_size != 0]
  log2clustersize <- log2(sce_$cluster_size + 1)

  # Define breaks
  breaks <- c(0, 2, 3, 4, 5, 6, 8, ceiling(max(log2clustersize)))
  # Cut the vector into intervals
  vec_cut <- cut(log2clustersize, breaks = breaks, include.lowest = T, right = F)
  vec_group_labels <- as.character(vec_cut)

  markers <- metaMarkers

  ## Get certain marker expression
  for (marker in markers) {
    value <- assay(sce_[marker, ])
    # value <- rescale(value,new_range = c(0,1))

    TC_Marker <- cbind(t(value), vec_group_labels)
    TC_Marker <- as.data.frame(TC_Marker)
    colnames(TC_Marker) <- c("marker", "log2ClusterSize")

    ## downsample
    set.seed(1)
    TC_Marker <- TC_Marker[sample.int(nrow(TC_Marker), size = 5000), ]
    TC_Marker[, 1] <- as.numeric(TC_Marker[, 1])

    p <- ggplot(TC_Marker, aes(x = log2ClusterSize, y = marker)) +
      geom_boxplot(fill = "darkgreen") +
      theme_bw() +
      xlab("Log2 Cluster Size") +
      ylab(paste0(marker, " Expression Value")) +
      stat_compare_means(aes(group = log2ClusterSize), method = "anova", label.x.npc = "middle")


    pdf(paste0(savePath_, marker, " change and inhibit size.pdf"), width = 5, height = 3)
    print(p)
    dev.off()
  }
}

# Immune
if (F) {
  sce <- readRDS(paste0(savePath, "TAT_sce_with_Immune_stucture.rds"))
  savePath_ <- paste0(savePath, "Immune/")

  ## Investigate the marker expression and cluster size
  sce_ <- sce[, sce$SubType %in% centralTypes]
  # sce_ <- sce[,sce$cluster_size != 0]
  log2clustersize <- log2(sce_$cluster_size + 1)

  # Define breaks
  breaks <- c(0, 2, 3, 4, 5, 6, 8, ceiling(max(log2clustersize)))
  # Cut the vector into intervals
  vec_cut <- cut(log2clustersize, breaks = breaks, include.lowest = T, right = F)
  vec_group_labels <- as.character(vec_cut)

  markers <- metaMarkers

  ## Get certain marker expression
  for (marker in markers) {
    value <- assay(sce_[marker, ])

    TC_Marker <- cbind(t(value), vec_group_labels)
    TC_Marker <- as.data.frame(TC_Marker)
    colnames(TC_Marker) <- c("marker", "log2ClusterSize")

    ## downsample
    set.seed(1)
    TC_Marker <- TC_Marker[sample.int(nrow(TC_Marker), size = 5000), ]
    TC_Marker[, 1] <- as.numeric(TC_Marker[, 1])

    p <- ggplot(TC_Marker, aes(x = log2ClusterSize, y = marker)) +
      geom_boxplot(fill = "darkgreen") +
      theme_bw() +
      xlab("Log2 Cluster Size") +
      ylab(paste0(marker, " Expression Value")) +
      stat_compare_means(aes(group = log2ClusterSize), method = "anova", label.x.npc = "middle")


    pdf(paste0(savePath_, marker, " change and immune size.pdf"), width = 5, height = 3)
    print(p)
    dev.off()
  }
}

## Compare inhibit and immune
if (F) {
  inhibitSCE <- readRDS("/mnt/data/lyx/IMC/analysis/NewStructure/TAT_sce_with_Inhibit_stucture.rds")
  immuneSCE <- readRDS("/mnt/data/lyx/IMC/analysis/NewStructure/TAT_sce_with_Immune_stucture.rds")

  inhibitSCE <- inhibitSCE[, inhibitSCE$SubType %in% centralTypes]
  immuneSCE <- immuneSCE[, immuneSCE$SubType %in% centralTypes]

  sce_ <- inhibitSCE
  colnames(colData(sce_))[49:50] <- c("Inhibit_cluster_size", "Inhibit_InterestypeDensity")
  colData(sce_)$Immune_cluster_size <- immuneSCE$cluster_size
  colData(sce_)$Immune_InterestypeDensity <- immuneSCE$InterestypeDensity

  head(colData(sce_))

  test <- as.data.frame(colData(sce_)[, c("PID", "RFS_time", "RFS_status", "Inhibit_cluster_size", "Inhibit_InterestypeDensity", "Immune_cluster_size", "Immune_InterestypeDensity")])
  head(test)

  test2 <- test %>%
    group_by(PID) %>%
    summarise(across(c(1:(ncol(test) - 1)), mean, na.rm = TRUE))

  df <- as.data.frame(test2)
  df$TestLabel <- df$Inhibit_cluster_size / df$Immune_cluster_size
  df <- df[,c(1,2,3,8)]
  df <- df[order(df$TestLabel,decreasing = T),]
  
  cutpoint <- surv_cutpoint(data = df, time = "RFS_time", event = "RFS_status", variables = "TestLabel")
  cutpoint <- summary(cutpoint)$cutpoint

  df$Label <- ifelse(df$TestLabel > cutpoint, 1, 0)

  df$RFS_status <- as.numeric(df$RFS_status)
  df$RFS_time <- as.numeric(df$RFS_time)

  ## km curve
  fit <- survfit(Surv(RFS_time, RFS_status) ~ Label, data = df)
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

  pdf(paste0(savePath_, "KM for Label (Patient).pdf"), width = 8, height = 6)
  print(p)
  dev.off()

  test2$Diff <- test2$Inhibit_cluster_size - test2$Immune_cluster_size
  test2 <- test2[order(test2$Diff), ]
}

  ## Using the mean size of cluster to survival analysis
if (T) {
  DF <- as.data.frame(colData(sce_)[, match(c("PID", "RFS_time", "RFS_status"), colnames(colData(sce_)))])
  DF$log2clustersize <- log2clustersize

  MeanDF <- DF %>%
    group_by(PID) %>%
    summarise(across(c(1:(ncol(DF) - 1)), mean, na.rm = TRUE))

  MeanDF <- as.data.frame(MeanDF[order(MeanDF$log2clustersize,decreasing = T),])

  ## KM
  if (F) {
    write.table(DF, "EpCAM cell neighbors density in TAT for model construction.csv", sep = ",", row.names = T, col.names = T)
  }

  cutpoint <- surv_cutpoint(data = MeanDF, time = "RFS_time", event = "RFS_status", variables = "log2clustersize")
  cutpoint <- summary(cutpoint)$cutpoint

  df <- MeanDF[, c("log2clustersize", "RFS_status", "RFS_time")]
  df[, 1] <- ifelse(df[, 1] > cutpoint, "high", "low")

  df$RFS_status <- as.numeric(df$RFS_status)
  df$RFS_time <- as.numeric(df$RFS_time)

  ## km curve
  fit <- survfit(Surv(RFS_time, RFS_status) ~ log2clustersize, data = df)
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

  pdf(paste0(savePath_, "KM for log2clustersize (Patient).pdf"), width = 8, height = 6)
  print(p)
  dev.off()
}
