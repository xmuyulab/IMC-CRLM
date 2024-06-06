## Analysis the difference between Relapse and non-relapse

library(SingleCellExperiment)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggforce)
library(ggsignif)

library(survival)
library(survminer)
library(pROC)
library(dbscan)

library(dplyr)
library(tidyr)

source("./structural_analysis_functions.r")
source("./GA_functions.r")

## load IM sce object
### Distance to Tumor
sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")

savePath <- "/mnt/data/lyx/IMC/analysis/GradientAna/"
if (!dir.exists(savePath)) {
  dir.create(savePath, recursive = T)
}

metaMarkers <- c(
  "PRPS1", "FASN", "GLUT1", "HK2", "Ki67", ## Cell Growth and Division
  "VEGF", "CAIX" ## Hypoxia
)

immuneMarkers <- c(
  "CD274", "CD279", "TIGIT", "CD366", ## Immuno-checkpoint
  "CD127", "CD27", "CD80" ## Immuno-activation
)

sce_ <- sce[, sce$Tissue %in% c("IM")]

targetType <- "Tumor"
targetTypeCol <- "MajorType"

sce_ <- Dis2Boundary(sce_, targetType = targetType, targetTypeCol = targetTypeCol, imageIDcol = "ID", coorCol = "Position", DisColName = "Dis2Tumor", k = 10)
saveRDS(sce_, paste0(savePath, "sce_IM_with_distance2Tumor.rds"))

## Gradient analysis
sce_ <- readRDS(paste0(savePath, "sce_IM_with_distance2Tumor.rds"))

head(colData(sce_))

df <- as.data.frame(colData(sce_))

## Seperate distance to tumor analysis
if (T) {
  ## Define distance breaks
  maxDis <- as.integer(max(df$Dis2Tumor) + 1)
  breaks <- c(0, 50, 100, 150, 200, 300, 500, maxDis)

  df_ranges <- df
  df_ranges <- df_ranges %>%
    mutate(distance_range = cut(Dis2Tumor, breaks = breaks)) %>%
    # First, compute the counts for each combination of ID, distance_range, SubType, and RFS_status
    group_by(ID, distance_range, SubType, RFS_status) %>%
    summarise(n = n(), .groups = "drop") %>%
    # Now, group by ID and SubType to get the total number of cells for each cell type
    group_by(ID, SubType) %>%
    mutate(total_cells_cell_type = sum(n), cell_fraction = n / total_cells_cell_type)


  # plotType <- c("B", "CD4T", "CD8T", "NK", "Treg", "Macro_CD163", "Macro_CD169", "Mono_CD11c", "Macro_HLADR")
  plotType <- names(table(df$SubType))
  plotType <- plotType[c(1:15,20)]


  savePathTemp <- paste0(savePath, "Type density to tumor (seperate)/")
  if (!dir.exists(savePathTemp)) {
    dir.create(savePathTemp, recursive = T)
  }

  for (type_ in plotType) {
    df_type <- df_ranges[df_ranges$SubType %in% type_, ]
    df_type <- as.data.frame(df_type)
    df_type <- df_type[, -c(1, 3)]

    plotdf <- df_type %>%
      group_by(RFS_status, distance_range) %>%
      summarise(across(c(1:(ncol(df_type) - 2)), mean, na.rm = TRUE))

    plotdf$RFS_status <- as.factor(plotdf$RFS_status)

    p <- ggplot(plotdf, aes(x = distance_range, y = cell_fraction, group = RFS_status, color = RFS_status)) +
      geom_line(aes(linetype = RFS_status), size = 1) +
      geom_point(size = 3) +
      scale_color_manual(values = ggsci::pal_jco("default")(2)) +
      labs(
        title = "Cell Fraction Change Along Distance Range",
        x = "Distance Range",
        y = "Cell Fraction",
        color = "RFS_status",
        linetype = "RFS_status"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      stat_compare_means(aes(group = RFS_status),
        data = df_type,
        method = "t.test",
        hide.ns = FALSE,
        label = "p.format"
      )

    pdf(paste0(savePathTemp, type_, " Density to Tumor Boundary.pdf"), width = 8, height = 6)
    print(p)
    dev.off()
  }
}

## Continuous distance to tumor analysis
if (T) {
  ## Define distance breaks
  maxDis <- as.integer(500)

  df_ranges <- df[df$Dis2Tumor <= maxDis, ]


  plotType <- c("B", "CD4T", "CD8T", "NK", "Treg", "Macro_CD163", "Macro_CD169", "Mono_CD11c", "Macro_HLADR")


  savePathTemp <- paste0(savePath, "Type density to tumor (continuous)/")
  if (!dir.exists(savePathTemp)) {
    dir.create(savePathTemp, recursive = T)
  }

  for (type_ in plotType) {
    df_type <- df_ranges[df_ranges$SubType %in% type_, ]
    df_type <- as.data.frame(df_type)
    df_type <- df_type[, -c(1, 3)]

    df_stat <- df[df$SubType %in% type_, ]
    df_stat <- df_stat[, c("ID", "RFS_status", "Dis2Tumor")]

    ## Statistical Testing
    df_stat <- df_stat %>%
      group_by(ID) %>%
      summarise(across(c(1:(ncol(df_stat) - 1)), mean, na.rm = TRUE))

    wilcox_test_result <- wilcox.test(Dis2Tumor ~ RFS_status, data = df_stat)
    p_value <- wilcox_test_result$p.value

    ## Plot figure
    df_type$RFS_status <- as.factor(df_type$RFS_status)

    p <- ggplot(df_type, aes(x = Dis2Tumor, fill = RFS_status, color = RFS_status)) +
      geom_line(stat = "density", size = 1) +
      scale_color_brewer(palette = "Set1") +
      labs(subtitle = type_, x = "Distance to Tumor Boundary", y = "Density") +
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

    pdf(paste0(savePathTemp, type_, " Density to Tumor Boundary.pdf"), width = 8, height = 6)
    print(p)
    dev.off()
  }
}


### Distance to EpCAM+ cell in TAT
sce_ <- sce[, sce$Tissue %in% c("TAT")]
sce_$MajorType[sce_$MajorType == "Tumor"] <- "EpCAM+ cell"

targetType <- "EpCAM+ cell"
targetTypeCol <- "MajorType"

sce_ <- Dis2Boundary(sce_, targetType = targetType, targetTypeCol = targetTypeCol, imageIDcol = "ID", coorCol = "Position", DisColName = "Dis2BD", k = 10)
saveRDS(sce_, paste0(savePath, "sce_TAT_with_distance2BileDuct.rds"))

## Gradient analysis
sce_ <- readRDS(paste0(savePath, "sce_TAT_with_distance2BileDuct.rds"))

head(colData(sce_))

df <- as.data.frame(colData(sce_))

## Do not seperate into breaks
plotType <- c("B", "CD4T", "CD8T", "NK", "Treg", "Macro_CD163", "Macro_CD169", "Mono_CD11c", "Macro_HLADR", "SC_COLLAGEN", "SC_Vimentin", "SC_FAP")

## Seperate distance to bile duct analysis
if (T) {
  ## Define distance breaks
  breaks <- c(0, 50, 100, 200, 300, 500, 1000)

  ## Calculate the cell subpopulation fraction across the distance
  maxDis <- max(breaks)

  df_ranges <- df
  df_ranges <- df_ranges %>%
    mutate(distance_range = cut(Dis2BD, breaks = breaks)) %>%
    # First, compute the counts for each combination of ID, distance_range, SubType, and RFS_status
    group_by(ID, distance_range, SubType, RFS_status) %>%
    summarise(n = n(), .groups = "drop") %>%
    # Now, group by ID and SubType to get the total number of cells for each cell type
    group_by(ID, SubType) %>%
    mutate(total_cells_cell_type = sum(n), cell_fraction = n / total_cells_cell_type)

  savePathTemp <- paste0(savePath, "Type density to bile duct (seperate)/")
  if (!dir.exists(savePathTemp)) {
    dir.create(savePathTemp, recursive = T)
  }

  for (type_ in plotType) {
    df_type <- df_ranges[df_ranges$SubType %in% type_, ]
    df_type <- as.data.frame(df_type)
    df_type <- df_type[, -c(1, 3)]

    plotdf <- df_type %>%
      group_by(RFS_status, distance_range) %>%
      summarise(across(c(1:(ncol(df_type) - 2)), mean, na.rm = TRUE))

    plotdf$RFS_status <- as.factor(plotdf$RFS_status)

    p <- ggplot(plotdf, aes(x = distance_range, y = cell_fraction, group = RFS_status, color = RFS_status)) +
      geom_line(aes(linetype = RFS_status), size = 1) +
      geom_point(size = 3) +
      scale_color_manual(values = ggsci::pal_jco("default")(2)) +
      labs(
        title = "Cell Fraction Change Along Distance Range",
        x = "Distance Range",
        y = "Cell Fraction",
        color = "RFS_status",
        linetype = "RFS_status"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      stat_compare_means(aes(group = RFS_status),
        data = df_type,
        method = "t.test",
        hide.ns = FALSE,
        label = "p.format"
      )

    pdf(paste0(savePathTemp, type_, " Density to Bile Duct.pdf"), width = 8, height = 6)
    print(p)
    dev.off()
  }
}

## Continuous distance to bile duct analysis
if (T) {
  ## Define distance breaks
  maxDis <- as.integer(500)

  df_ranges <- df[df$Dis2BD <= maxDis, ]

  savePathTemp <- paste0(savePath, "Type density to bile duct (continuous)/")
  if (!dir.exists(savePathTemp)) {
    dir.create(savePathTemp, recursive = T)
  }

  for (type_ in plotType) {
    df_type <- df_ranges[df_ranges$SubType %in% type_, ]
    df_type <- as.data.frame(df_type)
    df_type <- df_type[, -c(1, 3)]

    df_stat <- df[df$SubType %in% type_, ]
    df_stat <- df_stat[, c("ID", "RFS_status", "Dis2BD")]

    ## Statistical Testing
    df_stat <- df_stat %>%
      group_by(ID) %>%
      summarise(across(c(1:(ncol(df_stat) - 1)), mean, na.rm = TRUE))

    wilcox_test_result <- wilcox.test(Dis2BD ~ RFS_status, data = df_stat)
    p_value <- wilcox_test_result$p.value

    ## Plot figure
    df_type$RFS_status <- as.factor(df_type$RFS_status)

    p <- ggplot(df_type, aes(x = Dis2BD, fill = RFS_status, color = RFS_status)) +
      geom_line(stat = "density", size = 1) +
      scale_color_brewer(palette = "Set1") +
      labs(subtitle = type_, x = "Distance to Bile Duct", y = "Density") +
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

    pdf(paste0(savePathTemp, type_, " Density to Bile Duct.pdf"), width = 8, height = 6)
    print(p)
    dev.off()
  }
}
