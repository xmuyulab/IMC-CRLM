# clustering functions

## library
library(pheatmap)
library(FlowSOM)
library(Rphenograph)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Rtsne)
library(cowplot)
library(RColorBrewer)
library(SingleCellExperiment)

library(randomForest)
library(caret)

## load data
load_RawData <- function(raw_csv_dir, meta_csv_dir, Area_cutoff = 5) {
  file_name <- list.files(raw_csv_dir, pattern = ".csv$", full = TRUE)

  total.res <- data.frame()
  for (each in file_name) {
    tmp <- read.csv(each, row.names = 1, check.names = FALSE)
    tmp$sample <- gsub(".csv", "", basename(each))

    position <- tmp$Position
    area <- tmp$Area

    removeCol <- c("ObjectNumber", "Position", "Area")
    tmp <- tmp[, !(colnames(tmp) %in% removeCol)]

    tmp <- cbind(tmp, "Position" = position)
    tmp <- cbind(tmp, "Area" = area)

    ## modify the colnames
    for (i in 1:length(colnames(tmp))) {
      if (colnames(tmp)[i] == "SMA") {
        colnames(tmp)[i] <- "AlphaSMA"
      }
    }

    colnames(tmp) <- gsub(pattern = "-", replacement = "", colnames(tmp))

    total.res <- rbind(total.res, tmp)
  }

  meta <- read.csv(meta_csv_dir)

  meta <- (subset(meta, normal == 1 | tumor == 1 | invasive == 1 | TLS == 1))
  meta$ID <- paste0(meta$sample, "_", meta$ROI)

  total.res$ID <- sapply(total.res$sample, function(x) {
    temp <- strsplit(x, split = "_")[[1]]
    return(paste0(temp[1], "_", temp[2]))
  })

  total.res2 <- total.res[total.res$ID %in% meta$ID, ]
  total.res2$Batch <- meta[match(total.res$ID, meta$ID), ]$Batch

  total.res2.qc <- total.res2[total.res2$Area > Area_cutoff, ]

  return(total.res2.qc)
}

## load marker
load_Marker <- function(panel) {
  allMarker <- panel[panel$full == 1, ]$marker
  allMarker <- allMarker[!is.na(allMarker)]

  ## Major markers
  MajorIdMarker <- panel[panel$MajorIden == 1, ]$marker
  MajorIdMarker <- MajorIdMarker[!is.na(MajorIdMarker)]

  ## Minor markers
  ### Lympocytes
  LymphocyteIdMarker <- panel[panel$LymIden == 1, ]$marker
  LymphocyteIdMarker <- LymphocyteIdMarker[!is.na(LymphocyteIdMarker)]
  ### Myeloids
  MyeloidsIdMarker <- panel[panel$MyeIden == 1, ]$marker
  MyeloidsIdMarker <- MyeloidsIdMarker[!is.na(MyeloidsIdMarker)]
  ### Stromal
  StromalIdMarker <- panel[panel$StrIden == 1, ]$marker
  StromalIdMarker <- StromalIdMarker[!is.na(StromalIdMarker)]
  ### Tumor
  TumorIdMarker <- panel[panel$TumIden == 1, ]$marker
  TumorIdMarker <- TumorIdMarker[!is.na(TumorIdMarker)]
  ### Reclstering
  ReclusterMarker <- panel[panel$Reclustering == 1, ]$marker
  ReclusterMarker <- ReclusterMarker[!is.na(ReclusterMarker)]

  MarkerList <- list()
  MarkerList[["All_Marker"]] <- allMarker
  MarkerList[["Major_Marker"]] <- MajorIdMarker
  MarkerList[["Lymphocyte_Marker"]] <- LymphocyteIdMarker
  MarkerList[["Myeloid_Marker"]] <- MyeloidsIdMarker
  MarkerList[["Stromal_Marker"]] <- StromalIdMarker
  MarkerList[["Tumor_Marker"]] <- TumorIdMarker
  MarkerList[["Recluster_Marker"]] <- ReclusterMarker

  return(MarkerList)
}

## normalize data
normData <- function(sce_data, marker_total, censor_val = NULL, arcsinh = FALSE, is.Batch = F, norm_method = "0-1") {
  #' censor_val=0.999'
  # sce_data <- sce_snf_FDR0.1_nozero
  # sce_data <- sce_p95_FDR0.01_median_0712
  # sce_data <- sce_p85_FDR0.01_median_0713
  if (is.Batch) {
    batches <- sce_data$Batch
  }
  if (is.data.frame(sce_data)) {
    exp_data <- sce_data[, c(marker_total)]
  } else {
    exp_data <- data.frame(t(assay(sce_data)), check.names = TRUE)
  }

  if (arcsinh) {
    exp_data <- asinh(exp_data)
  }

  if (is.data.frame(sce_data)) {
    exp_data$filelist <- sce_data$filelist
  } else {
    exp_data$filelist <- sce_data@metadata$filelist
  }
  if (is.Batch) {
    batcheNames <- names(table(batches))
    dat <- as.data.frame(matrix(nrow = 0, ncol = 5))
    for (batch in batcheNames) {
      dattemp <- data.table((exp_data[batches %in% batch, ]) %>% pivot_longer(-filelist, names_to = "channel", values_to = "mc_counts"))
      dattemp$mc_counts <- as.numeric(dattemp$mc_counts)

      if (!is.null(censor_val)) {
        dattemp[, c_counts := censor_dat(mc_counts, censor_val), by = channel]
        if (norm_method == "0-1") {
          dattemp[, c_counts_scaled := ((c_counts - min(c_counts)) / (max(c_counts) - min(c_counts))), by = channel]
          dattemp[c_counts_scaled < 0, c_counts_scaled := 0, by = channel]
        } else if (norm_method == "znorm") {
          dattemp[, c_counts_scaled := ((c_counts - mean(c_counts)) / sd(c_counts)), by = channel]
        } else if (norm_method == "null") {
          dattemp[, c_counts_scaled := c_counts, by = channel]
        }
      } else {
        if (norm_method == "0-1") {
          dattemp[, c_counts_scaled := ((mc_counts - min(mc_counts)) / (max(mc_counts) - min(mc_counts))), by = channel]
          dattemp[c_counts_scaled < 0, c_counts_scaled := 0, by = channel]
        } else if (norm_method == "znorm") {
          dattemp[, c_counts_scaled := ((mc_counts - mean(mc_counts)) / sd(mc_counts)), by = channel]
        }
      }
      dat <- rbind(dat, dattemp)
    }
  } else {
    dat <- data.table(exp_data %>% pivot_longer(-filelist, names_to = "channel", values_to = "mc_counts"))
    dat$mc_counts <- as.numeric(dat$mc_counts)

    if (!is.null(censor_val)) {
      dat[, c_counts := censor_dat(mc_counts, censor_val), by = channel]
      if (norm_method == "0-1") {
        dat[, c_counts_scaled := ((c_counts - min(c_counts)) / (max(c_counts) - min(c_counts))), by = channel]
        dat[c_counts_scaled < 0, c_counts_scaled := 0, by = channel]
      } else if (norm_method == "znorm") {
        dat[, c_counts_scaled := ((c_counts - mean(c_counts)) / sd(c_counts)), by = channel]
      } else if (norm_method == "null") {
        dat[, c_counts_scaled := c_counts, by = channel]
      }
    } else {
      if (norm_method == "0-1") {
        dat[, c_counts_scaled := ((mc_counts - min(mc_counts)) / (max(mc_counts) - min(mc_counts))), by = channel]
        dat[c_counts_scaled < 0, c_counts_scaled := 0, by = channel]
      } else if (norm_method == "znorm") {
        dat[, c_counts_scaled := ((mc_counts - mean(mc_counts)) / sd(mc_counts)), by = channel]
      }
    }
  }

  a <- data.frame(dat$channel, dat$c_counts_scaled) %>%
    group_by(dat.channel) %>%
    dplyr::mutate(index = row_number()) %>%
    pivot_wider(
      names_from = dat.channel,
      values_from = dat.c_counts_scaled
    ) %>%
    dplyr::select(-index)

  histdata <- reshape2::melt(a, variable.name = "marker", value.name = "expression")
  p1 <- ggplot(data = histdata, aes(x = expression)) +
    geom_histogram(bins = 60, colour = "black", fill = "blue", alpha = 0.5) +
    facet_wrap(~marker, scale = "free")

  pdf("Marker Expression level.pdf", height = 8, width = 12)
  print(p1)
  dev.off()
  return(a)
}

censor_dat <- function(x, quant = 0.999, ignorezero = TRUE, symmetric = F) {
  if (symmetric) {
    lower_quant <- (1 - quant) / 2
    quant <- quant + lower_quant
  }

  if (ignorezero) {
    q <- stats::quantile(x[x > 0], quant)
  } else {
    q <- stats::quantile(x, quant)
  }
  x[x > q] <- q

  if (symmetric) {
    q <- stats::quantile(x, lower_quant)
    x[x < q] <- q
  }
  return(x)
}

## FloSOM clustering
runFlowsomPheno <- function(norm_exp_df, markers, xdim = 100, ydim = 100, phenok = 15) {
  set.seed(123)
  Do_FlowSOM <- TRUE
  if (Do_FlowSOM) {
    timestart <- Sys.time()

    #  (A).Input：
    #    1) 设定xdim和ydim的数值
    #      FloSOM第一步SOM聚类，需要预先指定cluster数目，最终会生成xdim*ydim个cluster，
    # xdim=100    #<-- SOM聚类会生成xdim*ydim个cluster
    # ydim=100   #<-- 建议xdim和ydim 取数值相等（或相近）的自然数


    #  (B).运行SOM聚类：
    FlowSOM_input_data <- norm_exp_df[, markers]
    print(dim(FlowSOM_input_data))
    map <- SOM(
      data = as.matrix(FlowSOM_input_data),
      xdim = xdim,
      ydim = ydim,
      silent = F
    )

    SOM_result <- map$mapping[, 1]
    print(length(unique(SOM_result)))

    #  (C).运行metacluster
    #     FlSOM第二步是对获得的SOM cluster再进行聚类，产生metacluster，数目可人为设定或者通过elbow_test得到

    FlowSOM_combined <- data.frame(
      FlowSOM = SOM_result,
      FlowSOM_input_data
    )

    metacluster_result <- metaClustering2(FlowSOM_combined,
      clustername = "FlowSOM",
      metaClustering_method = "metaClustering_PhenoGraph", #<-- 聚类方法
      k_value = phenok, #<-- 聚类方法中的k值
      elbow_test = F, #<-- 决定是否进行elbow test，当kvalue=NULL时，会自行选择进行elbowtest
      seed = 123
    ) #<-- 设定seed可以保持tsne图保持一致

    timeend <- Sys.time()
    runningtime <- timeend - timestart
    print(runningtime)
  }
  return(metacluster_result)
}

FlowSOM_clustering <- function(exp_df, norm_exp_df, markers, phenographOnly = F, xdim = 100, ydim = 100, Type, method, savePath, phenok = 30, using.kmenas = F) {
  if (using.kmenas) {
    df <- norm_exp_df[, markers]
    result_test <- kmeans(df, centers = phenok, iter.max = 1000000)
    norm_exp_df[paste0(Type, "_flowsom100pheno15")] <- result_test$cluster
    meta <- exp_df[setdiff(colnames(exp_df), colnames(norm_exp_df))]
    if (!dir.exists(savePath)) {
      dir.create(savePath)
    }
    write.csv(cbind(norm_exp_df, meta), file = paste0(savePath, Type, "_", method, "_flowsom_pheno_k", phenok, ".csv"), row.names = FALSE)
    return(norm_exp_df)
  }
  if (!phenographOnly) {
    result_test <- runFlowsomPheno(norm_exp_df, markers, xdim = xdim, ydim = ydim, phenok = phenok) # nolint
    norm_exp_df[paste0(Type, "_flowsom100pheno15")] <- (result_test[["metacluster_result"]])
    # tsneplot <- (result_test[["plotdf"]])
    meta <- exp_df[setdiff(colnames(exp_df), colnames(norm_exp_df))]
    if (!dir.exists(savePath)) {
      dir.create(savePath)
    }
    write.csv(cbind(norm_exp_df, meta), file = paste0(savePath, Type, "_", method, "_flowsom_pheno_k", phenok, ".csv"), row.names = FALSE)
    # write.csv(tsneplot, file = paste0(savePath, Type,"_",method, "_flowsom_pheno_tsne_k", phenok, ".csv"), row.names = FALSE)
  } else {
    tmp.pheno_out <- cytofkit::Rphenograph(norm_exp_df[, markers], k = phenok)
    norm_exp_df$cluster <- igraph::membership(tmp.pheno_out)
    meta <- exp_df[setdiff(colnames(exp_df), colnames(norm_exp_df))]
    if (!dir.exists(savePath)) {
      dir.create(savePath)
    }
    write.csv(cbind(norm_exp_df, meta), file = paste0(savePath, Type, "_", method, "_flowsom_pheno_k", phenok, ".csv"), row.names = FALSE)
  }

  print(paste0("FlowSOM clustering is done! savePath is ", savePath))
  return(norm_exp_df)
}

metaClustering2 <- function(indataframe = NULL,
                            clustername = NULL,
                            metaClustering_method = NULL,
                            usecol = NULL,
                            elbow_test = T,
                            k_value = NULL,
                            # tsne parameters：
                            view_tsne = T,
                            seed = NULL,
                            perplexity = 15,
                            max_iter = 1500,
                            ...) {
  cat("Get expression matrix of cluster centers...\n")
  if (is.null(usecol)) {
    usecol <- c(1:ncol(indataframe))
  }
  clusterparaid <- which(colnames(indataframe) == clustername)
  usecol <- union(usecol, c(clusterparaid))
  cluster_center_expr <- data.frame(indataframe[, usecol]) %>%
    dplyr::group_by_at(clustername) %>%
    dplyr::summarise_if(is.numeric, median)

  cluster_abundance <- data.frame(indataframe[, usecol]) %>%
    dplyr::group_by_at(clustername) %>%
    dplyr::summarise(num = n())

  cat("MetaClustering using method:", metaClustering_method, "...\n")
  if (elbow_test == T) {
    cat("Drawing elow curve...\n")
  }
  # Findvalue<-DetermineNumberOfClusters(data=cluster_center_expr,max=20,kstep = 2,method=metaClustering_method,plot = T)

  cc_metacluster <- MetaClustering(
    data = cluster_center_expr[, -clusterparaid],
    method = metaClustering_method,
    k_value = k_value,
    elbow_test = elbow_test
  )
  cat("\nMetaclustering cluster centers is finished.\n")
  cat("Start to mapping metaclusters to single cells...\n")

  Cluster_arrange <- data.frame(
    cluster = cluster_center_expr[, clustername],
    metacluster = cc_metacluster
  )

  cluster_arrange_fun <- function(cluster_id) {
    cellcluster <- subset(Cluster_arrange, Cluster_arrange[, 1] == cluster_id)$metacluster
    return(cellcluster)
  }
  metacluster_result <- apply(as.matrix(indataframe[, colnames(indataframe) == clustername]), 1, cluster_arrange_fun)

  if (view_tsne == T) {
    cat("Summarise metacluster information...\n")
    cat("Start to visualise metaclusters with tSNE...\n")
    # tsne (可以调节十多个参数，最重要的两个：perplexity 和 max_iter)
    if (is.null(seed)) {
      seed <- ceiling(runif(1, min = 0, max = 1) * 10000)
      cat("Seed is not specified, randomly set to: ", seed, "\n")
      set.seed(seed)
    } else {
      cat("Seed is set to: ", seed, ".\n")
      set.seed(seed)
    }


    tsne_result <- Rtsne(cluster_center_expr[, -c(1)],
      initial_dims = ncol(cluster_center_expr[, -c(1)]),
      dims = 2, check_duplicates = FALSE, pca = F, perplexity = 15, max_iter = 1500
    )$Y
    colnames(tsne_result) <- c("tsne_1", "tsne_2")

    combine_data_plot <- data.frame(cluster_center_expr,
      metacluster = cc_metacluster,
      tsne_result,
      num = cluster_abundance$num
    )

    centers <- combine_data_plot %>%
      dplyr::group_by(metacluster) %>%
      dplyr::summarise(tsne_1 = median(tsne_1), tsne_2 = median(tsne_2))


    ## visualization using ggplot2

    combine_data_plot$metacluster <- as.factor(combine_data_plot$metacluster)

    mytheme <- theme(
      panel.background = element_rect(fill = "white", colour = "black", size = 0.2), # 坐标系及坐标轴
      legend.key = element_rect(fill = "white", colour = "white"), # 图标
      legend.background = (element_rect(colour = "white", fill = "white"))
    )

    klab <- "Cluster number(k_value)"
    if (metaClustering_method == "metaClustering_PhenoGraph") klab <- "PhenoGraph_k(k_value)"
    ptsnemap <- ggplot(combine_data_plot) +
      geom_point(aes(x = tsne_1, y = tsne_2, colour = metacluster, size = num), alpha = 0.7) +
      guides(colour = guide_legend(ncol = 2, bycol = T)) +
      scale_size_continuous(range = c(0.1, 5)) +
      labs(title = paste0(klab, ": ", k_value)) +
      mytheme +
      geom_text(data = centers, aes(x = tsne_1, y = tsne_2), label = centers$metacluster, colour = "black", size = 5)

    png(filename = "./clustering_tsne.png")
    print(ptsnemap)
    dev.off()
  }

  cat("Metaclustering is finished successufully.\n")
  return(list(metacluster_result = metacluster_result, tsne = tsne_result, plotdf = combine_data_plot))
}

## Plot Heatmap
plotHeatmap <- function(a, marker_total, clustercol = "lcell_flowsom100pheno15", cutoff = FALSE, celltype = FALSE) {
  transdata <- a[marker_total]
  transdata["metacluster"] <- a[clustercol]
  colnames(transdata)

  tmp_per <- data.frame((table(transdata$metacluster) / dim(transdata)[1] * 100))
  tmp_per$percent <- round(tmp_per$Freq, 2)
  tmp_per <- tmp_per[, c("Var1", "percent")]
  tmp_per$type <- "percent"
  colnames(tmp_per)[2] <- "Freq"
  tmp_count <- data.frame((table(transdata$metacluster)))
  tmp_count$type <- "count"
  bardf <- rbind(tmp_per, tmp_count)

  g1 <- ggplot(bardf, aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity", width = 0.6, position = position_dodge(width = 0.8)) +
    geom_text(aes(label = Freq), position = position_dodge(width = 0.9), size = 4, vjust = 0.005, hjust = 0.005) +
    facet_wrap(~type, scales = "free") +
    coord_flip()

  print(g1)

  heatmap_data <- transdata %>%
    # dplyr::filter_at(vars(one_of(cluster_name)),all_vars(.%in% groups_to_show))%>%
    group_by_at(c("metacluster")) %>%
    summarise_if(is.numeric, mean, na.rm = TRUE) %>%
    data.frame(check.names = FALSE)

  if (cutoff != FALSE) {
    heatmap_data[, marker_total][(heatmap_data[, marker_total] > cutoff)] <- cutoff
  }

  htdf <- data.frame(t(heatmap_data[, marker_total]), check.names = FALSE)
  colnames(htdf) <- heatmap_data[, "metacluster"]
  p1 <- pheatmap(data.frame(htdf, check.names = FALSE),
    display_numbers = TRUE, cluster_rows = TRUE,
    clustering_method = "average", angle_col = "90"
  )
  return(p1)
}

## plot feature plot
plotFeaturePlot <- function(clusterCSVpath, markers = markers, savePath = savePath) {
  combine_data_plot <- read.csv(clusterCSVpath)
  head(combine_data_plot, 2)

  ## visualization using ggplot2

  combine_data_plot$metacluster <- as.factor(combine_data_plot$metacluster)

  mytheme <- theme(
    panel.background = element_rect(fill = "white", colour = "black", size = 0.2),
    legend.key = element_rect(fill = "white", colour = "white"),
    legend.background = (element_rect(colour = "white", fill = "white"))
  )

  p <- list()
  i <- 1

  for (marker in markers) {
    color <- as.numeric(unlist(combine_data_plot[marker]))
    color <- MaxMin_Nor(color)

    combine_data_plot_ <- combine_data_plot
    combine_data_plot_$color <- color

    ptsnemap <- ggplot(combine_data_plot_) +
      geom_point(aes(x = tsne_1, y = tsne_2, colour = color), size = 0.1) +
      guides(colour = guide_colorbar(title = "Expression")) +
      scale_size_continuous(range = c(0.1, 5)) +
      labs(title = paste0("Marker: ", marker)) +
      mytheme

    p[[i]] <- ptsnemap
    i <- i + 1
  }
  print(paste0("The number of marker is ", (i - 1)))
  if (i <= 16) {
    p_ <- plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], p[[7]], p[[8]], p[[9]], p[[10]], p[[11]], p[[12]], p[[13]], p[[14]], p[[15]], p[[16]],
      ncol = 4, nrow = 4
    )
    png(filename = paste0(savePath, "_1_featuresplot.png"), width = 960, height = 960)
    print(p_)
    dev.off()
  }

  if (i > 16 & i <= 32) {
    p_ <- plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], p[[7]], p[[8]], p[[9]], p[[10]], p[[11]], p[[12]], p[[13]], p[[14]], p[[15]], p[[16]],
      ncol = 4, nrow = 4
    )
    png(filename = paste0(savePath, "featuresplot_1.png"), width = 960, height = 960)
    print(p_)
    dev.off()

    p_ <- plot_grid(p[[17]], p[[18]], p[[19]], p[[20]], p[[21]], p[[22]], p[[23]], p[[24]], # p[[25]],p[[26]],p[[27]],p[[28]],p[[29]],p[[30]],p[[31]],p[[32]],
      ncol = 3, nrow = 3
    )
    png(filename = paste0(savePath, "featuresplot_2.png"), width = 960, height = 960)
    print(p_)
    dev.off()
  }
  return(NULL)
}

MaxMin_Nor <- function(vector) {
  max <- max(vector)
  min <- min(vector)

  return(vector - min) / (max - min)
}

## annotation
Annotate <- function(TsnecsvPath, ClustercsvPath, MajorTypeList, MinorTypeList, savePath) {
  num <- length(MajorTypeList)

  tsnecsv <- read.csv(TsnecsvPath)
  clustercsv <- read.csv(ClustercsvPath)

  tsnecsv$MajorType <- NA
  clustercsv$MajorType <- NA
  tsnecsv$SubType <- NA
  clustercsv$SubType <- NA

  for (i in 1:num) {
    clusterid <- as.character(i)
    tsnecsv[tsnecsv$metacluster == clusterid, "MajorType"] <- MajorTypeList[i]
    clustercsv[clustercsv$flowsom100pheno15 == clusterid, "MajorType"] <- MajorTypeList[i]
  }

  print(paste0("The annotate result of major type is:\n"))
  print(table(clustercsv["MajorType"]))

  for (i in 1:num) {
    clusterid <- as.character(i)
    tsnecsv[tsnecsv$metacluster == clusterid, "SubType"] <- MinorTypeList[i]
    clustercsv[clustercsv$flowsom100pheno15 == clusterid, "SubType"] <- MinorTypeList[i]
  }

  print(paste0("The annotate result of major types is:\n"))
  print(table(clustercsv["SubType"]))

  saveRDS(tsnecsv, paste0(savePath, "annotate_1wcells.rds"))
  saveRDS(clustercsv, paste0(savePath, "annotate_allcells.rds"))

  ## Plot the major celltypes of 10k cells
  tsnecsv$MajorType <- as.factor(tsnecsv$MajorType)

  mytheme <- theme(
    panel.background = element_rect(fill = "white", colour = "black", size = 0.2),
    legend.key = element_rect(fill = "white", colour = "white"),
    legend.background = (element_rect(colour = "white", fill = "white"))
  )

  ptsnemap <- ggplot(tsnecsv) +
    geom_point(aes(x = tsne_1, y = tsne_2, colour = MajorType), size = 0.2) +
    guides(colour = guide_legend(ncol = 2, bycol = T)) +
    scale_size_continuous(range = c(0.05, 0.2)) +
    labs(title = paste0("T-SNE of Major type")) +
    mytheme

  pdf(file = paste0(savePath, "Majortype_clustering.pdf"), width = 8, height = 6)
  print(ptsnemap)
  dev.off()

  ## Plot the cell subtypes of 10k cells
  tsnecsv$SubType <- as.factor(tsnecsv$SubType)

  mytheme <- theme(
    panel.background = element_rect(fill = "white", colour = "black", size = 0.2),
    legend.key = element_rect(fill = "white", colour = "white"),
    legend.background = (element_rect(colour = "white", fill = "white"))
  )

  ptsnemap <- ggplot(tsnecsv) +
    geom_point(aes(x = tsne_1, y = tsne_2, colour = SubType), size = 0.2) +
    guides(colour = guide_legend(ncol = 2, bycol = T)) +
    scale_size_continuous(range = c(0.1, 5)) +
    labs(title = paste0("T-SNE of subpopulation")) +
    mytheme

  pdf(file = paste0(savePath, "Subtype_clustering.pdf"), width = 10, height = 7.5)
  print(ptsnemap)
  dev.off()

  return(NULL)
}


## markers expression heatmap
MarkerHeatmap <- function(data, markers, savePath) {
  dataBack <- data

  MinorType <- names(table(data$SubType))
  MajorType <- c()

  for (minortype in MinorType) {
    x <- data[data$SubType == minortype, ]
    MajorType <- c(MajorType, x$MajorType[1])
  }

  typedf <- data.frame("MajorType" = MajorType, "SubType" = MinorType)
  typedf <- typedf[order(typedf[, 1]), ]
  rownames(typedf) <- 1:nrow(typedf)
  data <- data[, markers]

  plotdf <- data.frame(matrix(data = NA, nrow = length(typedf$SubType), ncol = length(markers)))
  rownames(plotdf) <- typedf$SubType
  colnames(plotdf) <- markers

  plotdf$fraction <- NA

  j <- 1
  for (i in typedf$SubType) {
    index <- (dataBack$SubType == i)
    x <- data[index, ]

    plotdf$fraction[j] <- paste0(typedf$SubType[j], "(", as.character(round(dim(x)[1] / dim(data)[1], 3)), ")")

    x <- apply(x, 2, mean)

    plotdf[j, 1:(ncol(plotdf) - 1)] <- x
    j <- j + 1
  }

  plotdf$MajorType <- typedf$MajorType
  plotdf$MinorType <- typedf$MinorType
  rownames(plotdf) <- plotdf$fraction

  annotationRow <- as.data.frame(plotdf[, ncol(plotdf)])
  colnames(annotationRow) <- "MajorType"
  ann_colors <- list(MajorType = c(Lymphocyte = "#66C2A5", Myeloid = "#FC8D62", Stromal = "#8DA0CB", Tumor = "#E78AC3", UNKNOWN = "#A6D854"))

  p <- pheatmap(plotdf[, 1:length(markers)],
    cluster_rows = F, # gaps_row = c(8, 25, 30, 34),
    cluster_cols = F, # gaps_col = c(7, 16, 20, 21),
    color = colorRampPalette(colors = c("white", "red"))(100),
    annotation_row = annotationRow,
    annotation_colors = ann_colors
  )

  pdf(file = paste0(savePath, "Cluster Markers Heatmap.pdf"), width = 12, height = 8)
  print(p)
  dev.off()

  return(NULL)
}

## Heatmap for step 7 reclustering
getMeanExp <- function(sce, label) {
  types <- names(table(colData(sce)[label]))

  ## result matrix
  meanExp <- matrix(data = NA, nrow = length(rownames(sce)), ncol = 0)
  meanExp <- as.data.frame(meanExp)
  rownames(meanExp) <- rownames(sce)

  fraction <- c()
  for (i in types) {
    sceTemp <- sce[, unlist(colData(sce)[label] %in% i)]
    expTemp <- assays(sceTemp)[[1]]
    expTemp <- apply(expTemp, MARGIN = 1, FUN = "mean")
    meanExp[, i] <- expTemp

    fraction <- c(fraction, round(ncol(sceTemp) / ncol(sce), 4))
  }
  colnames(meanExp) <- paste0(colnames(meanExp), " (", as.character(fraction), ")")

  return(meanExp)
}

SubtypeHeatmap <- function(sce, label, savePath) {
  meanExp <- getMeanExp(sce, label)
  meanExp <- t(meanExp)

  p <- pheatmap(
    meanExp,
    display_numbers = TRUE, cluster_rows = TRUE, cluster_cols = FALSE,
    clustering_method = "average", angle_col = "90"
  )

  pdf(savePath, width = 12, height = 3)
  print(p)
  dev.off()

  return(NULL)
}

## Access the clustering result
AssessClustering <- function(exp, markers, clusterCol, tFraction = 0.8, sampleSize = 50000, savePath) {
  set.seed(619)

  exp <- exp[sample.int(nrow(exp), size = sampleSize), ]

  y <- exp[, clusterCol]
  x <- exp[, markers]

  cat("The category are ", clusterCol, "\n")

  train <- sample.int(nrow(x), size = as.integer(nrow(x) * tFraction))
  Xtrain <- x[train, ]
  Xtest <- x[-train, ]

  Ytrain <- as.factor(y[train])
  Ytest <- y[-train]

  dfTrain <- cbind(Xtrain, "Label" = as.factor(Ytrain))

  RF <- randomForest(Label ~ ., data = dfTrain, importance = TRUE)
  predict <- predict(RF, Xtest)

  ## confusion matrix
  TrainConfuMat <- RF[["confusion"]]
  TrainConfuMat <- TrainConfuMat[, -ncol(TrainConfuMat)]
  TrainConfuMat <- apply(TrainConfuMat, MARGIN = 2, function(x) {
    return(round(x / sum(x), digits = 3))
  })

  TestConfuMat <- confusionMatrix(predict, as.factor(Ytest))[["table"]]
  TestConfuMat <- apply(TestConfuMat, MARGIN = 2, function(x) {
    return(round(x / sum(x), digits = 3))
  })

  ## visualize
  pTrain <- pheatmap(TrainConfuMat,
    color = brewer.pal(9, "Greens"),
    cellwidth = 25, cellheight = 15,
    cluster_row = F, cluster_col = F,
    angle_col = "90", display_numbers = TRUE, fontsize_number = 8
  )

  pTest <- pheatmap(TestConfuMat,
    color = brewer.pal(9, "Greens"),
    cellwidth = 25, cellheight = 15,
    cluster_row = F, cluster_col = F,
    angle_col = "90", display_numbers = TRUE, fontsize_number = 8
  )

  pdf(paste0(savePath, "Confusion Matrix of Traning.pdf"), width = 12, height = 8)
  print(pTrain)
  dev.off()
  pdf(paste0(savePath, "Confusion Matrix of Test.pdf"), width = 12, height = 8)
  print(pTest)
  dev.off()
  return(NULL)
}
