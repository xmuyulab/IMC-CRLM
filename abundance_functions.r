# functions for cell type abundance analysis
library(ggplot2)
library(ggrepel)
library(ggsci)
library(ggpubr)
library(survival)
library(tableone)
library(forestplot)
library(survminer)
library(SingleCellExperiment)
library(ggcor)

## transform annotated data into SingleCellExperiment object
SCE_Transform <- function(scedf, assay_col = c(1, 35), cellmeta_col = c(36, 49), clinical = NULL) {
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

## transfer sce to cell counts matrix
Transform_CellCountMat <- function(sceobj, group = c("IM", "CT"), clinicalFeatures, is.fraction = FALSE) {
  if (length(group) == 2) {
    sceobj <- sceobj[, sceobj$Tissue == group[1] | sceobj$Tissue == group[2]]
  }
  if (length(group) == 3) {
    sceobj <- sceobj[, sceobj$Tissue == group[1] | sceobj$Tissue == group[2] | sceobj$Tissue == group[3]]
  }

  cellMeta <- colData(sceobj)

  ## ROI, major celltype and cell subtype names and other clinical information
  ROIs <- names(table(cellMeta$filelist))

  # MajorTypes <- names(table(cellMeta$MajorType))
  SubTypes <- names(table(cellMeta$SubType))
  # alltypes <- unique(c(MajorTypes, SubTypes))
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
    # MajorTem <- as.data.frame(t(table(coldataTemp$MajorType)))
    # MajorTem$Freq <- round(MajorTem$Freq / cellnum, digits = 4)

    SubTem <- as.data.frame(t(table(coldataTemp$SubType)))
    # SubTem$Freq <- round(SubTem$Freq / sum(SubTem$Freq), digits = 4)

    # CellCountMat[match(ROI, rownames(CellCountMat)), match(MajorTem$Var2, colnames(CellCountMat))] <- MajorTem$Freq
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
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "Worse Prog versus Better"
    }
    if (feature == "RFS_status") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "Relapse versus Non-relapse"
    }
    if (feature == "Recurrence_site") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "Extrahepatic relapse versus Liver"
    }
    if (feature == "fong_score") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "High-FongScore versus Low-FongScore"
    }
    if (feature == "Age") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- ">=60 versus <60"
    }
    if (feature == "Gender") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "Female versus Male"
    }
    if (feature == "KRAS_mutation") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "KRAS Mut versus WT"
    }
    if (feature == "BRAF_mutation") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "BRAF Mut versus WT"
    }
    if (feature == "mTBS") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "High-mTBS versus Low-mTBS"
    }
    if (feature == "CRLM_number") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "High-CRLMNUM versus Low-CRLMNUM"
    }
    if (feature == "CRLM_size") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "High-CRLMSize versus Low-CRLMSize"
    }
    if (feature == "Live_involvement_num") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "Liver involvement 2 versus involvement 1"
    }
    if (feature == "CEA") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "Abnormal CEA versus Normal"
    }
    if (feature == "CA199") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "Abnormal CA199 versus Normal"
    }
    if (feature == "Pathology") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "Mucinous versus Adenocarcinoma"
    }
    if (feature == "Lymph_grade") {
      DF[which(DF$Feature == feature), ]$FeatureGroup <- "N-positive versus N-negative"
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
MultiCliDotplot <- function(plotdf, tissue, savePath){
  plotdf$P.label <- ifelse(plotdf$P.value <= 0.001, "***", ifelse(plotdf$P.value <= 0.01, "**", ifelse(plotdf$P.value <= 0.05, "*", "n.s.")))
  plotdf$lgPvalue <- ifelse(plotdf$P.value <= 0.001,4,ifelse(plotdf$P.value <= 0.01, 3, ifelse(plotdf$P.value <= 0.05, 2, 1)))
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
    scale_color_gradient2(low = "royalblue",mid = "white", high = "firebrick3")

    pdf(paste0(savePath,"Multiple clinical features abundance differnt of ",tissue,".pdf"), height = 6, width = 10)
    print(p)
    dev.off()
  return(NULL)
}

## abundance boxplot
abundanceBoxplot <- function(plotdf, celltypes2Plot, expCol = c(1, 27)) {
  plotdf2 <- CountMat2Plotdf(plotdf, expCol)
  plotdf3 <- plotdf2[plotdf2$Celltype %in% celltypes2Plot, ]

  p <- ggplot(data = plotdf3, aes(x = Tissue, y = Abundance, fill = Recurrence)) +
    geom_boxplot(alpha = 0.7) +
    scale_y_continuous(name = "Cell Abundance") +
    scale_x_discrete(name = "Cell Population") +
    ggtitle("Boxplot of cell type abundance") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      text = element_text(size = 12),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(size = 11, angle = 90)
    ) +
    facet_wrap(~Celltype, ncol = 3) +
    scale_fill_lancet() +
    stat_compare_means(aes(group = Recurrence), label.y = 0.5, method = "t.test")

  pdf("/mnt/data/lyx/IMC/analysis/abundance/abundance_analysis(rec).pdf", width = 15, height = 8)
  print(p)
  dev.off()

  if (F) {
    p <- ggplot(data = plotdf2, aes(x = Prognosis, y = Abundance, fill = Prognosis)) +
      geom_boxplot(alpha = 0.7) +
      scale_y_continuous(name = "Cell Abundance") +
      scale_x_discrete(name = "Cell Population") +
      ggtitle("Boxplot of cell type abundance") +
      theme_bw() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11, angle = 90)
      ) +
      facet_wrap(~Celltype) +
      scale_fill_lancet() +
      stat_compare_means(aes(group = Prognosis), p.adjust.method = "BH", label.y = 0.6)

    pdf("/mnt/data/lyx/IMC/abundance/abundance_analysis(prog).pdf", width = 12, height = 9)
    print(p)
    dev.off()
  }

  return(NULL)
}

## Convert count matrix into plot dataframe
CountMat2Plotdf <- function(plotdf, expCol = expCol) {
  ### matrix transform
  sum_ <- c()

  ## calculate fraction
  for (i in 1:nrow(plotdf)) {
    sumTemp <- sum(plotdf[i, (expCol[1]):expCol[2]])
    sum_ <- c(sum_, sumTemp)
  }
  for (i in 1:nrow(plotdf)) {
    plotdf[i, (expCol[1]):expCol[2]] <- round(plotdf[i, (expCol[1]):expCol[2]] / sum_[i], digits = 5)
  }
  plotdf2 <- data.frame(matrix(data = NA, nrow = nrow(plotdf) * expCol[2]))

  abundanceVec <- c()
  typeVec <- rep(colnames(plotdf)[expCol[1]:expCol[2]], each = nrow(plotdf))
  TissueVec <- rep(plotdf$Tissue, times = ncol(plotdf2))
  # progVec <- rep(plotdf$Prognosis, times = ncol(plotdf2))
  RFSsTimeVec <- rep(plotdf$RFS_time, times = ncol(plotdf2))
  RFSsVec <- rep(plotdf$RFS_status, times = ncol(plotdf2))


  for (i in 1:expCol[2]) {
    abundanceVec <- c(abundanceVec, as.numeric(plotdf[, i]))
  }

  plotdf2$Abundance <- abundanceVec
  plotdf2$Celltype <- typeVec
  plotdf2$Tissue <- TissueVec
  # plotdf2$Prognosis <- progVec
  plotdf2$RecurrenceTime <- RFSsTimeVec
  plotdf2$Recurrence <- RFSsVec
  plotdf2 <- plotdf2[, -1]

  plotdf2$Recurrence <- ifelse(plotdf2$Recurrence == 1, "Relapse", "Non Relapse")
  # plotdf2$Prognosis <- ifelse(plotdf2$Prognosis == 1, "Worse", "Better")

  plotdf2$Recurrence <- as.factor(plotdf2$Recurrence)
  # plotdf2$Prognosis <- as.factor(plotdf2$Prognosis)

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
MultipleUniCOX <- function(df) {
  result <- matrix(data = NA, nrow = 0, ncol = 6)
  result <- as.data.frame(result)

  features <- colnames(df)[1:(ncol(df) - 2)]
  cat("The features in multi-cox are: ", features, "\n")

  univ_formulas <- sapply(
    features,
    function(x) {
      as.formula(paste("Surv(RFS_time, RFS_status)~", x))
    }
  )

  univ_models <- lapply(univ_formulas, function(x) {
    coxph(x, data = df)
  })

  univ_results <- lapply(univ_models, function(x) {
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
}

abundanceMetaAnalysis <- function(plotdf, celltypes2Plot, clinical, features) {
  ## PIDs
  plotdf <- as.data.frame(plotdf)
  plotdf$PID <- rownames(plotdf)

  clinical2 <- clinical[match(plotdf$PID, clinical$PID), features]

  clinical_features <- colnames(clinical2)
  print("Clinical Features: ")
  print(clinical_features)

  rownames(clinical2) <- c(1:nrow(clinical2))

  ## bind more clinical information
  plotdf$zs_rec_riskmodel <- clinical2$zs_rec_riskmodel
  plotdf$fong_score <- clinical2$fong_score
  plotdf$Gender <- clinical2$Gender
  plotdf$Age <- clinical2$Age
  plotdf$KRAS_mutation <- clinical2$KRAS_mutation
  plotdf$BRAF_mutation <- clinical2$BRAF_mutation
  plotdf$mTBS <- clinical2$mTBS
  plotdf$CRLM_number <- clinical2$CRLM_number
  plotdf$CRLM_size <- clinical2$CRLM_size
  plotdf$Live_involvement_num <- clinical2$Live_involvement_num
  plotdf$Pathology <- clinical2$Pathology
  plotdf$Differential_grad <- clinical2$Differential_grad
  plotdf$T_grade <- clinical2$T_grade
  plotdf$RFS_time <- clinical[match(plotdf$PID, clinical$PID), "RFS_time"]
  plotdf$RFS_status <- clinical[match(plotdf$PID, clinical$PID), "RFS_status"]

  # cat("Number NA is:", table(is.na(plotdf)), "\n")

  plotdf <- plotdf[, c(celltypes2Plot, features, "RFS_time", "RFS_status")]

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

  pdf("test.pdf", width = 12, height = 9)
  print(p)
  dev.off()

  return(NULL)
}

## Kaplan-Meier curve
KMVisualize <- function(df, celltype, cutoff = "median", savePath = NULL) {
  df <- df[, c(celltype, "RFS_status", "RFS_time")]

  if (cutoff == "median") {
    df[, celltype] <- ifelse(df[, 1] >= median(df[, 1]), "high", "low")
  }
  if (cutoff == "mean") {
    df[, celltype] <- ifelse(df[, 1] >= mean(df[, 1]), "high", "low")
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

  pdf(savePath, height = 12, width = 12)
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
