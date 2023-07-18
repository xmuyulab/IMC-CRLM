## building signature
library(SingleCellExperiment)
library(pheatmap)

source("/home/lyx/project/IMC/structural_analysis_functions.r")
source("./signature_functions.r")

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce(new).rds")

table(sce$SubType)

savePath <- paste0("/mnt/data/lyx/IMC/analysis/signature/")
if (!dir.exists(savePath)) {
    dir.create(savePath, recursive = T)
}

### Using Cell subpopulation fraction ratio to build lasso-cox model
clinicalFeatures <- c("RFS_time", "RFS_status", "KRAS_mutation")
countdf <- GetAbundance(sceobj = sce, countcol = "SubType", clinicalFeatures = clinicalFeatures, is.fraction = T, is.reuturnMeans = T)
colnames(countdf)
clinical <- countdf[, c(25:27)]
countdf <- countdf[, -c(23:ncol(countdf))]

# featureTypes <- colnames(countdf)[1:24]
FeatureType1 <- c("B", "CD4T", "CD8T")
FeatureType2 <- c("Macro_CD163", "PDL1+ Mono_CD11c", "Mono_CD11c", "PD1+ Treg", "Treg")

# allTypes <- colnames(countdf)
# RatioMat <- TranferAbun2Ratio(countdf,CombFeatureTypes = allTypes)

RatioMat <- TranferAbun2Ratio(countdf, FeatureType1 = FeatureType1, FeatureType2 = FeatureType2)

RatioMat$ID <- rownames(RatioMat)
clinical$ID <- rownames(clinical)

RatioMat2 <- dplyr::left_join(RatioMat, clinical, by = "ID")
rownames(RatioMat2) <- rownames(RatioMat)

### Boxplot to show the cell subpopulation fraction difference in certain clinical group
colnames(RatioMat2)
p <- RatioBoxPlot(RatioMat2, expCol = c(1:30), clinicalGroupCol = "RFS_status", ClinicalGroupName = c("R", "NR"), do.scale = T, FilterNotSig = F)
pdf(paste0(savePath, "SubPopulation abundance Ratio Boxplot (Treg_CD11c).pdf"), width = 12, height = 6)
print(p)
dev.off()

## correlation between survival time and cell subpopulation ratio (Not sig)
if (F) {
    plotType <- c("CD4T:Macro_CD163", "CD8T:Macro_CD163", "CD4T:PDL1+ Mono_CD11c")
    p <- RatioCorPlot(RatioMat2, expCol = c(1:30), clinicalGroupCol = "RFS_status", ClinicalGroupName = c("R", "NR"), plotType = plotType, do.scale = T)
    pdf(paste0(savePath, "SubPopulation abundance Ratio Linechart (Treg_CD11c).pdf"), width = 12, height = 6)
    print(p)
    dev.off()
}


### Using uni-variable cox to filter out some features
RatioMat <- UniCoxFilter(RatioMat, do.padjust = F, do.scale = T)

set.seed(619)
list_ <- SplitData(RatioMat, trainginSize = 1, do.scale = T) ## spliting samples

df <- FitLassoCoxBoost(list_, times = 100) ## lasso-cox with boost

p <- AnalysisFeatureAfterBoost(df) ## analysis features
pdf("/mnt/data/lyx/IMC/analysis/signature/lassocox_boostrap features (PD1+ Treg)(PID).pdf", height = 6, width = 10)
print(p)
dev.off()

model <- FitLassoCox(list_, savePath) ## build final model
AnalysisFeature(model, savePath) ## analysis the features

## Predict
pred_train <- Inference(list_[["x_train"]], list_[["y_train"]], model)
pred_test <- Inference(list_[["x_test"]], list_[["y_test"]], model)

## K-M plot
TrainCPList <- KMplot(list_[["y_train"]], pred_train, set = "Training", savePath)
TestCPList <- KMplot(list_[["y_test"]], pred_test, set = "Test", savePath)

p.test <- TestCPList[["logrank_p"]]
cindex <- TestCPList[["c-index"]]

if (p.test <= 0.05) {
    cat("In iteration ", i, " the p-value in test set is ", p.test, " the c-index is ", cindex, "\n")
    list_Temp <- list()
    list_Temp[["list_"]] <- list_
    list_Temp[["model"]] <- model

    saveRDS(list_Temp, paste0(savePath, "all_list.rds"))
} else {
    system(paste0("rm -rf ", savePath))
}

### Using marker expression intensity to build lasso-cox model
if (F) {
    markers <- c("CD4", "CD3", "CD8a", "CD163", "CD279")
    list_ <- GetMeanExpressionProfile(sce, LevelCol = "PID", markers = markers, clinicalGroupCol = c("RFS_time", "RFS_status"))
    df <- list_[["df"]]
    clinicalGroup <- list_[["clinicalGroup"]]

    list_ <- SplitData(cbind(df, clinicalGroup), trainginSize = 1, do.scale = T)
    model <- FitLassoCox(list_, savePath) ## lasso-cox
    AnalysisFeature(model, savePath) ## analysis the features

    ## Predict
    pred_train <- Inference(list_[["x_train"]], list_[["y_train"]], model)

    ## K-M plot
    TrainCPList <- KMplot(list_[["y_train"]], pred_train, set = "Training", savePath)
    p.test <- TrainCPList[["logrank_p"]]
    cindex <- TrainCPList[["c-index"]]

    ## lasso cox with boostrapping
    ## Patients level
    list_ <- GetMeanExpressionProfile(sce, LevelCol = "PID", markers = markers, clinicalGroupCol = c("RFS_time", "RFS_status"))
    df <- list_[["df"]]
    clinicalGroup <- list_[["clinicalGroup"]]
    list_ <- SplitData(cbind(df, clinicalGroup), trainginSize = 1, do.scale = T)
    df <- FitLassoCoxBoost(list_, times = 5000)
    p <- AnalysisFeatureAfterBoost(df)
    pdf("/mnt/data/lyx/IMC/analysis/signature/lassocox_boostrap features (Patient).pdf", height = 6, width = 10)
    print(p)
    dev.off()

    ## Final model
    markers <- c(
        "CD279", "CD163", "CD11c", ## risk
        "CD80", "CD16", "CD4" ## protect
    )
    list_ <- GetMeanExpressionProfile(sce, LevelCol = "ID", markers = markers, clinicalGroupCol = c("RFS_time", "RFS_status"))
    df <- list_[["df"]]
    clinicalGroup <- list_[["clinicalGroup"]]
    list_ <- SplitData(cbind(df, clinicalGroup), trainginSize = 1, do.scale = T)
    model <- FitLassoCox(list_, savePath)

    ## analysis the features
    AnalysisFeature(model, savePath)
    ## Predict
    pred_train <- Inference(list_[["x_train"]], list_[["y_train"]], model)

    ## K-M plot
    TrainCPList <- KMplot(list_[["y_train"]], pred_train, set = "Training", savePath)
    p.test <- TrainCPList[["logrank_p"]]
    cindex <- TrainCPList[["c-index"]]
}


## Using ratio lasso-cox model to predict the clinical groups
clinicalFeatures <- c("RFS_time", "RFS_status", "KRAS_mutation", "fong_score", "Age", "TBS", "CRLM_number", "CRLM_size", "CEA", "CA199", "T_stage")
countdf <- GetAbundance(sceobj = sce, countcol = "SubType", clinicalFeatures = clinicalFeatures, is.fraction = T, is.reuturnMeans = T)

colnames(countdf)
clinical <- countdf[, c(25:ncol(countdf))]
countdf <- countdf[, -c(23:ncol(countdf))]

FeatureType1 <- c("B", "CD4T", "CD8T")
FeatureType2 <- c("Macro_CD163", "PDL1+ Mono_CD11c", "Mono_CD11c")

RatioMat <- TranferAbun2Ratio(countdf, FeatureType1 = FeatureType1, FeatureType2 = FeatureType2)
# RatioMat <- TranferAbun2Ratio(countdf, CombFeatureTypes = colnames(countdf))

RatioMat$ID <- rownames(RatioMat)
clinical$ID <- rownames(clinical)

RatioMat2 <- dplyr::left_join(RatioMat, clinical[, c("ID", "RFS_time", "RFS_status")], by = "ID")
RatioMat2 <- RatioMat2[, -which(colnames(RatioMat2) == "ID")]
rownames(RatioMat2) <- rownames(RatioMat)

list_ <- SplitData(RatioMat2, trainginSize = 1, do.scale = T, selectingFeature = c("CD4T:Macro_CD163", "CD4T:Mono_CD11c", "Macro_CD163:CD4T", "Macro_CD163:CD8T", "PDL1+ Mono_CD11c:CD4T", "PDL1+ Mono_CD11c:B"))

set.seed(619)

for (i in 1:1000) {
    model <- tryCatch(
        {
            FitLassoCox(list_, savePath)
        },
        error = function(e) {
            return("ERROR")
        }
    )

    if ((model == "ERROR")[1]) {
        next
    }
    cat("Iteration ", i, " cindex = ", model[["cindex"]], "\n")

    if (model[["cindex"]] >= 0.69) {
        saveRDS(model, "/mnt/data/lyx/IMC/analysis/signature/model.rds")

        ## analysis the features
        AnalysisFeature(model, savePath)

        ## Predict
        pred_train <- Inference(list_[["x_train"]], list_[["y_train"]], model)
        # pred_test <- Inference(list_[["x_test"]], list_[["y_test"]], model)

        ## K-M plot
        TrainCPList <- KMplot(list_[["y_train"]], pred_train[["riskLabel"]], set = "Training", savePath)

        Temp <- KMplot(list_[["y_train"]], as.numeric(pred_train[["riskScore"]]), set = "Training", plot = F)
        cindex_ <- Temp[["c-index"]]
        pvalue <- Temp[["logrank_p"]]

        # TestCPList <- KMplot(list_[["y_test"]], pred_test[["riskLabel"]], set = "Test", savePath)
        # p.test <- TestCPList[["logrank_p"]]
        # cindex <- TestCPList[["c-index"]]

        ## combine features
        # result <- MultipleUniCOX(compareDF[, -match("ID", colnames(compareDF))])
        # write.table(result,paste0(savePath,"MultiCOX_compare.csv"),sep = ',')
        # mul_cox <- coxph(Surv(RFS_time, RFS_status) ~ riskScore + riskLabel + fong_score + Age + TBS + CRLM_number + CRLM_size + CEA + CA199 + T_stage,data = compareDF)
        break
    }
}

## Compare features
model <- readRDS("/mnt/data/lyx/IMC/analysis/signature/model.rds")

### Predict
pred_train <- Inference(list_[["x_train"]], list_[["y_train"]], model)

## combine features
compareDF <- cbind("riskScore" = pred_train[[1]], clinical)
pcoxDF <- matrix(data = NA, nrow = 0, ncol = 3)
pcoxDF <- as.data.frame(pcoxDF)
compareFeatures <- colnames(compareDF)[c(1, 4:12)]

### Certain cell subpopulation abundance
countdf$ID <- rownames(countdf)
compareDF <- dplyr::left_join(compareDF, countdf[, unique(c(FeatureType1, FeatureType2, "ID"))], by = "ID")

### Ratio
RatioMat2$ID <- rownames(RatioMat2)
selectingRatio <- c("CD4T:Macro_CD163", "Macro_CD163:CD4T", "Macro_CD163:CD8T", "PDL1+ Mono_CD11c:B", "PDL1+ Mono_CD11c:CD4T")
compareDF <- dplyr::left_join(compareDF, RatioMat2[, unique(c(selectingRatio, "ID"))], by = "ID")

### Intensity
df <- as.data.frame(df)
df$ID <- rownames(df)
compareDF <- dplyr::left_join(compareDF, df[, unique(c(colnames(df), "ID"))], by = "ID")

compareDF <- compareDF[, !(colnames(compareDF) == "ID")]

## Get C-index
pcoxDF <- matrix(data = NA, nrow = 0, ncol = 3)
pcoxDF <- as.data.frame(pcoxDF)

compareFeatures <- colnames(compareDF)[c(1, 4:ncol(compareDF))]
for (feature in compareFeatures) {
    Temp <- KMplot(list_[["y_train"]], as.numeric(compareDF[, feature]), set = "Training", plot = F)

    Tempdf <- c(feature, Temp[[1]], as.numeric(Temp[[2]]))
    pcoxDF <- rbind(pcoxDF, Tempdf)
}

colnames(pcoxDF) <- c("Features", "Pvalue", "Cindex")

## visualize
p <- LollipopForMultiFeatures(pcoxDF)
pdf(paste0(savePath, "Multifeature Lollipop plots.pdf"), width = 6, height = 8)
print(p)
dev.off()
