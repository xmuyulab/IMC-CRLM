## building signature
library(SingleCellExperiment)
library(pheatmap)

source("/home/lyx/project/IMC/structural_analysis_functions.r")
source("./signature_functions.r")

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
names(table(sce$Tissue))

## omit Tissue associated Tumor
sce <- sce[, sce$Tissue == "IM"]
sce

savePath <- paste0("/mnt/data/lyx/IMC/analysis/signature/")
if (!dir.exists(savePath)) {
    dir.create(savePath, recursive = T)
}

## Using Cell subpopulation fraction ratio to regression
clinicalFeatures <- c("RFS_time", "RFS_status")
countdf <- GetAbundance(sce, "SubType", is.fraction = T, is.reuturnMeans = T)
clinical <- countdf[, c(29, 28)]
countdf <- countdf[, -c(26:30)]

countdf <- cbind(countdf, clinical)

featureTypes <- colnames(countdf)[1:24]

RatioMat <- TranferAbun2Ratio(countdf, featureTypes)

## Uni-variable cox to fileter feature
RatioMat <- UniCoxFilter(RatioMat, do.padjust = F)

set.seed(619)
## spliting samples
list_ <- SplitData(RatioMat, trainginSize = 1, do.scale = T)

## lasso-cox
model <- FitLassoCox(list_, savePath)

## analysis the features
AnalysisFeature(model, savePath)

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

## using xgboost to seperate patients into relapse and non-relapse
