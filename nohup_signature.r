## building signature
library(SingleCellExperiment)
library(pheatmap)

source("/home/lyx/project/IMC/structural_analysis_functions.r")
source("./signature_functions.r")

if (F) {
    ## IM
    if (T) {
        sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce(Treg).rds")
        names(table(sce$Tissue))
        table(sce$SubType)

        savePath <- paste0("/mnt/data/lyx/IMC/analysis/signature/IM/")
        if (!dir.exists(savePath)) {
            dir.create(savePath, recursive = T)
        }

        ## Using Cell subpopulation fraction ratio to regression
        clinicalFeatures <- c("RFS_time", "RFS_status")
        countdf <- GetAbundance(sce, "SubType", is.fraction = T, is.reuturnMeans = T)
        clinical <- countdf[, c(26, 25)]
        countdf <- countdf[, -c(23:ncol(countdf))]

        countdf <- cbind(countdf, clinical)

        FeatureType1 <- c("B", "CD4T", "CD8T")
        FeatureType2 <- c("Macro_CD163", "PD1+ Treg", "Treg")

        RatioMat <- TranferAbun2Ratio(countdf, CombFeatureTypes = c(FeatureType1, FeatureType2))

        set.seed(619)
        ## spliting samples
        list_ <- SplitData(RatioMat, trainginSize = 1, do.scale = T)

        ## lasso-cox
        df <- FitLassoCoxBoost(list_, times = 1000)
        p <- AnalysisFeatureAfterBoost(df)
        pdf(paste0(savePath, "lassocox_boostrap features (PD1+ Treg)(PID).pdf"), height = 6, width = 10)
        print(p)
        dev.off()

        ## ROI level
        ## Using Cell subpopulation fraction ratio to regression
        countdf <- GetAbundance(sce, "SubType", is.fraction = T, is.reuturnMeans = F)
        clinical <- countdf[, c(26, 25)]
        countdf <- countdf[, -c(23:ncol(countdf))]

        countdf <- cbind(countdf, clinical)
        RatioMat <- TranferAbun2Ratio(countdf, CombFeatureTypes = c(FeatureType1, FeatureType2))

        set.seed(619)
        ## spliting samples
        list_ <- SplitData(RatioMat, trainginSize = 1, do.scale = T)

        ## lasso-cox
        df <- FitLassoCoxBoost(list_, times = 1000)
        p <- AnalysisFeatureAfterBoost(df)
        pdf(paste0(savePath, "lassocox_boostrap features (PD1+ Treg)(ROI).pdf"), height = 6, width = 10)
        print(p)
        dev.off()

        ## Alltypes
        clinicalFeatures <- c("RFS_time", "RFS_status")
        countdf <- GetAbundance(sce, "SubType", is.fraction = T, is.reuturnMeans = T)
        clinical <- countdf[, c(26, 25)]
        countdf <- countdf[, -c(23:ncol(countdf))]

        countdf <- cbind(countdf, clinical)
        CombFeatureTypes <- colnames(countdf)[1:21]

        RatioMat <- TranferAbun2Ratio(countdf, CombFeatureTypes = CombFeatureTypes)
        ## Uni-variable cox to fileter feature
        RatioMat <- UniCoxFilter(RatioMat, do.padjust = F, do.scale = T)

        set.seed(619)
        ## spliting samples
        list_ <- SplitData(RatioMat, trainginSize = 1, do.scale = F)

        ## lasso-cox
        df <- FitLassoCoxBoost(list_, times = 1000)
        p <- AnalysisFeatureAfterBoost(df)
        pdf(paste0(savePath, "lassocox_boostrap features (PD1+ Treg)(PID)(ALL).pdf"), height = 6, width = 10)
        print(p)
        dev.off()

        ## ROI level
        ## Using Cell subpopulation fraction ratio to regression
        countdf <- GetAbundance(sce, "SubType", is.fraction = T, is.reuturnMeans = F)
        clinical <- countdf[, c(26, 25)]
        countdf <- countdf[, -c(23:ncol(countdf))]

        countdf <- cbind(countdf, clinical)
        RatioMat <- TranferAbun2Ratio(countdf, CombFeatureTypes = CombFeatureTypes)
        ## Uni-variable cox to fileter feature
        RatioMat <- UniCoxFilter(RatioMat, do.padjust = F, do.scale = T)

        set.seed(619)
        ## spliting samples
        list_ <- SplitData(RatioMat, trainginSize = 1, do.scale = F)

        ## lasso-cox
        df <- FitLassoCoxBoost(list_, times = 1000)
        p <- AnalysisFeatureAfterBoost(df)
        pdf(paste0(savePath, "lassocox_boostrap features (PD1+ Treg)(ROI)(ALL).pdf"), height = 6, width = 10)
        print(p)
        dev.off()
    }

    ### paralle
    if (F) {
        library(parallel)

        perFunc <- function(list_) {
            source("/home/lyx/project/IMC/structural_analysis_functions.r")
            source("./signature_functions.r")
            sce <- list_[[1]]
            tissue <- list_[[2]]

            cat("Performing phenograph clustering with phenok = ", phenok, "\n")
            sce <- sce[, sce$Tissue == tissue]

            savePath <- paste0("/mnt/data/lyx/IMC/analysis/signature/", tissue, "/")
            if (!dir.exists(savePath)) {
                dir.create(savePath, recursive = T)
            }

            ## Using Cell subpopulation fraction ratio to regression
            clinicalFeatures <- c("RFS_time", "RFS_status")
            countdf <- GetAbundance(sce, "SubType", is.fraction = T, is.reuturnMeans = T)
            clinical <- countdf[, c(25, 24)]
            countdf <- countdf[, -c(21:ncol(countdf))]

            countdf <- cbind(countdf, clinical)

            CombFeatureTypes <- colnames(countdf)[1:20]

            RatioMat <- TranferAbun2Ratio(countdf, CombFeatureTypes = CombFeatureTypes)
            ## Uni-variable cox to fileter feature
            RatioMat <- UniCoxFilter(RatioMat, do.padjust = F, do.scale = T)

            set.seed(619)
            ## spliting samples
            list_ <- SplitData(RatioMat, trainginSize = 1, do.scale = F)

            ## lasso-cox
            df <- FitLassoCoxBoost(list_, times = 1000)
            p <- AnalysisFeatureAfterBoost(df)
            pdf(paste0(savePath, "lassocox_boostrap features (PID)(ALL).pdf"), height = 6, width = 10)
            print(p)
            dev.off()

            ## ROI level
            ## Using Cell subpopulation fraction ratio to regression
            countdf <- GetAbundance(sce, "SubType", is.fraction = T, is.reuturnMeans = F)
            clinical <- countdf[, c(25, 24)]
            countdf <- countdf[, -c(21:ncol(countdf))]

            countdf <- cbind(countdf, clinical)
            RatioMat <- TranferAbun2Ratio(countdf, CombFeatureTypes = CombFeatureTypes)
            ## Uni-variable cox to fileter feature
            RatioMat <- UniCoxFilter(RatioMat, do.padjust = F, do.scale = T)

            set.seed(619)
            ## spliting samples
            list_ <- SplitData(RatioMat, trainginSize = 1, do.scale = F)

            ## lasso-cox
            df <- FitLassoCoxBoost(list_, times = 1000)
            p <- AnalysisFeatureAfterBoost(df)
            pdf(paste0(savePath, "lassocox_boostrap features (ROI)(ALL).pdf"), height = 6, width = 10)
            print(p)
            dev.off()
            return(NULL)
        }

        ## form targets to multiple precess
        cat("Form the multi-process targets", "\n")
        sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
        tissue <- c("CT", "TAT")

        targets <- list()
        z <- 1
        for (i in 1:length(tissue)) {
            targets[[z]] <- list(sce, tissue[z])
            z <- z + 1
        }
        cat("Form the multi-process targets done.", "\n")
        cl <- makeCluster(12)
        cat("Start multi-process.", "\n")
        results <- parLapply(cl, targets, perFunc)
        stopCluster(cl)
    }
}

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce(Mono_CD11c).rds")

savePath <- paste0("/mnt/data/lyx/IMC/analysis/signature/")
if (!dir.exists(savePath)) {
    dir.create(savePath, recursive = T)
}

## Using ratio lasso-cox model to predict the
clinicalFeatures <- c("RFS_time", "RFS_status", "KRAS_mutation", "fong_score", "Age", "TBS", "CRLM_number", "CRLM_size", "CEA", "CA199", "T_stage")
countdf <- GetAbundance(sceobj = sce, countcol = "SubType", clinicalFeatures = clinicalFeatures, is.fraction = T, is.reuturnMeans = T)

colnames(countdf)
clinical <- countdf[, c(24:ncol(countdf))]
countdf <- countdf[, -c(22:ncol(countdf))]

FeatureType1 <- c("B", "CD4T", "CD8T")
FeatureType2 <- c("Macro_CD163", "PDL1+ Mono_CD11c", "Mono_CD11c", "Treg")

RatioMat <- TranferAbun2Ratio(countdf, FeatureType1 = FeatureType1, FeatureType2 = FeatureType2)

RatioMat$ID <- rownames(RatioMat)
clinical$ID <- rownames(clinical)

RatioMat2 <- dplyr::left_join(RatioMat, clinical[, c("ID", "RFS_time", "RFS_status")], by = "ID")
RatioMat2 <- RatioMat2[, -which(colnames(RatioMat2) == "ID")]
rownames(RatioMat2) <- rownames(RatioMat)

list_ <- SplitData(RatioMat2, trainginSize = 1, do.scale = T)
x <- list_[["x_train"]]
y <- list_[["y_train"]]

x <- as.matrix(x)
y <- as.matrix(y)
colnames(y) <- c("time", "status")

set.seed(619)

for (i in 1:1000) {
    cindex_lamda_set <- cv.glmnet(x, y, family = "cox", type.measure = "C", alpha = 1, nfolds = 4, maxit = 1e6)
    feature_wights <- coef(cindex_lamda_set, s = "lambda.min")
    feature_wights <- as.data.frame(as.matrix(feature_wights))
    weights <- as.numeric(feature_wights[, 1])
    weightname <- rownames(feature_wights)

    model <- list()
    model[["model"]] <- cindex_lamda_set
    model[["weights"]] <- feature_wights
    model[["cindex"]] <- cindex_lamda_set$cvm[cindex_lamda_set$index["min", ]]

    model[["active_features"]] <- weightname[weights != 0]
    if (length(model[["active_features"]]) == 0) {
        next
    }

    x_ <- x[, model[["active_features"]]]
    coxph_model <- coxph(Surv(time = as.double(y[, 1]), event = as.double(y[, 2])) ~ x_)
    model[["coxph_model"]] <- coxph_model

    # AnalysisFeature(model, savePath)

    pred_train <- Inference(list_[["x_train"]], list_[["y_train"]], model)
    compareDF <- cbind("riskScore" = pred_train[[1]], "riskLabel" = pred_train[[2]], clinical)

    mul_cox <- coxph(Surv(RFS_time, RFS_status) ~ riskScore + riskLabel + fong_score + Age + TBS + CRLM_number + CRLM_size + CEA + CA199 + T_stage, data = compareDF)
    mul_cox1 <- summary(mul_cox)
    multi1 <- as.data.frame(round(mul_cox1$conf.int[, c(1, 3, 4)], 2))
    multi2 <- ShowRegTable(mul_cox,
        exp = TRUE, digits = 2, pDigits = 3,
        printToggle = TRUE, quote = FALSE, ciFun = confint
    )

    result <- cbind(multi1, multi2)

    if (as.numeric(result[2, 5]) <= 0.05) {
        saveRDS(model, "/mnt/data/lyx/IMC/analysis/signature/model.rds")
        break
    }
    if (i %% 50 == 0) {
        print(i)
    }
}
