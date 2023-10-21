# These functions for building signature
library(glmnet)
library(survival)
library(survminer)
library(ggthemes)
library(tableone)
library(SingleCellExperiment)
library(boot)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(forestplot)
library(Hmisc)

## Transfer Cell subpopulation abundance into ratio
TranferAbun2Ratio <- function(countdf, FeatureType1 = NULL, FeatureType2 = NULL, CombFeatureTypes = NULL) {
    ## Number of feature
    if (is.null(CombFeatureTypes)) {
        Features <- rbind(expand.grid(FeatureType1, FeatureType2), expand.grid(FeatureType2, FeatureType1))
    }
    if (!is.null(CombFeatureTypes)) {
        Features <- expand.grid(CombFeatureTypes, CombFeatureTypes)
        idx <- c()
        for (i in 1:nrow(Features)) {
            if (Features[i, 1] == Features[i, 2]) {
                idx <- c(idx, i)
            }
        }
        Features <- Features[!(rownames(Features) %in% idx), ]
    }
    numFeature <- nrow(Features)

    ## calcualte ratio
    RatioMat <- matrix(data = NA, nrow = nrow(countdf), ncol = numFeature)
    featureName <- c()

    z <- 1
    for (i in 1:nrow(Features)) {
        c1 <- as.character(Features[i, 1])
        c2 <- as.character(Features[i, 2])
        Vec1 <- countdf[, c1]
        Vec2 <- countdf[, c2]

        featureName <- c(featureName, paste0(c1, ":", c2))
        RatioTemp <- Vec1 / Vec2
        max_ <- max(RatioTemp[!is.infinite(RatioTemp)])
        RatioTemp[is.infinite(RatioTemp)] <- max_
        RatioMat[, z] <- RatioTemp
        z <- z + 1
    }

    RatioMat <- ifelse(is.na(RatioMat), 0, RatioMat)

    colnames(RatioMat) <- featureName
    rownames(RatioMat) <- rownames(countdf)

    RatioMat <- as.data.frame(RatioMat)

    return(RatioMat)
}

## Using univariable cox to filter features
UniCoxFilter <- function(RatioMat, do.padjust = F, do.scale = T) {
    numFeature <- length(RatioMat) - 2

    x <- RatioMat[, 1:numFeature]
    y <- RatioMat[, (numFeature + 1):ncol(RatioMat)]

    if (do.scale) {
        x <- apply(x, MARGIN = 2, function(a) {
            return(scale(a))
        })
    }

    ## Iteractive COX
    result <- matrix(data = NA, nrow = 0, ncol = 6)
    result <- as.data.frame(result)

    colnames(x) <- sapply(1:ncol(x), function(a) {
        return(paste0("feature_", a))
    })
    df <- cbind(x, y)

    univ_formulas <- sapply(
        colnames(x),
        function(a) {
            as.formula(paste0("Surv(RFS_time, RFS_status)~", a))
        }
    )

    univ_models <- lapply(univ_formulas, function(a) {
        coxph(a, data = df)
    })

    univ_results <- lapply(univ_models, function(a) {
        mul_cox1_result <- summary(a)
        multi1 <- round(mul_cox1_result$conf.int[, c(1, 3, 4)], 4)
        multi1 <- as.data.frame(t(multi1))

        multi2 <- ShowRegTable(
            a,
            exp = TRUE,
            digits = 2, pDigits = 3,
            printToggle = TRUE, quote = FALSE, ciFun = confint
        )

        result <- cbind(multi1, multi2)
        result <- cbind(Features = rownames(result), result)

        return(result)
    })

    for (a in univ_results) {
        result <- rbind(result, a)
    }
    pvalues <- result$p
    pvalues <- ifelse(pvalues == "<0.001", 0.001, pvalues)
    if (do.padjust) {
        qvalues <- p.adjust(as.numeric(pvalues), method = "BH")
        idx <- (qvalues <= 0.05)
    } else {
        idx <- (as.numeric(pvalues) <= 0.05)
    }

    returnMat <- cbind(RatioMat[, idx], y)
    return(returnMat)
}

## Spliting data
SplitData <- function(RatioMat, trainginSize = 0.8, selectingFeature = NULL, do.scale = T) {
    if (!is.null(selectingFeature)) {
        clinical <- RatioMat[, ((ncol(RatioMat) - 1):ncol(RatioMat))]
        RatioMat <- cbind(RatioMat[, colnames(RatioMat) %in% selectingFeature], clinical)
    }

    n <- nrow(RatioMat)
    n_feature <- ncol(RatioMat) - 2
    n_train <- round(n * trainginSize, 0)

    id_train <- sample(1:n, n_train)
    id_test <- (1:n)[-id_train]

    ## scaling
    if (do.scale) {
        RatioMat[, 1:n_feature] <- apply(RatioMat[, 1:n_feature], MARGIN = 2, function(x) {
            return(scale(x))
        })
    }

    x_train <- RatioMat[id_train, (1:n_feature)]
    y_train <- RatioMat[id_train, ((ncol(RatioMat) - 1):ncol(RatioMat))]

    x_test <- RatioMat[id_test, (1:n_feature)]
    y_test <- RatioMat[id_test, ((ncol(RatioMat) - 1):ncol(RatioMat))]

    list_ <- list()
    list_[["NumFeature"]] <- n_feature
    list_[["TrainSize"]] <- length(id_train)
    list_[["TestSize"]] <- length(id_test)
    list_[["x_train"]] <- x_train
    list_[["y_train"]] <- y_train
    list_[["x_test"]] <- x_test
    list_[["y_test"]] <- y_test

    return(list_)
}

## Using lasso-cox to fit the model
FitLassoCox <- function(list_, savePath) {
    x <- list_[["x_train"]]
    y <- list_[["y_train"]]

    x <- as.matrix(x)
    y <- Surv(y[, 1], y[, 2])

    ## max-normalization
    ## maxium C-index
    cindex_lamda_set <- cv.glmnet(x, y, family = "cox", type.measure = "C", alpha = 1, nfolds = 4, maxit = 1e6)

    pdf(file = paste0(savePath, "lasso-cox model cindex lamda_set.pdf"), width = 8, height = 6)
    print(plot(cindex_lamda_set))
    dev.off()

    best_lambda <- cindex_lamda_set$lambda.min

    feature_wights <- coef(cindex_lamda_set, s = "lambda.min")
    feature_wights <- as.data.frame(as.matrix(feature_wights))
    weights <- as.numeric(feature_wights[, 1])
    weightname <- rownames(feature_wights)

    model <- list()
    model[["model"]] <- cindex_lamda_set
    model[["weights"]] <- feature_wights

    # Compute predicted risk scores using the best lambda value
    risk_scores <- predict(cindex_lamda_set, newx = x, s = best_lambda, type = "link")
    c_index <- 1 - (rcorr.cens(risk_scores, y)[[1]])
    model[["cindex"]] <- c_index

    model[["active_features"]] <- weightname[weights != 0]
    if (length(model[["active_features"]]) == 0) {
        next
    }

    x_ <- x[, model[["active_features"]]]
    coxph_model <- coxph(Surv(time = as.double(y[, 1]), event = as.double(y[, 2])) ~ x_)
    model[["coxph_model"]] <- coxph_model

    return(model)
}

## Boostrapping with lasso cox
lassocox_boost <- function(data, indices) {
    # Generate a bootstrap sample using the provided indices
    sample_data <- data[indices, ]

    # Prepare the response variable (Surv object) and covariates
    y <- Surv(sample_data$time, sample_data$status)
    x <- as.matrix(sample_data[, -c(1, 2)])

    # Fit the Lasso-Cox model using cross-validation to select the optimal lambda
    cv_fit <- cv.glmnet(x, y, family = "cox", nfolds = 3, type.measure = "C")

    # Extract the coefficients for the optimal lambda
    coef <- as.vector(coef(cv_fit, s = "lambda.min"))

    return(coef)
}

FitLassoCoxBoost <- function(list_, times = 1000) {
    ## Traning data
    data <- cbind(list_[["y_train"]], list_[["x_train"]])
    colnames(data) <- c("time", "status", colnames(data)[3:ncol(data)])
    data <- as.data.frame(data)
    ## Boostrap
    bootstrap_results <- boot(data, lassocox_boost, R = times)

    ## Features
    features <- colnames(list_[["x_train"]])
    # Calculate the selection frequencies
    non_zero_coefs <- apply(bootstrap_results$t, 2, function(x) x != 0)
    selected_features <- colSums(non_zero_coefs) / nrow(non_zero_coefs)
    # Calculate the average coefficient values
    avg_coefs <- colMeans(bootstrap_results$t)

    retudf <- matrix(data = NA, nrow = length(features), ncol = 3)
    retudf <- as.data.frame(retudf)
    retudf[, 1] <- features
    retudf[, 2] <- avg_coefs
    retudf[, 3] <- selected_features
    colnames(retudf) <- c("Feature", "Avg_coef", "Freaquency")
    return(retudf)
}

## Visualize Feature after boostrapping selection
AnalysisFeatureAfterBoost <- function(df) {
    df <- df[order(df$Avg_coef), ]
    df$Feature <- factor(df$Feature, levels = df$Feature)

    p <- ggplot(df, aes(x = Feature, y = Avg_coef, fill = Avg_coef > 0, alpha = Freaquency)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = c("steelblue", "darkorange"), labels = c("Protect", "Risk"), name = "RFS-associated coef") +
        labs(x = "Marker", y = "Average Coefficient") +
        theme_minimal() +
        theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 12, face = "bold"),
            axis.text = element_text(size = 10),
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10),
            legend.position = "bottom"
        )
    return(p)
}

## analysis the features
AnalysisFeature <- function(model, savePath) {
    weights <- model[["weights"]]
    features <- model[["active_features"]]

    weights <- weights[which(weights != 0), ]
    coxmodel <- model[["coxph_model"]]

    coxmodel <- summary(coxmodel)

    model_coefficients <- as.data.frame(matrix(data = NA, nrow = length(weights), ncol = 4))
    model_coefficients[, 1] <- features
    model_coefficients[, 2] <- weights
    model_coefficients[, 3] <- as.numeric(coxmodel$coef[, 1])
    model_coefficients[, 4] <- as.numeric(coxmodel$coef[, 5])

    colnames(model_coefficients) <- c("Feature", "Coefficient", "lnHR", "pvalue")
    model_coefficients$lnHR <- ifelse(model_coefficients$lnHR > 10, 10, model_coefficients$lnHR)
    model_coefficients$lnHR <- ifelse(model_coefficients$lnHR < (-10), -10, model_coefficients$lnHR)


    p <- ggplot(data = model_coefficients, aes(x = reorder(Feature, -Coefficient), y = Coefficient, fill = Coefficient > 0, alpha = pvalue)) +
        geom_bar(stat = "identity", width = 0.7, color = "black") +
        scale_fill_manual(values = c("darkred", "darkblue"), guide = FALSE) +
        theme_minimal() +
        theme(
            text = element_text(size = 14, family = "Helvetica"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 8, margin = margin(r = 3, l = 3)),
            axis.title.x = element_text(size = 16, margin = margin(t = 10)),
            axis.title.y = element_text(size = 16, margin = margin(r = 10)),
            plot.title = element_text(size = 20, margin = margin(b = 20)),
            plot.subtitle = element_text(size = 14, margin = margin(b = 10)),
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA),
            panel.grid.major = element_line(color = "grey80"),
            panel.grid.minor = element_blank(),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.key = element_rect(fill = "transparent", color = NA)
        ) +
        labs(
            x = "Features",
            y = "Coefficients"
        ) +
        coord_flip()

    # Add HR as point size
    p <- p + geom_point(aes(x = Feature, y = Coefficient, size = abs(lnHR)), shape = 21, color = "black")

    # Adjust point size scale and add alpha scale
    p <- p + scale_size_continuous(range = c(2, 8), name = "lnHR") +
        scale_alpha_continuous(range = c(1, 0.1), name = "pvalue", limits = c(0, 1))

    pdf(paste0(savePath, "Diverging Barplot of Feature Coefficients.pdf"), , height = 6, width = 8)
    print(p)
    dev.off()

    return(NULL)
}

## Inference
Inference <- function(x, y, model) {
    ## Extract data and model
    x <- as.matrix(x)

    weights <- model[["weights"]]
    weights <- weights[which(weights != 0), ]
    weights <- as.matrix(weights, ncol = 1)
    features <- model[["active_features"]]

    ## Inference
    x <- x[, features]

    riskScore <- x %*% weights

    ## Assign label
    df <- cbind(y, riskScore)
    colnames(df) <- c("RFS_time", "RFS_status", "RiskScore")
    if (!is.data.frame(df)) {
        df <- as.data.frame(df)
    }
    cutpoint <- surv_cutpoint(data = df, time = "RFS_time", event = "RFS_status", variables = "RiskScore")
    cutpoint <- summary(cutpoint)$cutpoint
    label <- ifelse(df[, 3] >= cutpoint, "high", "low")

    list_ <- list()
    list_[["riskScore"]] <- as.numeric(riskScore)
    list_[["riskLabel"]] <- as.character(label)

    return(list_)
}

# Define a function for plotting the Kaplan-Meier curve
KMplot <- function(clinical, label, set, savePath, plot = T) {
    # Convert the first two columns of the clinical data to numeric
    clinical[, 1] <- as.numeric(clinical[, 1])
    clinical[, 2] <- as.numeric(clinical[, 2])

    # Combine clinical data and risk labels into a single data frame
    data <- cbind(clinical, score = label)
    if (!is.data.frame(data)) {
        data <- as.data.frame(data)
    }

    # Perform survival analysis using the Surv() function and fit the data using surv_fit()
    fit <- surv_fit(Surv(clinical[, 1], clinical[, 2]) ~ label, data = data)

    # Perform log-rank test using survdiff()
    logrank_pvalue <- survdiff(Surv(clinical[, 1], clinical[, 2]) ~ label, data = data)
    logrank_pvalue <- 1 - pchisq(logrank_pvalue$chisq, length(logrank_pvalue$n) - 1)

    # Calculate the concordance index (C-index)
    c_index <- 1 - (rcorr.cens(label, Surv(clinical[, 1], clinical[, 2]))[[1]])
    c_index <- round(digits = 3, c_index)

    if (plot) {
        # Create a Kaplan-Meier plot with risk table and confidence intervals
        p <- ggsurvplot(fit,
            data = data,
            linetype = c("solid", "solid"),
            surv.median.line = "hv",
            surv.scale = "percent",
            pval = T,
            risk.table = T,
            conf.int = T,
            conf.int.alpha = 0.1,
            conf.int.style = "ribbon",
            risk.table.y.text = T,
            palette = c("#CC3300", "#3300CC"),
            xlab = "Relapse-Free survival time(month)"
        )

        # Customize the plot appearance
        p$plot <- p$plot +
            theme_minimal() +
            theme(
                plot.title = element_text(face = "bold", size = 20),
                plot.subtitle = element_text(size = 16),
                plot.caption = element_text(size = 12),
                axis.title.x = element_text(face = "bold", size = 16),
                axis.title.y = element_text(face = "bold", size = 16),
                axis.text.x = element_text(size = 14),
                axis.text.y = element_text(size = 14),
                legend.text = element_text(size = 14),
                legend.title = element_text(face = "bold", size = 16),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_line(color = "grey90")
            ) +
            labs(caption = paste0("C-index: ", c_index))

        # Save the plot as a PDF file
        pdf(paste0(savePath, "KM of RFS in ", set, ".pdf"), height = 6, width = 8)

        # Print the plot
        print(p)
        dev.off()
    }

    # Return cindex and p value
    return_list <- list()
    return_list[["logrank_p"]] <- logrank_pvalue
    return_list[["c-index"]] <- c_index

    return(return_list)
}

## Get mean expression profile
GetMeanExpressionProfile <- function(sce, LevelCol, markers, clinicalGroupCol) {
    ids <- names(table(colData(sce)[, LevelCol]))

    df <- matrix(data = NA, nrow = length(ids), ncol = length(markers))
    clinicalMat <- matrix(data = NA, nrow = 0, ncol = length(clinicalGroupCol))
    for (i in 1:length(ids)) {
        idx <- ids[i]
        ## Subset the sce object
        sceTemp <- sce[, colData(sce)[, LevelCol] %in% idx]
        numCell <- ncol(sceTemp)
        ## Get the expression profile
        expTemp <- assay(sceTemp)
        expTemp <- expTemp[match(markers, rownames(expTemp)), ]
        df[i, ] <- apply(expTemp, MARGIN = 1, FUN = "mean")
        clinicalTemp <- colData(sceTemp)[1, match(clinicalGroupCol, colnames(colData(sceTemp)))]
        clinicalMat <- rbind(clinicalMat, as.matrix(clinicalTemp))
    }

    ## rename
    rownames(df) <- ids
    colnames(df) <- markers
    colnames(clinicalMat) <- clinicalGroupCol
    rownames(clinicalMat) <- ids

    ## return list
    list_ <- list()
    list_[["df"]] <- df
    list_[["clinicalGroup"]] <- clinicalMat
    return(list_)
}

## Boxplot for cell subpopulations under certain clinical groups
RatioBoxPlot <- function(RatioMat2, expCol, clinicalGroupCol, ClinicalGroupName, do.scale = F, FilterNotSig = F, p_threshold = 0.05) {
    metadf <- RatioMat2[, -expCol]
    metadf <- metadf[, !colnames(metadf) %in% "ID"]
    expdf <- as.matrix(RatioMat2[, expCol])
    if (do.scale) {
        expdf <- apply(expdf, MARGIN = 2, function(x) {
            return(scale(x, center = F, scale = T))
        })
    }
    rownames(expdf) <- rownames(RatioMat2)

    AbunBoxDF <- as.data.frame(matrix(data = NA, ncol = (3 + ncol(metadf)), nrow = nrow(expdf) * ncol(expdf)))
    AbunBoxDF[, 1] <- as.numeric(as.matrix(expdf))
    AbunBoxDF[, 2] <- rep(colnames(expdf), each = nrow(RatioMat2))
    AbunBoxDF[, 3] <- rep(rownames(RatioMat2), times = ncol(expdf))
    for (i in 4:ncol(AbunBoxDF)) {
        AbunBoxDF[, i] <- rep(metadf[, i - 3], times = ncol(expdf))
    }

    colnames(AbunBoxDF) <- c("Ratio", "Pair", "ID", colnames(metadf))

    AbunBoxDF[, clinicalGroupCol] <- as.factor(ifelse(AbunBoxDF[, clinicalGroupCol] == 1, ClinicalGroupName[1], ClinicalGroupName[2]))
    colnames(AbunBoxDF)[which(colnames(AbunBoxDF) == clinicalGroupCol)] <- "ClinicalGroup"

    if (FilterNotSig) {
        AllPair <- names(table(AbunBoxDF[, "Pair"]))
        remainType <- c()
        for (pair_ in AllPair) {
            AbunBoxDF_ <- AbunBoxDF[AbunBoxDF$Pair %in% pair_, ]
            g1 <- AbunBoxDF_[AbunBoxDF_$ClinicalGroup %in% ClinicalGroupName[1], "Ratio"]
            g2 <- AbunBoxDF_[AbunBoxDF_$ClinicalGroup %in% ClinicalGroupName[2], "Ratio"]

            p.sig <- t.test(g1, g2)$p.value
            if (p.sig <= p_threshold) {
                remainType <- c(remainType, pair_)
            }
        }
        AbunBoxDF <- AbunBoxDF[AbunBoxDF$Pair %in% remainType, ]
    }

    p <- ggplot(AbunBoxDF, aes(x = Pair, y = Ratio, fill = ClinicalGroup)) +
        geom_boxplot(alpha = 0.7) +
        scale_y_continuous(name = "Cell Abundance") +
        scale_x_discrete(name = "Cell Population") +
        theme_bw() +
        theme(
            plot.title = element_text(size = 14, face = "bold"),
            text = element_text(size = 12),
            axis.title = element_text(face = "bold"),
            axis.text.x = element_text(size = 11, angle = 90)
        ) +
        scale_fill_lancet() +
        stat_compare_means(aes(group = ClinicalGroup), label.y = max(expdf), method = "t.test", label = "p.signif")

    return(p)
}

## Ratio correlation
RatioCorPlot <- function(RatioMat2, expCol, clinicalGroupCol, ClinicalGroupName, plotType, do.scale = F) {
    metadf <- RatioMat2[, -expCol]
    metadf <- metadf[, !colnames(metadf) %in% "ID"]
    expdf <- as.matrix(RatioMat2[, expCol])
    if (do.scale) {
        expdf <- apply(expdf, MARGIN = 2, function(x) {
            return(scale(x, center = F, scale = T))
        })
    }
    rownames(expdf) <- rownames(RatioMat2)

    AbunBoxDF <- as.data.frame(matrix(data = NA, ncol = (3 + ncol(metadf)), nrow = nrow(expdf) * ncol(expdf)))
    AbunBoxDF[, 1] <- as.numeric(as.matrix(expdf))
    AbunBoxDF[, 2] <- rep(colnames(expdf), each = nrow(RatioMat2))
    AbunBoxDF[, 3] <- rep(rownames(RatioMat2), times = ncol(expdf))
    for (i in 4:ncol(AbunBoxDF)) {
        AbunBoxDF[, i] <- rep(metadf[, i - 3], times = ncol(expdf))
    }

    colnames(AbunBoxDF) <- c("Ratio", "Pair", "ID", colnames(metadf))

    AbunBoxDF[, clinicalGroupCol] <- as.factor(ifelse(AbunBoxDF[, clinicalGroupCol] == 1, ClinicalGroupName[1], ClinicalGroupName[2]))
    colnames(AbunBoxDF)[which(colnames(AbunBoxDF) == clinicalGroupCol)] <- "ClinicalGroup"

    AbunBoxDF <- AbunBoxDF[AbunBoxDF$Pair %in% plotType, ]

    p <- ggplot(AbunBoxDF, aes(x = RFS_time, y = Ratio)) +
        geom_point(aes(color = ClinicalGroup)) +
        geom_smooth(method = "lm", formula = y ~ x, se = T) +
        labs(x = "RFS_time", y = "Ratio") +
        theme_minimal() +
        scale_color_manual(
            values = c("R" = "#BB0021FF", "NR" = "#3B4992FF"),
            name = "Clinical Group",
            labels = c("R" = "R", "NR" = "NR")
        ) +
        stat_cor(label.x = 0.1, color = "black") +
        facet_wrap(~Pair)


    return(p)
}

## Compare signatures via multi-variable cox
MultipleUniCOX <- function(df, UniCOX = TRUE) {
    result <- matrix(data = NA, nrow = 0, ncol = 6)
    result <- as.data.frame(result)

    colnames(df) <- sapply(colnames(df), function(x) {
        sub(pattern = "\\+", replacement = "", x)
    })
    colnames(df) <- sapply(colnames(df), function(x) {
        sub(pattern = " ", replacement = "_", x)
    })

    features <- colnames(df)[1:(ncol(df) - 2)]
    cat("The features in multi-cox are: ", features, "\n")

    # features <- features[-(match(c("RFS_time", "RFS_status"), features))]

    if (UniCOX) {
        univ_formulas <- sapply(
            features,
            function(x) {
                as.formula(paste0("Surv(RFS_time, RFS_status)~", x))
            }
        )

        univ_models <- lapply(univ_formulas, function(x) {
            coxph(x, data = df)
        })

        univ_results <- lapply(univ_models, function(x) {
            ftest <- cox.zph(x)
            ftest
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
    } else {
        formula_ <- paste0("Surv(RFS_time, RFS_status)~", features[1])
        for (i in 2:length(features)) {
            formula_ <- paste0(formula_, "+", features[i])
        }
        mul_cox <- coxph(as.formula(formula_), data = df)
        mul_cox1 <- summary(mul_cox)

        multi1 <- as.data.frame(round(mul_cox1$conf.int[, c(1, 3, 4)], 2))
        multi2 <- ShowRegTable(mul_cox,
            exp = TRUE,
            digits = 2,
            pDigits = 3,
            printToggle = TRUE,
            quote = FALSE,
            ciFun = confint
        )
        result <- cbind(Features = rownames(multi2), multi1, multi2)
        return(result)
    }
}

CompareSig <- function(mul_cox) {
    mul_cox1 <- summary(mul_cox)
    colnames(mul_cox1$conf.int)
    multi1 <- as.data.frame(round(mul_cox1$conf.int[, c(1, 3, 4)], 2))
    multi2 <- ShowRegTable(mul_cox,
        exp = TRUE,
        digits = 2,
        pDigits = 3,
        printToggle = TRUE,
        quote = FALSE,
        ciFun = confint
    )
    result <- cbind(multi1, multi2)
    result <- tibble::rownames_to_column(result, var = "Characteristics")
}

## Visualize the c-index of each features
LollipopForMultiFeatures <- function(pcoxDF) {
    pcoxDF[, 1] <- as.character(pcoxDF[, 1])
    pcoxDF[, 2] <- as.numeric(pcoxDF[, 2])
    pcoxDF[, 3] <- as.numeric(pcoxDF[, 3])

    p <- ggdotchart(pcoxDF,
        x = "Features", y = "Cindex",
        color = "#E377C2FF",
        sorting = "ascending", add = "segments", rotate = F, dot.size = 12,
        label = round(as.numeric(pcoxDF[, "Cindex"]), digits = 3),
        font.label = list(color = "black", size = 10, vjust = 0.5),
        ggtheme = theme_classic()
    ) +
        ylab("C-index") + xlab("Clinical Features") +
        coord_flip()

    return(p)
}

## Get the cell subpopulation abundance fraction from certain tissue
GetFractionFromTissue <- function(df, tissue, subtype) {
    df2 <- subset(df, Tissue == tissue)
    df2 <- df2[, match(c("PID", subtype), colnames(df2))]
    colnames(df2)[2] <- paste0(tissue, "_", subtype)

    return(df2)
}

## forest plot to visualize the cox results
ForestPlot <- function(multicox_result, savePath) {
    ## Definate space
    ins <- function(x) {
        c(as.character(x), rep(NA, ncol(multicox_result) - 1))
    }
    numfeatures <- nrow(multicox_result)

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

    ## forestplot
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
        col = fpColors(
            box = "#021eaa",
            lines = "#021eaa",
            zero = "black"
        )
    )

    pdf(savePath, width = 8, height = 6)
    print(p)
    dev.off()
}
