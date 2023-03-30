# These functions for building signature
library(glmnet)
library(survival)
library(survminer)
library(ggthemes)
library(tableone)

## Transfer Cell subpopulation abundance into ratio
TranferAbun2Ratio <- function(countdf, featureTypes) {
    ## Number of feature
    numType <- length(featureTypes)
    numFeature <- numType * (numType - 1)

    ## calcualte ratio
    RatioMat <- matrix(data = NA, nrow = nrow(countdf), ncol = numFeature)
    featureName <- c()

    z <- 1
    for (i in 1:numType) {
        c1 <- featureTypes[i]
        for (j in 1:numType) {
            if (i == j) {
                next
            }

            c2 <- featureTypes[j]

            Vec1 <- countdf[, c1]
            Vec2 <- countdf[, c2]

            featureName <- c(featureName, paste0(c1, "-", c2))
            RatioMat[, z] <- Vec1 / Vec2
            z <- z + 1
        }
    }

    RatioMat <- ifelse(is.na(RatioMat), 0, RatioMat)
    RatioMat <- ifelse(is.infinite(RatioMat), 0, RatioMat)

    colnames(RatioMat) <- featureName
    rownames(RatioMat) <- rownames(countdf)

    RatioMat <- as.data.frame(RatioMat)
    RatioMat$RFS_time <- countdf$RFS_time
    RatioMat$RFS_status <- countdf$RFS_status

    return(RatioMat)
}

## Using univariable cox to filter features
UniCoxFilter <- function(RatioMat, do.padjust=F) {
    numFeature <- length(RatioMat) - 2

    x <- RatioMat[, 1:numFeature]
    y <- RatioMat[, (numFeature + 1):ncol(RatioMat)]

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
    if(do.padjust){
        qvalues <- p.adjust(as.numeric(pvalues),method="BH")
        idx <- (qvalues <= 0.05)
    }
    else{
        idx <- (as.numeric(pvalues) <= 0.05)
    }
    
    returnMat <- cbind(RatioMat[, idx],y)
    return(returnMat)
}

## Spliting data
SplitData <- function(RatioMat, trainginSize = 0.8, do.scale = T) {
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

    y <- as.matrix(y)
    colnames(y) <- c("time", "status")

    ## max-normalization
    ## maxium C-index
    cindex_lamda_set <- cv.glmnet(x, y, family = "cox", type.measure = "C", alpha = 1, nfolds = 10, maxit = 100000)

    pdf(file = paste0(savePath, "lasso-cox model cindex lamda_set.pdf"), width = 8, height = 6)
    print(plot(cindex_lamda_set))
    dev.off()

    feature_wights <- coef(cindex_lamda_set, s = "lambda.min")

    model <- list()
    model[["model"]] <- cindex_lamda_set
    model[["weights"]] <- feature_wights
    model[["cindex"]] <- cindex_lamda_set$cvm[cindex_lamda_set$index["min", ]]
    model[["active_features"]] <- names(feature_wights[which(feature_wights != 0), ])

    x_ <- x[, model[["active_features"]]]
    coxph_model <- coxph(Surv(time = as.double(y[, 1]), event = as.double(y[, 2])) ~ x_)
    model[["coxph_model"]] <- coxph_model

    return(model)
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
    model_coefficients[, 3] <- as.numeric(coxmodel$coef[, 2])
    model_coefficients[, 4] <- as.numeric(coxmodel$coef[, 5])

    colnames(model_coefficients) <- c("Feature", "Coefficient", "HR", "pvalue")


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
    p <- p + geom_point(aes(x = Feature, y = Coefficient, size = HR), shape = 21, color = "black")

    # Adjust point size scale and add alpha scale
    p <- p + scale_size_continuous(range = c(2, 8), name = "HR") +
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
    cutpoint <- surv_cutpoint(data = df, time = "RFS_time", event = "RFS_status", variables = "RiskScore")
    cutpoint <- summary(cutpoint)$cutpoint
    label <- ifelse(df[, 3] >= cutpoint, "high", "low")

    list_ <- list()
    list_[["riskScore"]] <- as.numeric(riskScore)
    list_[["riskLabel"]] <- as.character(label)

    return(list_)
}

# Define a function for plotting the Kaplan-Meier curve
KMplot <- function(clinical, pred_list, set, savePath) {
    # Extract risk labels from the prediction list
    label <- pred_list[["riskLabel"]]

    # Convert the first two columns of the clinical data to numeric
    clinical[, 1] <- as.numeric(clinical[, 1])
    clinical[, 2] <- as.numeric(clinical[, 2])

    # Combine clinical data and risk labels into a single data frame
    data <- cbind(clinical, score = label)

    # Perform survival analysis using the Surv() function and fit the data using surv_fit()
    fit <- surv_fit(Surv(clinical[, 1], clinical[, 2]) ~ label, data = data)

    # Perform log-rank test using survdiff()
    logrank_pvalue <- survdiff(Surv(clinical[, 1], clinical[, 2]) ~ label, data = data)
    logrank_pvalue <- 1 - pchisq(logrank_pvalue$chisq, length(logrank_pvalue$n) - 1)

    # Fit the data using the Cox proportional hazards model (coxph())
    fit2 <- coxph(Surv(clinical[, 1], clinical[, 2]) ~ label, data = data)
    fit2_sum <- summary(fit2)

    # Calculate the concordance index (C-index)
    c_index <- fit2_sum$concordance[1]
    c_index <- round(digits = 3, c_index)

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

    # Return cindex and p value
    return_list <- list()
    return_list[["logrank_p"]] <- logrank_pvalue
    return_list[["c-index"]] <- c_index

    return(return_list)
}
