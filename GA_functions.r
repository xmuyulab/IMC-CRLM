## Functions for Tumor-Budding like structure analysis
library(ggplot2)
library(SingleCellExperiment)

## Convert to adjacency matrix
Convert2AdjMat <- function(edges) {
    # Create a dataframe of edges and their lengths
    df <- data.frame(
        from = edges$ind1,
        to = edges$ind2,
        weight = sqrt((edges$x1 - edges$x2)^2 + (edges$y1 - edges$y2)^2)
    )

    # create a square matrix with zero values
    vertices <- unique(c(df[, 1], c(df[, 2])))
    vertices <- sort(vertices)
    n <- length(vertices)
    adjacency_matrix <- matrix(1000, nrow = n, ncol = n)
    for (i in 1:nrow(df)) {
        adjacency_matrix[df[i, 1], df[i, 2]] <- df[i, 3]
        adjacency_matrix[df[i, 1], df[i, 2]] <- df[i, 3] # as the graph is undirected
    }

    return(adjacency_matrix)
}

## Random Walk
RandomWalk <- function(adjacency_matrix, startVertex, steps = 20) {
    ## Record
    PathDistance <- 0
    EachEdgeDistance <- c()
    WalkingNode <- c()

    ## Find the next point via minimum edge
    probaEdge <- adjacency_matrix[startVertex, ]
    NextNode <- which.min(probaEdge)
    PathDistance <- PathDistance + probaEdge[NextNode]
    WalkingNode <- c(WalkingNode, startVertex)
    EachEdgeDistance <- c(EachEdgeDistance, probaEdge[NextNode])

    ## Loop for find path
    for (i in 2:steps) {
        adjacency_matrixTemp <- adjacency_matrix
        ## Forbid to pass same node
        adjacency_matrixTemp[WalkingNode[-1], ] <- 1000
        adjacency_matrixTemp[, WalkingNode[-1]] <- 1000

        lastNode <- WalkingNode[i - 1]

        NowNode <- NextNode
        adjacency_matrixTemp[NowNode, lastNode] <- 1000

        probaEdge <- adjacency_matrixTemp[NowNode, ]
        NextNode <- which.min(probaEdge)
        PathDistance <- PathDistance + probaEdge[NextNode]
        EachEdgeDistance <- c(EachEdgeDistance, probaEdge[NextNode])
        WalkingNode <- c(WalkingNode, NowNode)

        ## Stop conditions
        ## Form Loop
        if (NextNode == startVertex) {
            break
        }
        ### Isolate Point
        # if (1000 %in% EachEdgeDistance) {
        #    break
        # }
    }

    returenList <- list()
    returenList[["PassNodes"]] <- WalkingNode
    returenList[["PassEdgesLength"]] <- EachEdgeDistance
    returenList[["PathLength"]] <- PathDistance
    returenList[["MeanStepLength"]] <- PathDistance / length(WalkingNode)

    return(returenList)
}

## Filtering Path
FilterPath <- function(results, MaxMeanStepLength, MaxPassNodes, MinPassNodes) {
    for (i in length(results):1) {
        if (results[[i]][["MeanStepLength"]] > MaxMeanStepLength) {
            results[[i]] <- NULL
        }
    }
    for (i in length(results):1) {
        if (length(results[[i]][["PassNodes"]]) > MaxPassNodes) {
            results[[i]] <- NULL
        }
    }
    for (i in length(results):1) {
        if (length(results[[i]][["PassNodes"]]) < MinPassNodes) {
            results[[i]] <- NULL
        }
    }
    return(results)
}

## Match the vertex index on coordinate
MatchPointOnCor <- function(results, edges) {
    x <- c()
    y <- c()
    Nodes <- c()

    edges_1 <- edges[, c(1, 2, 5)]
    edges_2 <- edges[, c(3, 4, 6)]
    colnames(edges_2) <- colnames(edges_1)

    edges_ <- rbind(edges_1, edges_2)

    for (i in 1:length(results)) {
        VertexTemp <- results[[i]][["PassNodes"]]
        Nodes <- c(Nodes, VertexTemp)
    }
    Nodes <- unique(Nodes)

    ## Find Coordinate
    for (i in 1:length(Nodes)) {
        idxTemp <- match(Nodes[i], edges_$ind1)
        x <- c(x, edges_[idxTemp, 1])
        y <- c(y, edges_[idxTemp, 2])
    }

    plotdf <- cbind(x = x, y = y)
    plotdf <- as.data.frame(plotdf)

    p <- ggplot(plotdf, aes(x = x, y = y)) +
        geom_point(size = 1) +
        scale_x_continuous(limits = c(0, 1000)) +
        scale_y_continuous(limits = c(0, 1000)) +
        theme_test()

    pdf("test.pdf", height = 6, width = 8)
    print(p)
    dev.off()
}

## Preprocessing the dataframe to plot Volcano graph
PreprocessingVolcan <- function(mat, nameCol, FC_threshold = 1.4, Qvalue_threshold = 0.05) {
    mat$Foldchange <- as.numeric(mat$Foldchange)
    mat$P.value <- as.numeric(mat$P.value)
    mat$Q.value <- as.numeric(mat$Q.value)

    mat$dir <- ""
    for (i in 1:nrow(mat)) {
        if ((mat[i, "Foldchange"] >= FC_threshold) & (mat[i, "Q.value"] <= Qvalue_threshold)) {
            mat$dir[i] <- "up-regulated"
        }
        if ((mat[i, "Foldchange"] <= (1 / FC_threshold)) & (mat[i, "Q.value"] <= Qvalue_threshold)) {
            mat$dir[i] <- "down-regulated"
        }
    }

    mat$label <- ifelse(mat$dir != "", mat[, nameCol], "")

    mycol <- c(ggsci::pal_npg("nrc", alpha = 0.8)(length(names(table(mat$dir))) - 1), "gray")
    names(mycol) <- c(sort(unique(mat$dir))[3:2], "NOT")

    returnList <- list()
    returnList[["mat"]] <- mat
    returnList[["mycol"]] <- mycol

    return(returnList)
}

## Get the cell subpopulation abundance fraction from certain tissue
GetFractionFromTissue <- function(df, tissue, subtype) {
    df2 <- subset(df, Tissue == tissue)
    df2 <- df2[, match(c("PID", subtype), colnames(df2))]
    colnames(df2)[2] <- paste0(tissue, "_", subtype)

    return(df2)
}

## Calcualte the x and y from a position vector
GetXandY <- function(posVec) {
    xVec <- c()
    yVec <- c()

    spliteRes <- sapply(posVec, function(x) {
        return(strsplit(x, split = ",")[[1]])
    })

    for (i in seq(1, length(spliteRes), 2)) {
        xTemp <- strsplit(spliteRes[i], "\\(")[[1]][2]
        yTemp <- strsplit(spliteRes[i + 1], ")")[[1]][1]

        xVec <- c(xVec, as.numeric(xTemp))
        yVec <- c(yVec, as.numeric(yTemp))
    }

    cor <- as.data.frame(matrix(data = NA, ncol = 2, nrow = length(xVec)))
    cor[, 1] <- xVec
    cor[, 2] <- yVec

    colnames(cor) <- c("cor_x", "cor_y")
    return(cor)
}

# Define the dist2 function
dist2 <- function(mat1, mat2) {
    sqrt(outer(mat1[, 1], mat2[, 1], "-")^2 + outer(mat1[, 2], mat2[, 2], "-")^2)
}

# Define the distance to boundary
dist2boundary <- function(x) {
    # Compute the IQR.
    IQR <- IQR(x)

    # Compute the lower and upper bounds for what we consider non-outliers.
    lower_bound <- quantile(x, 0.25) - 1.5 * IQR
    upper_bound <- quantile(x, 0.75) + 1.5 * IQR

    # Remove outliers.
    a_without_outliers <- x[x > lower_bound & x < upper_bound]

    return(mean(a_without_outliers))
}

## Get the k-nn distance to boundary
get_nearest_points <- function(df, targetType, n = 10) {
    # Split the dataframe into target and non-target points
    df_target <- df[df$Type == targetType, ]

    # Check if there's no point of target type
    if (nrow(df_target) == 0) {
        stop("No points of target type found in the dataframe")
    }

    # Calculate the Euclidean distance matrix between target and non-target points
    # Create matrices for the coordinates
    mat_target <- as.matrix(df_target[, c("cor_x", "cor_y")])
    mat_all <- as.matrix(df[, c("cor_x", "cor_y")])

    # Compute distances between every point in df_target and every point in df_other
    distances <- as.matrix(dist2(mat_target, mat_all))

    # Find the indices of the 'n' non-target points with the shortest minimum distances
    sort_distances <- apply(distances, 2, function(x) {
        return(sort(x)[1:n])
    })

    # Get the distance to boundary
    disValue <- apply(sort_distances, 2, function(x) {
        return(dist2boundary(x))
    })

    return(disValue)
}

## Determine the distance to certain border
Dis2Boundary <- function(sce_, targetType, targetTypeCol, imageIDcol = "ID", coorCol = "Position", DisColName = "Dis2Tumor", k = 10) {
    colData(sce_)[, DisColName] <- 0

    ## Get images ID
    ROIs <- names(table(colData(sce_)[, imageIDcol]))

    for (ROI in ROIs) {
        sceTemp <- sce_[, colData(sce_)[, imageIDcol] == ROI]

        ## Calculate each cell distance to boundary
        coorTemp <- GetXandY(colData(sceTemp)[, coorCol])
        coorTemp <- cbind(coorTemp, Type = colData(sceTemp)[, targetTypeCol])
        
        disValue <- get_nearest_points(df = coorTemp, targetType = targetType, n = k)
        colData(sce_)[match(colData(sceTemp)[, imageIDcol], colData(sce_)[, imageIDcol]), ] <- disValue
    }

    return(sce_)
}