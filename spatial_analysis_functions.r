## functions for permutation test in R
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(SingleCellExperiment)
library(survival)
library(rtiff)
library(tiff)

### Make groupInfo file
GetGroupInfo <- function(sce, clinical){
    ROIs <- names(table(sce$ID))

    returnDF <- matrix(data = NA, nrow = length(ROIs), ncol = 2)
    colnames(returnDF) <- c("RFS_time", "RFS_status")
    rownames(returnDF) <- ROIs

    for(i in 1:nrow(returnDF)){
        PIDTemp <- strsplit(rownames(returnDF)[i], "_")[[1]][1]
        returnDF[i, 1] <- clinical[clinical$PID == PIDTemp, "RFS_time"]
        returnDF[i,2] <- clinical[clinical$PID==PIDTemp,"RFS_status"]
    }

    return(as.data.frame(returnDF))
}

### get the results of permutation test
getResult <- function(ResultPath, GroupInfo, clinicalGroup, celltypes, p_threshold = 0.05) {
    ## Interaction
    numtypes <- length(celltypes)

    interDF <- matrix(data = 0, nrow = numtypes, ncol = numtypes)
    interDF <- as.data.frame(interDF)
    rownames(interDF) <- celltypes
    colnames(interDF) <- celltypes

    IDsTemp <- rownames(GroupInfo[GroupInfo$RFS_status == clinicalGroup, ])
    numID <- length(IDsTemp)

    for (ID in IDsTemp) {
        filepathTemp <- paste0(ResultPath, ID, "/ClosePvalue/ClosePvalue.csv")
        csvTemp <- read.csv(filepathTemp)
        rownames(csvTemp) <- csvTemp[, 1]
        csvTemp <- csvTemp[, -1]

        csvTemp <- ifelse(csvTemp <= p_threshold, 1, 0)

        colname_ <- colnames(csvTemp)
        rowname_ <- rownames(csvTemp)

        for (i in 1:nrow(csvTemp)) {
            for (j in 1:ncol(csvTemp)) {
                if (csvTemp[i, j] == 1) {
                    interDF[rowname_[i], colname_[j]] <- interDF[rowname_[i], colname_[j]] + 1
                }
            }
        }
    }
    interDF <- round(interDF / numID,4)  

    ## Avoideness
    avoidDF <- matrix(data = 0, nrow = numtypes, ncol = numtypes)
    avoidDF <- as.data.frame(avoidDF)
    rownames(avoidDF) <- celltypes
    colnames(avoidDF) <- celltypes

    IDsTemp <- rownames(GroupInfo[GroupInfo$RFS_status == clinicalGroup, ])

    for (ID in IDsTemp) {
        filepathTemp <- paste0(ResultPath, ID, "/AvoidPvalue/AvoidPvalue.csv")
        csvTemp <- read.csv(filepathTemp)
        rownames(csvTemp) <- csvTemp[, 1]
        csvTemp <- csvTemp[, -1]

        csvTemp <- ifelse(csvTemp <= p_threshold, 1, 0)

        colname_ <- colnames(csvTemp)
        rowname_ <- rownames(csvTemp)

        for (i in 1:nrow(csvTemp)) {
            for (j in 1:ncol(csvTemp)) {
                if (csvTemp[i, j] == 1) {
                    avoidDF[rowname_[i], colname_[j]] <- avoidDF[rowname_[i], colname_[j]] + 1
                }
            }
        }
    }

    avoidDF <- round(avoidDF / numID,4)  

    ## Merge
    labelDF <- matrix(data = 0, nrow = numtypes, ncol = numtypes)
    labelDF <- as.data.frame(labelDF)
    rownames(labelDF) <- celltypes
    colnames(labelDF) <- celltypes

    MergeDF <- labelDF

    for (i in 1:nrow(labelDF)) {
        for (j in 1:ncol(labelDF)) {
            if (interDF[i, j] > avoidDF[i, j]) {
                MergeDF[i, j] <- interDF[i, j]
                labelDF[i, j] <- 1
            } else {
                MergeDF[i, j] <- avoidDF[i, j]
                labelDF[i, j] <- -1
            }
        }
    }

    return(list("MergeDF" = MergeDF, "LabelDF" = labelDF))
}

### Transfrom matrix into plot dataframe
BubblePlot <- function(MergeDF,LabelDF,savePath){

    if (class(MergeDF) == "list"){
        MergeDF <- MergeDF[[1]]
        LabelDF <- LabelDF[[1]]
    }

    nrow_ <- nrow(MergeDF)
    ncol_ <- ncol(MergeDF)

    plotdf <- matrix(data = NA, nrow = nrow_ * ncol_, ncol = 4)
    plotdf <- as.data.frame(plotdf)

    c1 <- rep(rownames(MergeDF),each = nrow_)
    c2 <- rep(colnames(MergeDF),time = ncol_)
    num <- c()
    label <- c()

    for (i in 1:nrow_){
        num <- c(num,as.numeric(MergeDF[i,]))
        label <- c(label,as.numeric(LabelDF[i,]))
    }
    label <- ifelse(label == 1, "Interaction","Avoidence")
    
    colnames(plotdf) <- c("Celltype1", "Celltype2", "ROI Fraction", "Spatial Status")

    plotdf["Celltype1"] <- c1
    plotdf["Celltype2"] <- c2
    plotdf["ROI Fraction"] <- num
    plotdf["Spatial Status"] <- as.factor(label) 

    ### bubble plot
    p <- ggplot(plotdf, aes(x = Celltype1, y = Celltype2, size = `ROI Fraction`, color = `Spatial Status`)) +
        geom_point() +
        theme(
            panel.background = element_blank(),
            panel.grid.major = element_line(colour = "white"),
            panel.border = element_rect(colour = "white", fill = NA),
        axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.25))

    pdf(savePath, width = 8, height = 8)
    print(p)
    dev.off()

    return(NULL)

}

## t-test for interaction number between groups
LoadAnalysisResult <- function(IDs1,celltypes,sce){

    numtypes <- length(celltypes)

    array1 <- array(data = NA,dim = c(numtypes,numtypes,length(IDs1)))
    rownames(array1) <- celltypes
    colnames(array1) <- celltypes

    i = 1
    for (ID in IDs1) {
        sceTemp <- sce[, sce$ID == ID]
        
        ## load result dataframe
        filepathTemp <- paste0(ResultPath, ID, "/TureInteraction/TureInteraction.csv")
        csvTemp <- read.csv(filepathTemp)
        rownames(csvTemp) <- csvTemp[, 1]
        csvTemp <- csvTemp[, -1]
        csvTemp <- as.matrix(csvTemp)

        rownameTemp <- rownames(csvTemp)
        colnameTemp <- colnames(csvTemp)


        TempMat <- matrix(data = NA, nrow = numtypes, ncol = numtypes)
        colnames(TempMat) <- celltypes
        rownames(TempMat) <- celltypes

        ## sort the interaction number
        if (nrow(csvTemp) != numtypes) {
            for (j in 1:nrow(csvTemp)) {
                for (z in 1:ncol(csvTemp)) {
                    c1 <- rownameTemp[j]
                    c2 <- colnameTemp[z]

                    numc1 <- ncol(sceTemp[, sceTemp$SubType == c1])
                    numc2 <- ncol(sceTemp[, sceTemp$SubType == c2])

                    TempMat[c1, c2] = round(csvTemp[j, z] / (numc1+numc2), 6) 
                }
            }
        } 
        else {
            csvTemp <- csvTemp[match(celltypes, rownames(csvTemp)), ]
            TempMat <- csvTemp[, match(celltypes, colnames(csvTemp))]
        }
        
        ## calculate the cell number of certain cell types
        cellnumMat <- matrix(data = 0, nrow = numtypes, ncol = numtypes)
        rownames(cellnumMat) <- celltypes
        colnames(cellnumMat) <- celltypes

        for(j in 1:ncol(cellnumMat)){
            c1 <- colnames(cellnumMat)[j]
            c1num <- ncol(sceTemp[,sceTemp$SubType==c1])
            cellnumMat[,j] <- cellnumMat[,j] + c1num
        }
        for(j in 1:nrow(cellnumMat)){
            c2 <- rownames(cellnumMat)[j]
            c2num <- ncol(sceTemp[, sceTemp$SubType == c2])
            cellnumMat[j,] <- cellnumMat[j,] + c2num
        }

        TempMat <- round(TempMat / cellnumMat,6) 

        array1[, , i] <- TempMat
        i = i + 1
    }
    
    return(array1)
}

TestDiff <- function(array1,array2,celltypes,savepath){
    numtypes <- length(celltypes)
    ResultMat <- matrix(data = NA, nrow = numtypes, ncol = numtypes)
    rownames(ResultMat) <- celltypes
    colnames(ResultMat) <- celltypes

    MaskMat <- ResultMat

    for (i in 1:nrow(ResultMat)) {
        for (j in 1:ncol(ResultMat)) {
            value1 <- as.numeric(na.omit(array1[i, j, ]))
            value2 <- as.numeric(na.omit(array2[i, j, ]))

            ResultMat[i, j] <- t.test(value1, value2)$p.value

            if(mean(value1) < mean(value2)){
                MaskMat[i, j] <- 1
            }
            else {
               MaskMat[i, j] <- -1
            }
        }
    }
    ResultMat <- -log10(ResultMat)
    HeatmapForDiff(ResultMat,MaskMat,savepath)

    return(NULL)
}

getSig <- function(dc) {
  sc <- ''
  if (dc >= 3) sc <- '***'
  else if (dc >= 2) sc <- '**'
  else if (dc >= 1.3) sc <- '*'
  return(sc)
} 

HeatmapForDiff <- function(ResultMat,MaskMat,savepath){
    sig_mat <- matrix(sapply(ResultMat, getSig), nrow = nrow(ResultMat))

    plotdf <- ResultMat * MaskMat
    
    p <- pheatmap(
        plotdf,
        cellwidth = 20, cellheight = 15,
        cluster_row = F, cluster_col = F,
        legend_labels = c("Non-Rec","Rec"),legend_breaks = c(-1.3,1.3),
        angle_col = '45', display_numbers = sig_mat, fontsize_number = 15
    )
    pdf(savepath, width = 15, height = 15)
    print(p)
    dev.off()

    return(NULL)
}

getInteracDiff <- function(ResultPath, sce, GroupInfo, groups, celltypes,savepath){

    IDs1 <- rownames(GroupInfo[GroupInfo$RFS_status == groups[1], ])
    IDs2 <- rownames(GroupInfo[GroupInfo$RFS_status == groups[2], ])

    array1 <- LoadAnalysisResult(IDs1,celltypes, sce)
    array2 <- LoadAnalysisResult(IDs2,celltypes, sce)

    TestDiff(array1,array2,celltypes,savepath)

    return(NULL)
}

## cox test
cox_test <- function(ResultPath, GroupInfo, groups, celltypes,savepath){
    IDs1 <- rownames(GroupInfo[GroupInfo$RFS_status == groups[1], ])
    IDs2 <- rownames(GroupInfo[GroupInfo$RFS_status == groups[2], ])

    array1 <- LoadAnalysisResult(IDs1,celltypes)
    array2 <- LoadAnalysisResult(IDs2,celltypes)

    GroupInfo <- GroupInfo[c(IDs1,IDs2),]

    time <- GroupInfo$RFS_time
    status <- GroupInfo$RFS_status
    value <- c(array1[1, 1, ], array2[1, 1, ])

    data <- matrix(data = NA, nrow = length(time), ncol = 3)
    data <- as.data.frame(data)
    colnames(data) <- c("Number","RFS_time","RFS_status")
    data$Number <- value
    data$RFS_time <- as.numeric(time)
    data$RFS_status <- ifelse(status=="Recurrence",1,0)
    data <- na.omit(data)

    data$Number[1:61] <- sample(1:100,size = 61)
    data$Number[62:115] <- sample(100:200,size = 54)

    fit <- coxph(Surv(RFS_time, RFS_status) ~ Number, data = data)
    summary(fit)

    return(NULL)
}

## visualize functions for sce object

## get cell coordinate
getCoordinate <- function(sce_){
    if(!"Position"%in%colnames(colData(sce_))){
        cat("There is not column named Position!", "\n")
        return(NULL)
    }

    position <- sapply(sce_$Position,function(a){
        strsplit(a,",")
    })
    names(position) <- NULL

    x <- lapply(position,function(b){
        return(as.numeric(strsplit(b[1],"\\(")[[1]][2])) 
    })
    y <- lapply(position,function(b){
        return(as.numeric(strsplit(b[2],"\\)")[[1]][1])) 
    })

    x <- unlist(x)
    y <- unlist(y)

    return(list(x,y))
}

## plot Marker expression value on cell-level 
PlotMarker <- function(sce, ROI, Marker, SavePath){
    if(!Marker %in% rownames(sce)){
        cat(Marker, " is no in sce object!","\n")
        return(NULL)
    }

    sce_ <- sce[,sce$ID == ROI]

    coorList <- getCoordinate(sce_)

    value <- sce_[Marker,]
    value <- as.vector(assay(value))

    plotdf <- as.data.frame(matrix(data = NA, nrow = length(value), ncol = 3))
    colnames(plotdf) <- c("x", "y", "Expression")
    plotdf["x"] <- coorList[[1]]
    plotdf["y"] <- coorList[[2]]
    plotdf["Expression"] <- value

    myPalette <- brewer.pal(2,"Paired")
    p <- ggplot(plotdf, aes(x = x, y = y, color = Expression)) +
        geom_point(size = 1) +
        scale_colour_gradientn(colours = myPalette, limits = c(0, 1)) +
        labs(title = paste0(Marker," expression level in ",ROI))+
        theme_test()
        
    pdf(SavePath, height = 6, width = 8)
    print(p)
    dev.off()

    return(p)
}

## plot celltypes on cell-level 
PlotCelltypes <- function(sce, ROI, selectCelltypes, SavePath){
    #colnames(colData(sce))
    sce_ <- sce[,sce$ID == ROI]

    coorList <- getCoordinate(sce_)

    celltypes <- sce_$SubType
    celltypes <- ifelse(celltypes %in% selectCelltypes, celltypes, "Background")
    celltypes <- as.factor(celltypes)

    idx <- match("Background",levels(celltypes))
    levels(celltypes) <- c('Background',levels(celltypes)[-idx])

    plotdf <- as.data.frame(matrix(data = NA, nrow = length(celltypes), ncol = 4))
    colnames(plotdf) <- c("x", "y", "Identity","value")
    plotdf["x"] <- coorList[[1]]
    plotdf["y"] <- coorList[[2]]
    plotdf["Identity"] <- celltypes 

    myPalette <- brewer.pal(length(selectCelltypes),"Dark2")

    p <- ggplot(plotdf, aes(x = x, y = y, color = Identity)) +
        geom_point(size = 1) +
        scale_colour_manual(values = c("grey",myPalette)) +
        labs(title = paste0(ROI))+
        theme_test()
        
    pdf(SavePath, height = 6, width = 8)
    #pdf("test.pdf", height = 6, width = 8)
    print(p)
    dev.off()

    return(NULL)
}


## calculate the cell distance
getDistance <- function(coorList1,coorList2){

    x1 <- round(coorList1[[1]],6)
    y1 <- round(coorList1[[2]],6)
    x2 <- round(coorList2[[1]],6)
    y2 <- round(coorList2[[2]],6)

    df <- matrix(data = NA, nrow = length(x2), ncol = length(x1))
    rownames(df) <- paste0(x2,",",y2)
    colnames(df) <- paste0(x1,",",y1)

    for(i in 1:ncol(df)){
        x1Temp <- rep(x1[i],times = nrow(df)) 
        y1Temp <- rep(y1[i],times = nrow(df)) 

        dist <- sqrt((x1Temp-x2)^2 + (y1Temp-y2)^2)
        dist <- ifelse(dist==0,999,dist)
        df[, i] <- dist
    }

    return(df)
}

## plot celltype interaction
PlotCelltypeInteraction <- function(sce, ROI, celltypes, radius = 22, savePath){
    sce_ <- sce[, sce$ID == ROI]
    c1 <- celltypes[1]
    c2 <- celltypes[2]

    allcoorList <- getCoordinate(sce_)
    x_all <- round(allcoorList[[1]],6)
    y_all <- round(allcoorList[[2]],6)

    sce_1 <- sce_[,sce_$SubType == c1]
    sce_2 <- sce_[,sce_$SubType == c2]

    if(dim(sce_1)[2] < 5){
        cat(c1," in ",ROI," were less than 5!",'\n')
    }
    if (dim(sce_2)[2] < 5) {
        cat(c2, " in ", ROI, " were less than 5!", "\n")
    }
    
    coorList1 <- getCoordinate(sce_1)
    coorList2 <- getCoordinate(sce_2)

    df <- getDistance(coorList1, coorList2)
    MASK <- ifelse(df<=radius,1,0)

    plotdf <- as.data.frame(matrix(data = NA, nrow = length(x_all), ncol = 3))
    colnames(plotdf) <- c("x", "y", "Subtype")
    plotdf["x"] <- x_all
    plotdf["y"] <- y_all
    plotdf["Subtype"] <- NA

    for (i in 1:nrow(MASK)) {
        for (j in 1:ncol(MASK)) {
            if (MASK[i, j] == 1) {
                c1Temp <- colnames(MASK)[j]
                c2Temp <- rownames(MASK)[i]

                c1x <- as.numeric(strsplit(c1Temp, ",")[[1]][1])
                c1y <- as.numeric(strsplit(c1Temp, ",")[[1]][2])
                c2x <- as.numeric(strsplit(c2Temp, ",")[[1]][1])
                c2y <- as.numeric(strsplit(c2Temp, ",")[[1]][2])

                index1 <- (plotdf$x == c1x) & (plotdf$y == c1y)
                plotdf$Subtype[index1] = c1

                index2 <- (plotdf$x == c2x) & (plotdf$y == c2y)
                plotdf$Subtype[index2] = c2
            }
        }
    }
    plotdf$Subtype[is.na(plotdf$Subtype)] = "Background"
    
    p <- ggplot(plotdf, aes(x = x, y = y, color = Subtype)) +
        geom_point(size = 1) +
        labs(title = paste0("Interaction of ",c1," and ",c2," in ",ROI))+
        theme_test()
        
    pdf(savePath, height = 6, width = 8)
    print(p)
    dev.off()

    return(p)
}

## load cell mask files and rename
LoadCellMask <- function(cellMaskPath, ROI) {

    files <- list.files(cellMaskPath)
    filesROI <- sapply(files,function(x){
        temp <- strsplit(x,split = "_")[[1]]
        paste0(temp[2],"_",temp[3]) 
    })

    MaskPath <- paste0(cellMaskPath,names(filesROI[filesROI==ROI]))

    ## read image
    CellMaskMat <- ReadTiff(MaskPath)

    return(CellMaskMat)
}

## load channel files
LoadChannelImage <- function(channelIamgePath, ROI, channel){
    
    pathTemp1 <- paste0(channelIamgePath, ROI,"/")
    channels <- list.files(pathTemp1)

    channelImages <- sapply(channels, function(x) {
        temp <- strsplit(x, split = "_")[[1]]
        return(temp[3])
    })

    channelsPath <- c()
    for(i in channel){
        pathTemp2 <- paste0(pathTemp1,names(channelImages[channelImages==i]))
        channelsPath <- c(channelsPath,pathTemp2)
    }

    ## Load channel images
    ImageList <- list()
    for (i in 1:length(channelsPath)) {
        imgTemp <- png::readPNG(channelsPath[i])
        imgTemp <- TransfromtoGray(imgTemp)
        #imgTemp <- PercentileFilter(imgTemp,cutoff=0.95)
        imgTemp <- MedianFilter(imgTemp)
        imgTemp <- Transfromto0255(imgTemp)

        ImageList[[i]] <- imgTemp

    }

    if(length(ImageList) == 2){
        ImageList[[3]] <- matrix(data = 0,nrow = nrow(ImageList[[2]]), ncol = ncol(ImageList[[2]]))
    }
    if(length(ImageList) == 1){
        ImageList[[2]] <- matrix(data = 0,nrow = nrow(ImageList[[1]]), ncol = ncol(ImageList[[1]]))
        ImageList[[3]] <- matrix(data = 0,nrow = nrow(ImageList[[1]]), ncol = ncol(ImageList[[1]]))
    }

    ## combine different channel
    imagearray <- array(data = NA, dim = c(nrow(ImageList[[1]]), ncol(ImageList[[1]]), length(ImageList)))

    for (i in 1:length(channelsPath)) {
        imagearray[, , i] <- ImageList[[i]]
    }
    
    #png::writePNG(imgTemp,"median.png")

    return(imagearray)
    
}

## read tiff
ReadTiff <- function(filePath,filter = F){
    tif <- readTIFF(filePath)
    tif <- Transfromto0255(tif)
    if(filter){
        tif <- ImageFilter(tif)
    }
    
    return(tif)
}

## 0 - 255 tranformation
Transfromto0255 <- function(mat){
    max_ <- max(mat)
    min_ <- min(mat)

    mat <- (mat-min_) / (max_-min_)
    return(mat * 255)
}

## Transform RGB image into gray-scale image
TransfromtoGray <- function(array) {
    mat <- array[,,1] + array[,,2] + array[,,3] / 3
    return(mat)
}

## Image filter, remove the value more than percentile 
PercentileFilter <- function(mat, cutoff = 0.99){
    cutoff <- quantile(mat, probs = cutoff)
    ifelse(mat >= cutoff,0,mat)
    return(mat)
}

## Image fileter, median filter
MedianFilter <- function(mat,filterSize = 3){
    rowNum <- nrow(mat)
    colNum <- ncol(mat)

    margin <- (filterSize-1) %/% 2

    for(i in (1+margin):(rowNum-margin)){
        for(j in (1+margin):(colNum-margin)){
            seq_ <- as.numeric(mat[(i - margin):(i + margin), (j - margin):(j + margin)])
            mat[i,j] = median(seq_)
        }
    }
    return(mat)
}

## Visualize the cellsubtype, cell mask and channel of ROI
VisTypeMaskChannel <- function(sce, ROI, celltypes, channel, maskPath, channelPath, SavePath){
    SavePath <- paste0(SavePath,ROI,"/")
    if (!dir.exists(SavePath)) {
        dir.create(SavePath)
    }
    
    MaskMat <- LoadCellMask(maskPath, ROI)
    ChannelArray <- LoadChannelImage(channelPath, ROI, channel)
    PlotCelltypes(sce, ROI, celltypes, paste0(SavePath,celltypes," on cell level.pdf"))

    png::writePNG(MaskMat, paste0(SavePath, "CellMask.png"), dpi = 100)
    png::writePNG(ChannelArray, paste0(SavePath, channel[1], "-", channel[2], "_channel.png"), dpi = 100)
    cat(ROI, ": CellMask, celltypes and channel image were done!",'\n')
    return(NULL)
}