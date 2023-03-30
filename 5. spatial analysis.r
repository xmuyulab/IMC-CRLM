# Spatial Analysis
library(SingleCellExperiment)

source("./spatial_analysis_functions.r")

## Permutation Result Analysis
FigurePath <- "/mnt/data/lyx/IMC/analysis/spatial/"

clinical <- read.csv("/mnt/data/lyx/IMC/clinical.csv")
head(clinical, 2)

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")


for (tissue in c("CT", "IM")) {
    ResultPath <- paste0("/mnt/data/lyx/IMC/analysis/spatial/permutation_", tissue, "/")

    sce_ <- sce[, sce$Tissue == tissue]

    GroupInfo <- GetGroupInfo(sce_, clinical)
    celltypes <- names(table(sce_$SubType))

    ### Relapse
    list_ <- getResult(ResultPath, GroupInfo, clinicalGroup = 1, celltypes = celltypes)
    MergeDF1 <- list_[[1]]
    labelDF1 <- list_[[2]]
    rm(list_)

    # BubblePlot(MergeDF1, labelDF1, savePath = paste0(FigurePath, "Recurrence_spatial_", tissue, ".pdf"))

    ### Non-relapse
    list_ <- getResult(ResultPath, GroupInfo, clinicalGroup = 0, celltypes = celltypes)
    MergeDF2 <- list_[[1]]
    labelDF2 <- list_[[2]]
    rm(list_)

    # BubblePlot(MergeDF2, labelDF2, savePath = paste0(FigurePath, "NonRecurrence_spatial_", tissue, ".pdf"))
    DoubleHeat(MergeDF1, labelDF1, group1 = "Relapse", MergeDF2, labelDF2, group2 = "Non-Relapse", plot = "circle", savePath = paste0(FigurePath, "Relapse_spatial_", tissue, ".pdf"))

    ### Different
    getInteracDiff(ResultPath, sce_, GroupInfo,
        groups = c(1, 0),
        celltypes = celltypes, savepath = paste0(FigurePath, "Relapse_Spatial_Diff_", tissue, ".pdf")
    )
}


## COX for interaction
## no necessary, cellular neighborhood pattern to analysis

## Visualize
if(F){
    sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")

tissue <- "IM"

sce_ <- sce[, sce$Tissue == tissue]

ROIs <- names(table(sce_$ID))

CellMaskPath <- "/mnt/data/lyx/IMC/unet/predictionProbability/ori_cellmask/"
ChannelPath <- "/mnt/data/lyx/IMC/IMCell_Output/save_img/raw/"

channels2plot <- c("CD57", "FoxP3", "CD279")
celltypes2plot <- c("Mono_CD57")

# ROIs <- c("B13_ROI7","B13_ROI8","B13_ROI9","B16_ROI5","B16_ROI9","B16_ROI12","B16_ROI13")

for (ROI in ROIs) {
    VisTypeMaskChannel(sce_, ROI, celltypes2plot, channels2plot, CellMaskPath, ChannelPath, SavePath = paste0("/mnt/data/lyx/IMC/analysis/spatial/imgOutput/", tissue, "/"))
}
}

