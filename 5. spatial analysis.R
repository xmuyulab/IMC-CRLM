# Spatial Analysis
library(SingleCellExperiment)

source("./spatial_analysis_functions.r")

## Permutation Result Analysis
ResultPath <- "/mnt/data/lyx/IMC/analysis/spatial/permutation/"
FigurePath <- "/mnt/data/lyx/IMC/analysis/spatial/"

clinical <- read.csv("/mnt/data/lyx/IMC/clinical.csv")
head(clinical, 2)

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
sce_ <- sce[, sce$Tissue == "IM"]

GroupInfo <- GetGroupInfo(sce_, clinical)
celltypes <- names(table(sce_$SubType))

### Relapse
list_ <- getResult(ResultPath, GroupInfo, clinicalGroup = 1, celltypes = celltypes)
MergeDF1 <- list_[[1]]
labelDF1 <- list_[[2]]
rm(list_)

BubblePlot(MergeDF1, LabelDF1, savePath = paste0(FigurePath, "Recurrence_spatial.pdf"))

### Non-relapse
list_ <- getResult(ResultPath, GroupInfo, clinicalGroup = 0, celltypes = celltypes)
MergeDF2 <- list_[[1]]
labelDF2 <- list_[[2]]
rm(list_)

BubblePlot(MergeDF2, LabelDF2, savePath = paste0(FigurePath, "NonRecurrence_spatial.pdf"))

### Different
TestDiff(MergeDF1, MergeDF2, celltypes, savepath = paste0(FigurePath, "Spatial_Diff.pdf"))

## COX for interactino

## Visualize
sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
sce_ <- sce[, sce$Tissue == "IM"]

ROIs <- names(table(sce_$ID))

CellMaskPath <- "/mnt/data/lyx/IMC/unet/predictionProbability/ori_cellmask/"
ChannelPath <- "/mnt/data/lyx/IMC/IMCell_Output/save_img/raw/"

channels2plot <- c("CD14", "CLEC9A", "CD11c")
celltypes2plot <- c("Mono_cDC_ITGAX")

for (ROI in ROIs) {
    VisTypeMaskChannel(sce_, ROI, celltypes2plot, channels2plot, CellMaskPath, ChannelPath, savePath = "/mnt/data/lyx/IMC/analysis/spatial/imgOutput/")
}
