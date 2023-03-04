library(SingleCellExperiment)

source("./spatial_analysis_functions.r")

## Permutation Result Analysis
ResultPath <- "/mnt/data/lyx/IMC/analysis/spatial/permutation/"
FigurePath <- "/mnt/data/lyx/IMC/analysis/spatial/"

## Visualize
sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")

tissue <- "IM"

sce_ <- sce[, sce$Tissue == tissue]

ROIs <- names(table(sce_$ID))

CellMaskPath <- "/mnt/data/lyx/IMC/unet/predictionProbability/ori_cellmask/"
ChannelPath <- "/mnt/data/lyx/IMC/IMCell_Output/save_img/raw/"

channels2plot <- c("CLEC9A", "CD57", "FoxP3")
celltypes2plot <- c("Mono_CLEC9A", "Mono_Multi")

for (ROI in ROIs) {
    VisTypeMaskChannel(sce_, ROI, celltypes2plot, channels2plot, CellMaskPath, ChannelPath, SavePath = paste0("/mnt/data/lyx/IMC/analysis/spatial/imgOutput/", tissue, "/"))
}
