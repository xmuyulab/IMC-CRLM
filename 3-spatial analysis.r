# Spatial Analysis
library(SingleCellExperiment)

source("./spatial_analysis_functions.r")

## Permutation Result Analysis
FigurePath <- "/mnt/data/lyx/IMC/analysis/spatial/"

clinical <- read.csv("/mnt/data/lyx/IMC/clinical.csv")
head(clinical, 2)

sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")


for (tissue in c("CT", "IM", "TAT")) {
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

    ### All ROI
    if (T) {
        list_ <- getResult(ResultPath, GroupInfo, clinicalGroup = c(0, 1), celltypes = celltypes)
        MergeDF2 <- list_[[1]]
        labelDF2 <- list_[[2]]
        rm(list_)

        ## Single group heatmap
        SingleHeatmap(MergeDF2, labelDF2, savePath = paste0(FigurePath, "All ROIs from ", tissue, " spatial heatmap.pdf"))
    }

    # BubblePlot(MergeDF2, labelDF2, savePath = paste0(FigurePath, "NonRecurrence_spatial_", tissue, ".pdf"))
    DoubleHeat(MergeDF1, labelDF1, group1 = "Relapse", MergeDF2, labelDF2, group2 = "Non-Relapse", plot = "circle", savePath = paste0(FigurePath, "Relapse_spatial_", tissue, ".pdf"))

    ### Different
    getInteracDiff(ResultPath, sce_, GroupInfo,
        groups = c(1, 0),
        celltypes = celltypes, savepath = paste0(FigurePath, "Relapse_Spatial_Diff_", tissue, ".pdf"),
        do.fdr = FALSE, FC_threshold = 0
    )
}
