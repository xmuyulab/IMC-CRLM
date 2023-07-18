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

## Plot certain cell subpopulations interaction stack plot
if (F) {
    celltypespair <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 2))

    celltypespair <- rbind(celltypespair, c("CD8T", "SC_Vimentin"), c("CD4T", "CD8T"), c("Macro_CD11b", "Macro_CD163")) # 3
    # celltypespair <- rbind(celltypespair, c("CD4T", "TC_Ki67"), c("TC_Ki67", "C_aSMA"), c("TC_Ki67", "TC_CAIX")) # 3
    # celltypespair <- rbind(celltypespair, c("CD8T", "Macro_CD169"), c("CD8T", "Macro_CD163"), c("CD8T", "Mono_CD11c"), c("Macro_CD11b", "CD8T")) # 4

    colnames(celltypespair) <- c("c1", "c2")

    celltypespair$Tissue <- c(
        rep("IM", 3) # , rep("CT", 3), rep("TAT", 4)
    )

    for (tissue in c("IM", "CT", "TAT")) { #
        ResultPath <- paste0("/mnt/data/lyx/IMC/analysis/spatial/permutation_", tissue, "/")
        sce_ <- sce[, sce$Tissue == tissue]
        GroupInfo <- GetGroupInfo(sce_, clinical)

        celltypespairTemp <- subset(celltypespair, Tissue == tissue)
        InterDiffStackplot(
            ResultPath, sce_, GroupInfo,
            groupCol = "RFS_status", groups = c(1, 0), groupsName = c("Relapse", "Non-Relapse"),
            celltypespair = celltypespairTemp, savepath = paste0(FigurePath, tissue, " Interaction difference between R and NR.pdf")
        )
    }
}

## Visualize
if (F) {
    sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
    tissues <- c("IM", "TAT", "CT")

    for (tissue in tissues) {
        celltypespairTemp <- subset(celltypespair, Tissue == tissue)
        if (nrow(celltypespairTemp) == 0) {
            next
        }

        Savepath2 <- paste0(FigurePath, tissue, " Spatial interaction/")
        if (!dir.exists(Savepath2)) {
            dir.create(Savepath2, recursive = T)
        }

        sce_ <- sce[, sce$Tissue == tissue]
        ROIs <- names(table(sce_$ID))

        # ROIs <- c("B13_ROI7","B13_ROI8","B13_ROI9","B16_ROI5","B16_ROI9","B16_ROI12","B16_ROI13")

        for (i in 1:nrow(celltypespairTemp)) {
            celltype <- c(celltypespairTemp[i, 1], celltypespairTemp[i, 2])
            for (ROI in ROIs) {
                # VisTypeMaskChannel(sce_, ROI, celltypes2plot, channels2plot, CellMaskPath, ChannelPath, SavePath = paste0("/mnt/data/lyx/IMC/analysis/spatial/imgOutput/", tissue, "/"))
                VisualizeMarkerAndType(sce = sce_, ROI = ROI, celltype = celltype, whichType = "SubType", marker = NULL, SavePath = Savepath2)
            }
        }
    }
}

## Visualize - KRAS
if (F) {
    sce <- readRDS("/mnt/data/lyx/IMC/analysis/allsce.rds")
    sce <- sce[, sce$Tissue %in% c("IM", "CT", "TAT")]
    ROIs <- names(table(sce$ID))

    marker <- c("CAIX")
    celltype <- c("Tumor")

    ROIList <- list()

    KRAS_Mut <- names(table(sce[, sce$KRAS_mutation == 1]$ID))
    KRAS_WT <- names(table(sce[, sce$KRAS_mutation == 0]$ID))

    CT <- names(table(sce[, sce$Tissue == "CT"]$ID))
    TAT <- names(table(sce[, sce$Tissue == "TAT"]$ID))
    IM <- names(table(sce[, sce$Tissue == "IM"]$ID))

    ROIList[["KRAS Mut"]] <- list("CT" = intersect(CT, KRAS_Mut), "TAT" = intersect(TAT, KRAS_Mut), "IM" = intersect(IM, KRAS_Mut))
    ROIList[["KRAS WT"]] <- list("CT" = intersect(CT, KRAS_WT), "TAT" = intersect(TAT, KRAS_WT), "IM" = intersect(IM, KRAS_WT))

    for (i in 1:length(ROIList)) {
        group <- names(ROIList)[i]
        Savepath <- paste0("/mnt/data/lyx/IMC/analysis/structure/", celltype, "_", marker, "/", group, "/")

        for (j in 1:length(ROIList[[i]])) {
            tissue <- names(ROIList[[i]])[j]
            ROIs <- ROIList[[i]][[j]]
            Savepath2 <- paste0(Savepath, tissue, "/")

            if (!dir.exists(Savepath2)) {
                dir.create(Savepath2, recursive = T)
            }
            for (ROI in ROIs) {
                VisualizeMarkerAndType(sce = sce, ROI = ROI, celltype = celltype, whichType = "MajorType", marker = marker, SavePath = Savepath2)
            }
        }
    }
}
