# scRNA-seq data (GSE178318) CRC-liver metastasis

library(Seurat)
library(clustree)
source("./scRNA_functions.r")

scRNAdataPath <- "/mnt/data/lyx/IMC/RNA/sc/"
savePath <- "/mnt/data/lyx/IMC/RNA/"

## Load data , save Seurat object
LoadandSave_GSE178318(scRNAdataPath, savePath)

GSE178318seu <- readRDS(paste0(savePath, "raw_scSeuratobj.rds"))

## CRC-LM
GSE178318seu <- subset(GSE178318seu, Tissue == "LM")
table(GSE178318seu@meta.data$PID) ## 6 patients
max(GSE178318seu@assays$RNA@data) ## raw counts

## Seurat Analysis Pipeline
### Preprocessing
GSE178318seu[["percent.mt"]] <- PercentageFeatureSet(GSE178318seu, pattern = "^MT-")

p <- VlnPlot(GSE178318seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pdf(paste0(savePath, "before_dataQC.pdf"), height = 6, width = 8)
print(p)
dev.off()

GSE178318seu <- subset(GSE178318seu, subset = nFeature_RNA > 500 & nFeature_RNA < 4000)
GSE178318seu <- subset(GSE178318seu, subset = percent.mt < 15)

p <- VlnPlot(GSE178318seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pdf(paste0(savePath, "after_dataQC.pdf"), height = 6, width = 8)
print(p)
dev.off()


### Normalizing , calculated hvgs
GSE178318seu <- NormalizeData(GSE178318seu, normalization.method = "LogNormalize", scale.factor = 10000)
GSE178318seu <- FindVariableFeatures(GSE178318seu, selection.method = "vst", nfeatures = 2000)

## dimention reduction, clustering , annotation
GSE178318seu <- ScaleData(GSE178318seu)
GSE178318seu <- RunPCA(GSE178318seu)
p <- ElbowPlot(GSE178318seu)
pdf(paste0(savePath, "ElbowPlot.pdf"), height = 6, width = 8)
print(p)
dev.off()

GSE178318seu <- FindNeighbors(GSE178318seu, dims = 1:15)

### clustree
GSE178318seu.clustree <- FindClusters(object = GSE178318seu, resolution = c(seq(0.2, 1.6, 0.2)))
p <- clustree(GSE178318seu.clustree@meta.data, prefix = "RNA_snn_res.")
pdf(paste0(savePath, "clustree.pdf"), height = 8, width = 10)
print(p)
dev.off()

GSE178318seu <- FindClusters(GSE178318seu, resolution = 0.6)
GSE178318seu <- RunUMAP(GSE178318seu, dims = 1:15)
p <- DimPlot(GSE178318seu, reduction = "umap")

pdf(paste0(savePath, "cluterid_UMAP.pdf"), height = 6, width = 8)
print(p)
dev.off()

saveRDS(GSE178318seu, paste0(savePath, "Seurat_CRCLM.rds"))

#### Markers for major type
if (F) {
    EPCs <- c("EPCAM")
    T_cells <- c("CD3D", "CD3G", "TRAC")
    B_cells <- c("CD19", "CD79A", , "MS4A1")
    Plasma_cells <- c("IGHG1", "IGHA1", "MZB1", "CD79A")
    Myeloid_cells <- c("CD68", "CD163", "CD14", "LYZ")
    NK_cells <- c("KLRF1", "KLRD1", "FGFBP2", "PRF1")
    CAFs <- c("FAP", "COL1A1", "COL3A1", "DCN", "ACTA2")
    Endothelial_cells <- c("CLDN5", "CDH5", "VMF")
    pDC <- c("LILRA4", "IL3RA")
    Mast_cells <- c("TPSAB1", "TPSB2", "MS4A2")
}


# plot makers
p <- DotPlot(seuratobj, features = c(
    "KIT", "CPA3", ## Mast cells
    "CD79A", "CD79B", ## B cells
    "CD3D", "CD3G", ## T cells
    "CD163", "CD68", ## Myeloid cells
    "EPCAM", ## Epithelial cells
    "COL3A1", ## Fibroblasts
    "PECAM1" ## Endothelial cells
))
pdf(paste0(savePath, "GSE178318 MajorType makers FeaturePlot.pdf"), width = 10, height = 7.5)
print(p)
dev.off()
