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
set.seed(619)
GSE178318seu <- RunPCA(GSE178318seu)
p <- ElbowPlot(GSE178318seu)
pdf(paste0(savePath, "ElbowPlot.pdf"), height = 6, width = 8)
print(p)
dev.off()

GSE178318seu <- FindNeighbors(GSE178318seu, dims = 1:15)

### clustree
GSE178318seu.clustree <- FindClusters(object = GSE178318seu, resolution = c(seq(0.2, 1.6, 0.2)))
p <- clustree(GSE178318seu.clustree@meta.data, prefix = "RNA_snn_res.")
pdf(paste0(savePath, "GSE178318 MajorType clustree.pdf"), height = 8, width = 10)
print(p)
dev.off()

GSE178318seu <- FindClusters(GSE178318seu, resolution = 0.4)
GSE178318seu <- RunUMAP(GSE178318seu, dims = 1:15)
p <- DimPlot(GSE178318seu, reduction = "umap")

ggsave(p, file = paste0(savePath, "GSE178318 clusterid.pdf"), width = 12, height = 10)
print(p)
dev.off()

saveRDS(GSE178318seu, paste0(savePath, "Seurat_CRCLM.rds"))

### Annotation
seuratobj <- readRDS(paste0(savePath, "Seurat_CRCLM.rds"))

#### Markers in original paper
if (F) {
    EPCs <- c("EPCAM")
    T_cells <- c("CD3D", "CD3G", "TRAC")
    B_cells <- c("CD19", "CD79A", "MS4A1")
    Plasma_cells <- c("IGHG1", "IGHA1", "MZB1", "CD79A")
    Myeloid_cells <- c("CD68", "CD163", "CD14", "LYZ")
    NK_cells <- c("KLRF1", "KLRD1", "FGFBP2", "PRF1")
    CAFs <- c("FAP", "COL1A1", "COL3A1", "DCN", "ACTA2")
    Endothelial_cells <- c("CLDN5", "CDH5", "VMF")
    pDC <- c("LILRA4", "IL3RA")
    Mast_cells <- c("TPSAB1", "TPSB2", "MS4A2")
}


#### plot makers
new.cluster.ids <- c(
    "T cell", "T cell", "Myeloid cell", "T cell", "T cell", "Epithelial cell",
    "T cell", "B cell", "B cell", "Stromal", "Unknown cell"
)

### add cluster id and cluster into metadata
seuratobj@meta.data$cell_typeid <- as.character(seuratobj@active.ident)
seuratobj@meta.data$Cell_subtype <- NA

names(new.cluster.ids) <- levels(seuratobj)
seuratobj <- RenameIdents(seuratobj, new.cluster.ids)

Cell_type <- as.character(seuratobj@active.ident)
seuratobj@meta.data$Cell_type <- Cell_type

### UMAP
p <- DimPlot(seuratobj, reduction = "umap", label = T, repel = F, pt.size = 0.5, label.size = 5, combine = T) +
    theme(plot.title = element_text(hjust = 0.5))
ggsave(p, file = paste0(savePath, "GSE178318 MajorType annotation.pdf"), width = 12, height = 10)

### Marker Dotplot
if (1) {
    markers_list <- list()
    markers_list[["B cells"]] <- c("MS4A1", "CD79A")
    markers_list[["T cells"]] <- c("CD3D", "CD3G")
    markers_list[["Myeloid cells"]] <- c("CD163", "LYZ")
    markers_list[["Epithelial cells"]] <- c("EPCAM")
    markers_list[["Fibroblasts"]] <- c("ACTA2", "COL3A1")
}
pdf(file = paste0(savePath, "GSE178318 MajorType maker dotplot.pdf"), width = 20, height = 7.5)
p <- DotPlot(seuratobj, features = markers_list) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()

### T cell
T_seurat <- subset(seuratobj, subset = Cell_type == "T cell")
set.seed(619)
T_seurat <- RunPCA(object = T_seurat)
T_seurat <- FindNeighbors(T_seurat, dims = 1:10)

test_seurat <- FindClusters(object = T_seurat, resolution = c(seq(0, 1.8, .2)))
p <- clustree(test_seurat@meta.data, prefix = "RNA_snn_res.")
pdf(paste0(savePath, "GSE178318 Tcells clustree.pdf"), width = 9, height = 12)
print(p)
dev.off()

T_seurat <- FindClusters(T_seurat, resolution = 0.4)
T_seurat <- RunUMAP(T_seurat, dims = 1:10)
## rename clusters
new.cluster.ids <- c(
    "TRM-like", "TRM-like", "TEM", "TEM", "CD8+ Effector", "Treg",
    "NK", "CD8+ Effector", "CD8+ Effector", "Tfh", "Treg"
)

names(new.cluster.ids) <- levels(T_seurat)
T_seurat <- RenameIdents(T_seurat, new.cluster.ids)
Cell_subtype <- as.character(T_seurat@active.ident)
T_seurat@meta.data$Cell_subtype <- Cell_subtype

## save the cell subtype information into whole seurat object
seuratobj@meta.data$Cell_subtype[match(rownames(T_seurat@meta.data), rownames(seuratobj@meta.data))] <- Cell_subtype

### umap plot
p1 <- DimPlot(T_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T) +
    theme(plot.title = element_text(hjust = 0.5))
pdf(paste0(savePath, "GSE178318 Tcells annotation.pdf"), width = 12, height = 7.5)
print(p1)
dev.off()

### makers dotplot
if (1) {
    markers_list <- list()
    markers_list[["CD4+"]] <- c("CD4")
    markers_list[["CD8+"]] <- c("CD8A", "CD8B")
    markers_list[["Treg"]] <- c("IL2RA", "FOXP3")
    markers_list[["Effector"]] <- c("GZMA", "GZMB", "GZMK", "IFNG")
    markers_list[["TRM-like"]] <- c("CD69", "IL7R", "CXCR4")
    markers_list[["Tfh"]] <- c("CXCL13", "PDCD1")
    markers_list[["NK"]] <- c("KLRF1", "KLRD1", "FGFBP2")
    markers_list[["Th17"]] <- c("IL17A", "IL23R", "RORA")
}
pdf(file = paste0(savePath, "GSE178318 Tcells maker dotplot.pdf"), width = 15, height = 7.5)
p <- DotPlot(T_seurat, features = markers_list) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()

### B cell
B_seurat <- subset(seuratobj, subset = Cell_type == "B cell")
set.seed(619)
B_seurat <- RunPCA(object = B_seurat)
B_seurat <- FindNeighbors(B_seurat, dims = 1:10)

test_seurat <- FindClusters(object = B_seurat, resolution = c(seq(0, 1.2, .2)))
p <- clustree(test_seurat@meta.data, prefix = "RNA_snn_res.")
pdf(paste0(savePath, "GSE178318 Bcells clustree.pdf"), width = 9, height = 12)
print(p)
dev.off()

B_seurat <- FindClusters(B_seurat, resolution = 0.2)
B_seurat <- RunUMAP(B_seurat, dims = 1:10)

## rename clusters
new.cluster.ids <- c(
    "IgG+ Plasma cell", "CD20+ B cell", "CD20+ B cell", "IgA+ Plasma cell", "Naive B cell", "CD20+ B cell"
)

names(new.cluster.ids) <- levels(B_seurat)
B_seurat <- RenameIdents(B_seurat, new.cluster.ids)
Cell_subtype <- as.character(B_seurat@active.ident)
B_seurat@meta.data$Cell_subtype <- Cell_subtype

## save the cell subtype information into whole seurat object
seuratobj@meta.data$Cell_subtype[match(rownames(B_seurat@meta.data), rownames(seuratobj@meta.data))] <- Cell_subtype

## umap plot
pdf(file = paste0(savePath, "GSE178318 Bcells annotation.pdf"), width = 12, height = 7.5)
p1 <- DimPlot(B_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T) +
    theme(plot.title = element_text(hjust = 0.5))
print(p1)
dev.off()

### makers dotplot
if (1) {
    markers_list <- list()
    markers_list[["CD20+"]] <- c("MS4A1")
    markers_list[["IgA+"]] <- c("IGHA1")
    markers_list[["IgG+"]] <- c("IGHG1")
    markers_list[["Naive"]] <- c("IGHD")
}

pdf(file = paste0(savePath, "GSE178318 Bcells maker dotplot.pdf"), width = 10, height = 7.5)
p <- DotPlot(B_seurat, features = markers_list) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()

### Stromal
Stromal_seurat <- subset(seuratobj, subset = Cell_type == "Stromal cell")
set.seed(619)

Stromal_seurat <- RunPCA(object = Stromal_seurat)
Stromal_seurat <- FindNeighbors(Stromal_seurat, dims = 1:15)

test_seurat <- FindClusters(object = Stromal_seurat, resolution = c(seq(0, 1.8, .2)))
p <- clustree(test_seurat@meta.data, prefix = "RNA_snn_res.")
pdf(paste0(savePath, "GSE178318 Stromal clustree.pdf"), width = 9, height = 12)
print(p)
dev.off()

Stromal_seurat <- FindClusters(Stromal_seurat, resolution = 0.6)
Stromal_seurat <- RunUMAP(Stromal_seurat, dims = 1:15)

## rename clusters
new.cluster.ids <- c(
    "Stromal", "Pericyte", "Endothelial cell", "myCAF", "mCAF", "Fibroblast"
)

names(new.cluster.ids) <- levels(Stromal_seurat)
Stromal_seurat <- RenameIdents(Stromal_seurat, new.cluster.ids)
Cell_subtype <- as.character(Stromal_seurat@active.ident)
Stromal_seurat@meta.data$Cell_subtype <- Cell_subtype

## save the cell subtype information into whole seurat object
seuratobj@meta.data$Cell_subtype[match(rownames(Stromal_seurat@meta.data), rownames(seuratobj@meta.data))] <- Cell_subtype

## umap plot
pdf(file = paste0(savePath, "GSE178318 Stromal annotation.pdf"), width = 10, height = 7.5)
p1 <- DimPlot(Stromal_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T) +
    theme(plot.title = element_text(hjust = 0.5))
print(p1)
dev.off()

## makers dotplot
if (1) {
    markers_list <- list()
    markers_list[["Fibroblast"]] <- c("OGN")
    markers_list[["myCAF"]] <- c("ACTA2", "MYLK")
    markers_list[["Pericyte"]] <- c("RGS5", "CSPG4", "ABCC9")
    markers_list[["mCAF"]] <- c("POSTN", "COL5A1")
    markers_list[["Endothelial"]] <- c("PECAM1", "PLVAP")
    markers_list[["other"]] <- c("FAP", "COL3A1", "VEGFA", "PRELID1")
}
p <- DotPlot(Stromal_seurat, features = markers_list) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
pdf(file = paste0(savePath, "GSE178318 Stromal maker dotplot.pdf"), width = 10, height = 7.5)

print(p)
dev.off()

### Myeloid
Myeloid_seurat <- subset(seuratobj, subset = Cell_type == "Myeloid cell")

Myeloid_seurat <- RunPCA(object = Myeloid_seurat)

Myeloid_seurat <- FindNeighbors(Myeloid_seurat, dims = 1:10)

test_seurat <- FindClusters(object = Myeloid_seurat, resolution = c(seq(0, 1.8, .2)))
p <- clustree(test_seurat@meta.data, prefix = "RNA_snn_res.")
pdf(paste0(savePath, "GSE178318 Myeloids clustree.pdf"), width = 9, height = 12)
print(p)
dev.off()

Myeloid_seurat <- FindClusters(Myeloid_seurat, resolution = 0.8)
Myeloid_seurat <- RunUMAP(Myeloid_seurat, dims = 1:10)

## rename clusters
new.cluster.ids <- c(
    "Mono_Claasic", "pDC", "Macro_CD14", "Mono_Claasic", "cDC_CD1C", "Macro_SIGLEC1",
    "Macro_CD14", "Macro_SPP1", "Macro_HLADA", "Mono_CD11c", "Macro_CD11b",
    "cDC_CLE79A"
)

names(new.cluster.ids) <- levels(Myeloid_seurat)
Myeloid_seurat <- RenameIdents(Myeloid_seurat, new.cluster.ids)
Cell_subtype <- as.character(Myeloid_seurat@active.ident)
Myeloid_seurat@meta.data$Cell_subtype <- Cell_subtype

### save the cell subtype information into whole seurat object
seuratobj@meta.data$Cell_subtype[match(rownames(Myeloid_seurat@meta.data), rownames(seuratobj@meta.data))] <- Cell_subtype

### umap plot
pdf(file = paste0(savePath, "GSE178318 Myeloid annotation.pdf"), width = 12, height = 7.5)
p1 <- DimPlot(Myeloid_seurat, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, combine = T) +
    theme(plot.title = element_text(hjust = 0.5))
print(p1)
dev.off()

## makers
pdf(file = paste0(savePath, "GSE178318 Myeloid maker dotplot.pdf"), width = 10, height = 7.5)
if (1) {
    makers_list <- list()
    makers_list[["Macrophages"]] <- c("CD68", "SPP1", "HLA-DRA")
    makers_list[["Monocytes"]] <- c("CD14", "FCGR3A")
    makers_list[["Myeloids"]] <- c("ITGAX", "ITGAM", "SIGLEC1")
    makers_list[["DCs"]] <- c("CLEC9A", "CD1C", "LILRA4", "IL3RA")
}
p <- DotPlot(Myeloid_seurat, features = makers_list) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()

### Epithelial
Epithelial_seurat <- subset(seuratobj, subset = Cell_type == "Epithelial cell")

Epithelial_seurat <- RunPCA(object = Epithelial_seurat)
Epithelial_seurat <- FindNeighbors(Epithelial_seurat, dims = 1:10)

test_seurat <- FindClusters(object = Epithelial_seurat, resolution = c(seq(0, 1.8, .2)))
p <- clustree(test_seurat@meta.data, prefix = "RNA_snn_res.")
pdf(paste0(savePath, "GSE178318 Epithelial clustree.pdf"), width = 9, height = 12)
print(p)
dev.off()

Epithelial_seurat <- FindClusters(Epithelial_seurat, resolution = 0.8)
Epithelial_seurat <- RunUMAP(Epithelial_seurat, dims = 1:10)

## rename clusters
new.cluster.ids <- c(
    "Mono_Claasic", "pDC", "Macro_CD14", "Mono_Claasic", "cDC_CD1C", "Macro_SIGLEC1",
    "Macro_CD14", "Macro_SPP1", "Macro_HLADA", "Mono_CD11c", "Macro_CD11b",
    "cDC_CLE79A"
)

names(new.cluster.ids) <- levels(Epithelial_seurat)
Epithelial_seurat <- RenameIdents(Epithelial_seurat, new.cluster.ids)
Cell_subtype <- as.character(Epithelial_seurat@active.ident)
Epithelial_seurat@meta.data$Cell_subtype <- Cell_subtype

### save the cell subtype information into whole seurat object
seuratobj@meta.data$Cell_subtype[match(rownames(Epithelial_seurat@meta.data), rownames(seuratobj@meta.data))] <- Cell_subtype

### umap plot
pdf(file = paste0(savePath, "GSE178318 Epithelial annotation.pdf"), width = 12, height = 7.5)
p1 <- DimPlot(Epithelial_seurat, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, combine = T) +
    theme(plot.title = element_text(hjust = 0.5))
print(p1)
dev.off()

## makers
pdf(file = paste0(savePath, "GSE178318 Epithelial maker dotplot.pdf"), width = 10, height = 7.5)
if (1) {
    makers_list <- list()
    makers_list[["Macrophages"]] <- c("EPCAM", "MKI67", "VEGFA", "CA9", "HK2", "FAP", "FASN", "CD80", "CD27", "PDCD1", "CD274")
}
p <- DotPlot(Epithelial_seurat, features = makers_list) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()

table(seuratobj@meta.data$Cell_subtype)
table(seuratobj@meta.data$Cell_type)
saveRDS(seuratobj, paste0(savePath, "Seurat_CRCLM_anno.rds"))

# c(
#   "1", "2", "3", "4", "5",
#   "6", "7", "8", "9", "10",
#  "11", "12", "13", "14", "15",
#   "16", "17", "18", "19", "20",
#  "21", "22", "23", "24", "25",
#  "26", "27", "28", "29"
# )

## Analysis
GSE178318seu <- readRDS(paste0(savePath, "Seurat_CRCLM_anno.rds"))
GSE178318seu
table(GSE178318seu@meta.data$Cell_type)

GSE178318seu_myeloid <- subset(GSE178318seu, Cell_type == "Myeloid cell")

## recualtering
GSE178318seu_myeloid <- RunPCA(GSE178318seu_myeloid)
GSE178318seu_myeloid <- FindNeighbors(GSE178318seu_myeloid, dims = 1:10)
GSE178318seu_myeloid <- FindClusters(GSE178318seu_myeloid, resolution = 0.8)
GSE178318seu_myeloid <- RunUMAP(GSE178318seu_myeloid, dims = 1:10)


p <- DimPlot(GSE178318seu_myeloid, reduction = "umap", label = TRUE)
ggsave(filename = paste0(savePath, "Myeloid reclutering id.pdf"), height = 6, width = 8)

## plot some marker
p <- VlnPlot(GSE178318seu_myeloid, features = c("CD274", "PDCD1", "CD80", "CD1C", "CLEC9A", "CD14", "FCGR3A", "CD68", "ITGAX"))
ggsave(filename = paste0(savePath, "Myeloid reclutering marker Vlnplot.pdf"), height = 9, width = 12)

p <- FeaturePlot(GSE178318seu_myeloid, features = c("CD274", "PDCD1", "CD80", "CD1C", "CLEC9A", "CD14", "FCGR3A", "CD68", "ITGAX"))
ggsave(filename = paste0(savePath, "Myeloid reclutering marker FeaturePlot.pdf"), height = 9, width = 12)
