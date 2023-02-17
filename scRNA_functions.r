# functions for scRNA-seq data analysis

library(Seurat)
library(Matrix)

LoadandSave_GSE178318 <- function(scRNAdataPath, savePath) {
    cat("Files in ", scRNAdataPath, ": ", list.files(scRNAdataPath), "\n")

    ## cell barcodes
    barcodes <- read.table(paste0(scRNAdataPath, "GSE178318_barcodes.tsv.gz"), header = FALSE, sep = "_")
    colnames(barcodes) <- c("CellID", "PID", "Tissue")
    barcodes <- as.data.frame(barcodes)
    barcodes$CellID <- sapply(c(1:nrow(barcodes)), function(i) {
        return(paste0(barcodes[i, 1], "_", barcodes[i, 2], "_", barcodes[i, 3]))
    })
    rownames(barcodes) <- barcodes$CellID

    ## Gene names
    genes <- read.table(paste0(scRNAdataPath, "GSE178318_genes.tsv.gz"), header = FALSE)
    colnames(genes) <- c("ENSEMBL_ID", "GENE_SYMBOL")
    genes <- as.data.frame(genes)

    ## Expression Profile
    exp <- Matrix::readMM(paste0(scRNAdataPath, "GSE178318_matrix.mtx.gz"))
    exp <- as.matrix(exp)

    ## remove replicate genes
    genesymbol_replicat <- as.data.frame(table(genes$GENE_SYMBOL))
    genesymbol_replicat <- subset(genesymbol_replicat, Freq > 1)$Var1
    genesymbol_replicat <- as.character(genesymbol_replicat)

    removeidx <- c()
    for (Repgene in genesymbol_replicat) {
        idx <- as.data.frame(genes$GENE_SYMBOL %in% Repgene)
        idx1 <- as.numeric(rownames(subset(idx, idx[, 1] == TRUE)))
        removeidx <- c(removeidx, idx1[2])
        exp[idx1[1], ] <- (exp[idx1[1], ] + exp[idx1[2], ]) / 2
    }


    genes <- genes[-removeidx, ]
    exp <- exp[-removeidx, ]

    ## remove PBMC samples
    idx_PBMC <- as.data.frame(barcodes$Tissue %in% "PBMC")
    idx1 <- as.numeric(rownames(subset(idx_PBMC, idx_PBMC[, 1] == TRUE)))

    exp <- exp[, -idx1]
    barcodes <- barcodes[-idx1, ]

    rownames(exp) <- genes$GENE_SYMBOL
    colnames(exp) <- barcodes$CellID

    ## transform exp into sparse matrix
    sparseExp <- as(exp, "dgCMatrix")

    ## clinical information
    clinical <- read.table(paste0(scRNAdataPath, "clinical.txt"), header = TRUE, sep = "\t")
    colnames(clinical) <- sapply(colnames(clinical), function(x) {
        return(gsub(x, pattern = "\\.", replacement = "_"))
    })

    barcodes <- dplyr::left_join(barcodes, clinical, "PID")
    rownames(barcodes) <- barcodes[, 1]

    seuobj <- CreateSeuratObject(counts = sparseExp, project = "CRCLM", meta.data = barcodes, min.cells = 3, min.features = 200)

    saveRDS(seuobj, paste0(savePath, "raw_scSeuratobj.rds"))
    return(NULL)
}
