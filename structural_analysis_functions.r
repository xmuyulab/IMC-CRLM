# Functions for structural analysis

## clinical information
load_clinical <- function(sce, clinicalFilePath) {
    clinical <- read.csv(clinicalFilePath)
    IDs <- names(table(sce$ID))

    ## match the IDs and response level
    GroupInfo <- data.frame(row.names = IDs)

    ### Patient ID
    GroupInfo["PID"] <- sapply(rownames(GroupInfo), function(x) {
        strsplit(x, split = "_")[[1]][1]
    })
    ### RFSS
    GroupInfo["RFS_status"] <- sapply(GroupInfo$PID, function(x) {
        clinical[clinical$PID == x, "RFS_status"]
    })
    GroupInfo["RFS_time"] <- sapply(GroupInfo$PID, function(x) {
        clinical[clinical$PID == x, "RFS_time"]
    })
    ### Mutation
    GroupInfo["KRAS_mutation"] <- sapply(GroupInfo$PID, function(x) {
        clinical[clinical$PID == x, "KRAS_mutation"]
    })
    GroupInfo["Pathology"] <- sapply(GroupInfo$PID, function(x) {
        clinical[clinical$PID == x, "Pathology"]
    })
    GroupInfo["CRC_site"] <- sapply(GroupInfo$PID, function(x) {
        clinical[clinical$PID == x, "CRC_site"]
    })

    return(GroupInfo)
}

## Merge abundance information
MergeAbundanceResult <- function(sce) {
    ## celltypes and ROIs
    celltypes <- names(table(sce$SubType))
    ROIs <- names(table(sce$ID))

    AbundanceDF <- matrix(data = 0, nrow = length(celltypes), ncol = length(ROIs))
    AbundanceDF <- as.data.frame(AbundanceDF)
    rownames(AbundanceDF) <- celltypes
    colnames(AbundanceDF) <- ROIs

    for (i in 1:length(ROIs)) {
        ROI <- ROIs[i]
        sceTemp <- sce[, sce$ID == ROI]
        abundanceTemp <- as.data.frame(table(sceTemp$SubType))

        for (j in 1:nrow(abundanceTemp)) {
            rowTemp <- as.character(abundanceTemp[j, 1])
            AbundanceDF[rowTemp, i] <- as.numeric(abundanceTemp[j, 2])
        }
    }

    return(AbundanceDF)
}
