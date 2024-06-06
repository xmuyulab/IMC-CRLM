
# Spatial Proteomic Decipher Treatment Response in Colorectal Cancer Liver Metastases

## Repository Overview
This repository contains the data analysis pipeline and scripts used in our study titled "Spatial Proteomic Decipher Treatment Response in Colorectal Cancer Liver Metastases." The analysis includes the construction of the spatial proteomic landscape of the tumor microenvironment, single-cell analysis, and the development of the BDME score.

## Table of Contents
1. [Installation](#installation)
2. [Data Preparation](#data-preparation)
3. [Analysis Pipeline](#analysis-pipeline)

## Installation
To replicate our analysis, ensure you have the following dependencies installed:

```r
# R packages
install.packages(c("SingleCellExperiment", "dplyr", "stringr", "reshape2", "parallel", 
                   "pheatmap", "ComplexHeatmap", "ggplot2", "ggpubr", "ggrepel", 
                   "cowplot", "RColorBrewer", "ggsci", "randomForest"))

# Install additional R packages if needed
# source("./clustering_functions.r")
# source("./FlowSOM_metaClustering.r")
```

```python
# Python packages
pip install scanpy pandas numpy scimap matplotlib scikit-learn joblib
```

## Data Preparation
1. **Protein Expression Matrix and Meta Data**:
   - Load the protein expression matrix and meta data from the specified directories.
   - Perform initial data quality checks.

2. **Marker Data**:
   - Load marker data used for normalization and clustering.

3. **Normalization**:
   - Normalize the data with and without arcsinh transformation.

## Analysis Pipeline
### Major Clustering
- Perform major clustering on the normalized data.
- Save clustering results and generate heatmaps.

### Minor Clustering
- Perform minor clustering for specific cell types (e.g., Myeloid cells, Tumor cells).
- Save minor clustering results and generate corresponding plots.

### Downstrem analysis

In this study, we conducted comprehensive downstream analysis to unravel the complexities of the tumor microenvironment (TME) in colorectal cancer liver metastasis (CRLM) patients undergoing adjuvant therapy (AT). Our analysis focused on identifying cellular and spatial features associated with early recurrence post-AT and developing predictive models to guide precision treatment strategies. Below are the key components and findings of our downstream analysis:

- Single-Cell Spatial Proteomic Landscape: Using multiplex single-cell image analysis, we created the first systematic single-cell landscape of the CRLM TME. This high-resolution analysis enabled us to identify cellular and spatial features contributing to early recurrence from AT in clinical biopsy samples.

- Early Recurrence and Immune Suppression: Our single-cell analysis identified a significant correlation between early relapse and a distinct subset of CD279+ regulatory T cells (Tregs) located primarily at the invasive margin of the tumor. These Tregs interact strongly with CD163+ macrophages and CD274+ dendritic cells, suggesting their critical role in maintaining an immune-suppressive microenvironment.

- Identification of Bile Duct Microenvironments (BDMEs): We discovered two distinct BDMEs in the peritumoral region of CRLM patients. The inflammatory BDME, characterized by increased size and upregulated expression of CAIX in bile duct cells, is indicative of hypoxic conditions and reduced AT efficacy. The ratio of the sizes of these BDME types was found to be associated with poor recurrence-free survival (RFS) in CRLM patients undergoing AT.

- BDME Score Development: Building on our spatial proteomic analysis, we developed a BDME score through computational analysis of whole slide images (WSI). This deep learning model serves as an alternative to imaging mass cytometry (IMC) for predicting outcomes in CRLM patients. The BDME score accurately stratifies CRLM patients likely to experience early relapse and distinguishes those who would benefit more from targeted therapy following chemotherapy.

Our downstream analysis provides a detailed understanding of the TME dynamics in CRLM patients, identifies individuals resistant to AT, and proposes more precise therapeutic strategies to improve patient outcomes.

## Contact
For any questions or contributions, please contact the study authors at **Issues**
