
# Spatial Proteomic Decipher Treatment Response in Colorectal Cancer Liver Metastases

## Abstract
Adjuvant therapy (AT) is a prevalent post-operative approach for colorectal cancer liver metastasis (CRLM), yet a significant number of patients could still experience early relapse after AT. In this study, we utilized imaging mass cytometry to construct a spatial proteomic landscape of the tumor microenvironment (TME) of 35 colorectal cancer liver metastasis (CRLM) patients, distinguished by their early relapse status post-adjuvant therapy (AT). Our single-cell analysis identified a significant correlation between early relapse and the enrichment of a distinct subset of CD279+ regulatory T cells at the invasive margin of the tumor. Additionally, spatial proteomic analysis has enabled us to identify two distinct bile duct microenvironments (BDMEs) within the peritumoral region of CRLM responding differently to AT. Leveraging on this finding, we developed a biomarker, termed ‘BDME score’, based on computational analysis on Hematoxylin and Eosin (H&E) stained pathological images. Our findings indicate that this score is positively associated with reduced recurrence-free time post-AT. Furthermore, patients with high BDME scores could potentially derive substantial benefit from a combined treatment strategy of chemotherapy and targeted therapy, indicating the potential of introducing more effective post-surgical treatments for CRLM patients based on improved patient stratification strategies.

## Repository Overview
This repository contains the data analysis pipeline and scripts used in our study titled "Spatial Proteomic Decipher Treatment Response in Colorectal Cancer Liver Metastases." The analysis includes the construction of the spatial proteomic landscape of the tumor microenvironment, single-cell analysis, and the development of the BDME score.

## Table of Contents
1. [Installation](#installation)
2. [Data Preparation](#data-preparation)
3. [Analysis Pipeline](#analysis-pipeline)
4. [Results](#results)
5. [References](#references)

## Installation
To replicate our analysis, ensure you have the following dependencies installed:

```r
# R packages
install.packages(c("SingleCellExperiment", "dplyr", "stringr", "reshape2", "parallel", 
                   "pheatmap", "ComplexHeatmap", "ggplot2", "ggpubr", "ggrepel", 
                   "cowplot", "RColorBrewer", "ggsci", "randomForest"))

# Install additional R packages if needed
# source("/home/lyx/project/IMC/clustering_functions.r")
# source("/home/lyx/project/IMC/FlowSOM_metaClustering.r")
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

### Random Forest Classifier
- Train a Random Forest classifier to predict patient relapse based on proteomic profiles.
- Evaluate the model using confusion matrix and save the trained model.

### T-SNE Visualization
- Perform T-SNE dimensionality reduction on the expression data.
- Generate and save T-SNE plots for visualizing different subtypes.

## Results
- The analysis identifies a correlation between early relapse and CD279+ regulatory T cells.
- Distinct bile duct microenvironments (BDMEs) are identified, leading to the development of the BDME score.
- High BDME scores are associated with reduced recurrence-free time post-AT.
- The BDME score suggests potential benefits from combined chemotherapy and targeted therapy for patients with high scores.

## References
For further details, refer to the original study and supplementary materials provided in this repository. Relevant scripts and data can be found in the respective directories.

## Contact
For any questions or contributions, please contact the study authors at [contact@example.com].
