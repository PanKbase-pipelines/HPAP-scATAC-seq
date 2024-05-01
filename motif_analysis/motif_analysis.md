# Single-Cell ATAC-seq Motif Analysis

This directory contains scripts for motif analysis of single-cell ATAC-seq data, using the chromVAR tool and related packages.

## Scripts Overview

1. **01_chromVAR_scATAC_Analysis.R**:
   - Description: This script runs chromVAR on peaks derived from scATAC-seq data stored in a Seurat object. It includes data preparation, motif analysis, and variability & deviation analysis.
   - Inputs:
     - Seurat object with peaks.
     - Motif information 1. data/jaspar_2022_object.Rdata
   - Outputs:
     - Deviation scores, variability metrics, and motifdata for downstream.

2. **02_motif_Analysis_celltype.R**:
   - Description: Analyzes transcription factor binding variability across different cell types using chromVAR, incorporating data manipulation, motif analysis.
   - Inputs:
     - Deviation scores, motif data, and sample metadata.
   - Outputs:
     - Filtered and aggregated data for visualization of variability of motifs across cell types (in tabular format & boxplot) .

3. **03_motif_Analysis_disease.R**:
   - Description: Analyzes barcode-level covariates in single-cell ATAC-seq data to understand the variability in transcription factor binding sites across different conditions(Linear Mixed Model).
   - Inputs:
     - Sample metadata and deviation scores.
   - Outputs:
     - Results from mixed-effect linear models, assessing the impact of covariates on motif enrichment.

## Prerequisites

Ensure you have R installed on your machine. The scripts are compatible with R version 4.0 and above.

## Installation

Before running the scripts, you need to install the required R packages. Run the `install_packages.R` script located in the repository:

```bash
Rscript install_packages.R
