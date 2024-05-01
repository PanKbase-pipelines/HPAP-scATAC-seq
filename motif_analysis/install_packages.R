# List of packages required by the scripts
packages <- c(
  "Seurat", "Signac", "JASPAR2020", "TFBSTools", "BSgenome.Hsapiens.UCSC.hg38",
  "patchwork", "readr", "stringr", "dplyr", "Matrix", "BiocParallel", "parallel",
  "plyr", "ggpubr", "car", "qvalue", "pheatmap", "RColorBrewer", "beeswarm",
  "chromVAR", "motifmatchr", "SummarizedExperiment", "lme4", "lmerTest", "vroom", "tibble"
)

# Function to install missing packages
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
    if (!require(package, character.only = TRUE, quietly = TRUE)) {
      message(paste("Failed to install", package))
    }
  }
}

# Apply the function to each package
sapply(packages, install_if_missing)

# Also, ensure Bioconductor packages are installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("JASPAR2020", "TFBSTools", "BSgenome.Hsapiens.UCSC.hg38", "motifmatchr", "SummarizedExperiment"))