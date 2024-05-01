# Load necessary libraries
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(readr)
library(stringr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(beeswarm)
library(parallel)
library(ggpubr)
library(car)
library(qvalue)

suppressMessages({
  library(chromVAR)
  library(motifmatchr)
  library(SummarizedExperiment)
})

# Initial setup
set.seed(1234)  # Ensure reproducibility

# Data Loading and JASPAR database
in_dir <- '../data/'
motifdata <- read.table(file=paste0(in_dir, 'motifdata.txt'), sep="\t")
info <- read.table(file=paste0(in_dir, 'info_table.txt'), sep="\t")
TFClass_Lookup <- read_csv("../data/Chromvar_to_Gene_By_Subfam_Complete_JAPRAR2022_TFClass.csv")
TFClass_Full <- read_csv("../data/Chromvar_to_Gene_Jaspar2022.csv")
variability <- read.table(file=paste0(in_dir, 'variability.txt'), sep="\t")

# Deviation score loading using vroom for faster I/O
dev_file <- paste0(in_dir, 'devscores.txt')
devscores <- vroom::vroom(file=dev_file, skip=1, col_names=FALSE)
devscores <- tibble::column_to_rownames(devscores, var="X1")
colnames(devscores) <- str_split(readLines(file(dev_file), n=1), " ")[[1]]

# Filter by variability
variable_to_keep <- rownames(filter(variability, variability > 1.2))
var_subset_devscores <- devscores[variable_to_keep,]

# Compute averages by cell type
mat_type = sapply(unique(info$cluster), 
                function(i) rowMeans(var_subset_devscores[info$cells[info$cluster==i]], na.rm=TRUE))
mat_df <- as.data.frame(mat_type)
mat_df_high_dev <- mat_df[apply(abs(mat_df) > 5, 1, na.rm=TRUE, any),]

# Organize motif database
colnames(motifdata) <- c('Motif','Name','jaspar_family')
motifdata_fixed <- left_join(motifdata, select(TFClass_Full, Name=jaspar_name_1, TFClass_family=lowest_level_fam_1))
motifdata_fixed <- motifdata_fixed[!duplicated(motifdata_fixed),]

# Read and Join with disease metadata
info$groups <- sub("Control", "ND", info$groups)
wd_meta='../data/'
meta.data =read.csv(paste0(wd_meta,'HPAP-scATAC-metadata-qc.csv'))
disease = meta.data[,c("Sample_Name","Disease","Age","Gender","BMI")]
disease_remove <- left_join(disease, select(info[!duplicated(info$samples),], Sample_Name=def, sample=samples)) %>%
                    filter(!is.na(sample))
disease_remove$sample <- as.character(disease_remove$sample)
info$cluster_samples <- paste0(info$cluster, ":", info$samples)

persamp_pertype =sapply(unique(info$cluster_samples), 
                function(i) rowMeans(as.data.frame(var_subset_devscores[info$cells[info$cluster_samples==i]]), na.rm=TRUE))
persamp_pertype_transform <- t(persamp_pertype)
persamp_pertype_transform <- as.data.frame(persamp_pertype_transform)
persamp_pertype_transform$Cell_Type <- str_split(rownames(persamp_pertype_transform), ":", simplify=TRUE)[,1]
persamp_pertype_transform$sample <- str_split(rownames(persamp_pertype_transform), ":", simplify=TRUE)[,2]
  
motif_matrix <- dplyr::inner_join(persamp_pertype_transform, disease_remove, by = 'sample')

# Data Export
out_dir <- "../data/"
write.table(motif_matrix, file=paste0(in_dir, 'motif_matrix_boxplot.txt'), quote=FALSE, sep='\t')
write.csv(var_subset_devscores, paste0(out_dir, "variable_devscores.csv"))
write.csv(mat_df, paste0(out_dir, "celltype_average_deviations.csv"))
