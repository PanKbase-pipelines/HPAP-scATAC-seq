# Load necessary libraries
library(patchwork)
library(readr)
library(stringr)
library(dplyr)
library(parallel)
library(ggpubr)
library(qvalue)
library(pheatmap)
library(RColorBrewer)
library(beeswarm)
library(lmerTest)
library(lme4)

# Ensure reproducibility
set.seed(1234)

# Load ATAC-seq data and metadata
atac <- readRDS('HPAP_atac_obj_withUnifiedpeaks.rds')
meta <- atac@meta.data

# Specify directory for data inputs and outputs
in_dir <- '../data'

# Export ATAC metadata as TSV
write.table(meta, file=paste0(in_dir, 'atac.metadata.tsv'), quote = FALSE, sep='\t')

# Read in barcode level covariates and additional metadata
info <- read.table(file=paste0(in_dir, 'info_table.txt'))
meta <- read.table(file=paste0(in_dir, 'atac.metadata.tsv'), sep='\t', header=TRUE)
wd_meta <- in_dir
meta.data <- read.csv(paste0(wd_meta, 'HPAP-scATAC-metadata-qc.csv'))
disease <- meta.data[, c("Sample_Name", "Disease", "AAB_Status", "Age", "Gender", "BMI", "Sample_Source")]

# Merge sample and barcode level covariates
disease.meta <- left_join(select(info, Sample_Name=def, groups, samples), disease)
disease.meta <- select(disease.meta, -Sample_Name, -Disease)
disease.meta <- disease.meta[!duplicated(disease.meta),]

# Load deviation scores
dev_file <- paste0(in_dir, 'devscores.txt')
devscores <- vroom::vroom(file=dev_file, skip=1, col_names=FALSE)
devscores <- tibble::column_to_rownames(devscores, var="X1")
colnames(devscores) <- str_split(readLines(file(dev_file),n=1), " ")[[1]]

# Filter devscores by a specific cell type (e.g., Q_Stellate)
cell_type <- 'Q_Stellate'
ct_info <- info[info$cluster == cell_type, ]
cell_type_barcodes <- rownames(ct_info)
ct_devscores <- devscores[cell_type_barcodes]
ct_meta <- tibble::rownames_to_column(select(filter(meta, Cell.Type == cell_type), TSS.enrichment, frac_reads_in_peaks, nCount_ATAC_peaks, nFeature_ATAC_peaks), var='bc')

# Join disease state info and sample info
barcode_motif_matrix <- left_join(tibble::rownames_to_column(as.data.frame(t(ct_devscores)), var="barcode"), 
                                  select(ct_info, groups, samples, barcode=cells))
barcode_motif_matrix <- select(barcode_motif_matrix, bc=barcode, everything())

# Prepare the data for mixed-effect linear modeling
new_df <- barcode_motif_matrix_bc_joined
new_df$samples <- as.character(new_df$samples) # Sample names were numbers; the formula tried to treat this as a continuous variable
new_df$Gender <- as.numeric(as.factor(new_df$Gender)) # Binarize, skip if not in model
new_df <- dplyr::filter(new_df, samples %in% names(which(table(barcode_motif_matrix$samples)  > 10))) # Filter for sufficient barcodes
new_df$T2D_weighted_barcode <- as.integer(new_df$groups=='T2D')
new_df$T2D_weighted_barcode[new_df$groups=='Control'] <- -(sum(new_df$groups=='T2D')/sum(new_df$groups=='Control'))

new_df$T1D_weighted_barcode <- as.integer(new_df$groups=='T1D')
new_df$T1D_weighted_barcode[new_df$groups=='Control'] <- -(sum(new_df$groups=='T1D')/sum(new_df$groups=='Control'))
new_df <- arrange(new_df, groups) %>%
  select(T2D_weighted_barcode, T1D_weighted_barcode, groups, samples, 
         Age, Gender, BMI, TSS.enrichment, frac_reads_in_peaks, nCount_ATAC_peaks, nFeature_ATAC_peaks, everything())

# Mixed-effect linear modeling to account for pseudoreplication
lmer_results <- data.frame(matrix(nrow=nrow(ct_devscores), ncol=7))
rownames(lmer_results) <- colnames(as.data.frame(t(ct_devscores)))
colnames(lmer_results) <- c('Motif', 'lmer.T2D.Effect', 'lmer.T2D.p.value', 'lmer.T2D.t.stat', 'lmer.T1D.Effect', 'lmer.T1D.p.value', 'lmer.T1D.t.stat')

# Iterate over motifs to fit models and store results
lmer_results = data.frame(matrix(nrow=nrow(ct_devscores), ncol = 7))
rownames(lmer_results)= colnames(as.data.frame(t(ct_devscores)))
# Create columns to save outputs. I took the effect estimate, the p value, and the t value (which is the test statistic)
colnames(lmer_results) <- c('Motif', 'lmer.T2D.Effect', 'lmer.T2D.p.value', 'lmer.T2D.t.stat',
                                          'lmer.T1D.Effect', 'lmer.T1D.p.value', 'lmer.T1D.t.stat')

for(i in 13:(ncol(new_df) - 2)){
#for(i in 17:23){
    # Get the motif name
    column <- names(new_df[i])

    # Run the model. Adjust with covariates, it is a good idea to scale and center them.
    # The random effect for sample is encoded as (1|sample), which indicates it allows the 
    # intercept but not the slope to vary
    mod6=lmerTest::lmer(new_df[,i] ~ T2D_weighted_barcode + T1D_weighted_barcode + scale(frac_reads_in_peaks) + scale(nCount_ATAC_peaks) + (1|samples), data = new_df)

    # Save outputs
    lmer_results[column,'Motif'] <- column

    # Check that you are getting the right rows of the model
    lmer_results[column,'lmer.T2D.Effect'] <- summary(mod6)[['coefficients']][2,1]
    lmer_results[column,'lmer.T2D.p.value'] <- summary(mod6)[['coefficients']][2,5]
    lmer_results[column,'lmer.T2D.t.stat'] <- summary(mod6)[['coefficients']][2,4]

    lmer_results[column,'lmer.T1D.Effect'] <- summary(mod6)[['coefficients']][3,1]
    lmer_results[column,'lmer.T1D.p.value'] <- summary(mod6)[['coefficients']][3,5]
    lmer_results[column,'lmer.T1D.t.stat'] <- summary(mod6)[['coefficients']][3,4]
}

# FDR calculations
lmer_results <- arrange(lmer_results, lmer.T2D.p.value)
lmer_results$lmer.T2D.q.value <- qvalue(lmer_results$lmer.T2D.p.value)$qvalues

lmer_results <- arrange(lmer_results, lmer.T1D.p.value)
lmer_results$lmer.T1D.q.value <- qvalue(lmer_results$lmer.T1D.p.value)$qvalues

# Save results and export for further analysis
write.table(lmer_results, paste0(in_dir, 'Q_Stellate_mixed_model_outs.tsv'), sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(new_df, paste0(in_dir, 'Q_Stellate_mixed_model_outs_boxplot.tsv'), sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
