# Load necessary libraries
suppressMessages({
  library(chromVAR)
  library(motifmatchr)
  library(SummarizedExperiment)
  library(Signac)
  library(Seurat)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(patchwork)
  library(readr)
  library(stringr)
  library(dplyr)
  library(Matrix)
  library(BiocParallel)
  library(parallel)
  library(plyr)
  library(ggpubr)
  library(car)
  library(qvalue)
})

# Set up environment
set.seed(1234)
register(MulticoreParam(8))

# Load Seurat object with unified peaks
atac <- readRDS('/path/to/your/HPAP_atac_obj_withUnifiedpeaks.rds')
DefaultAssay(atac) <- 'Unified_Peaks'

# Prepare inputs for SummarizedExperiment
# Extract matrix from Seurat object
sc.data <- GetAssayData(atac, slot='data', assay="Unified_Peaks")

# Extract peak locations and reformat into GRanges object
# This may need to be modified depending on your peak naming convention
bed <- str_split_fixed(rownames(sc.data), "\\-", 3)
bed[,1] <- paste0('chr', bed[,1])
gr <- GRanges(seqnames=bed[,1], ranges=IRanges(start=as.numeric(bed[,2]), end=as.numeric(bed[,3])))

# Create SummarizedExperiment object
fragment.counts <- SummarizedExperiment(assays=list(counts=sc.data), rowRanges=gr)

# Adjust to include any columns from the Seurat object metadata you will use downstream. Easier to include more now.
metrics <- select(atac[[]], orig.ident, nCount_ATAC_peaks, nFeature_ATAC_peaks, library, sex,
                  condition, seurat_clusters, Cell.Type) 
if (length(rownames(metrics)) == sum(rownames(colData(fragment.counts)) == rownames(metrics))) {
    print("Success adding meta data")
    colData(fragment.counts) <- cbind(colData(fragment.counts), metrics[rownames(colData(fragment.counts)),])
} else {
    print("Failed to add meta data, check column names.")
}

# Bias Correction and Motif Analysis
# Applies GC bias correction to the counts. Loads JASPAR motifs and matches them to the genomic ranges using motifmatchr, preparing for chromVAR analysis.
fragment.counts <- addGCBias(fragment.counts, genome=BSgenome.Hsapiens.UCSC.hg38)
fragment.counts

# Load and match motifs using JASPAR2020
jaspar.motifs <- readRDS(file ='../data/jaspar_2022_object.Rdata')
motif.ix <- matchMotifs(jaspar.motifs, fragment.counts, genome=BSgenome.Hsapiens.UCSC.hg38)

# Run chromVAR
dev <- computeDeviations(object=fragment.counts, annotations=motif.ix)
saveRDS(dev, "../data/ChromVAR_Object.RDS")

#This is a cell by motif deviation score (aka accessibility) matrix
devtab = deviationScores(dev)

# Variation of accessibiility across deviation scores
variability <- computeVariability(dev)

#Pull in some metadata for downstream scripts
info = data.frame(def   = atac$library,
                    cells = Cells(atac),
                    groups = atac$condition,
                    cluster = atac@active.ident, samples = atac$library)
write.table(info, file='../data/info_table_ctrlvdiabetes.txt',quote = FALSE, sep='\t')

#  create a two-column matrix where each row corresponds to a motif, the first column contains the motif names, and the second column contains their respective matrix classes
motifdata = cbind(sapply(jaspar.motifs, function(x) unlist(x@name)),sapply(jaspar.motifs, function(x) unlist(x@matrixClass )))                                                  

# Save results
output_dir <- '../data/'
write.table(devtab, file=paste0(output_dir,'devscores.txt'),quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(variability, file=paste0(output_dir, 'variability.txt'), quote = FALSE, col.names = TRUE, row.names = TRUE, sep='\t')
write.table(motifdata, file=paste0(output_dir,'motifdata.txt'),quote = FALSE, col.names = FALSE,sep='\t')
write.table(info, file=paste0(output_dir,'info_table.txt'),quote = FALSE, sep='\t')

# Print summary information
print("Analysis complete. Check output directory for results.")
