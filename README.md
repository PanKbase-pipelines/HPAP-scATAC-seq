# HPAP-scATAC-seq

### Sample bash script to run CellRanger on snATAC-seq data:
for SAMPLE in HPAP-109; do ~/scripts/cellranger-atac-2.0.0/cellranger-atac count --id ${SAMPLE} --fastqs ~/hpap/atac/${SAMPLE}/Upenn_scATACseq/fastq/ --sample ${SAMPLE} --reference ~/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/ --localcores 24 --disable-ui
