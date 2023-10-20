# MethylSeq_Pipeline
This repository contains the scripts used to analyze Methylseq datasets across all the Groups using MethylKit package in R

## The Input required for the Methyl-Seq Analysis is the Output froM EMseq pipeline.

# The Singularity image is added into the directory: /storage/colddata/basesolve/tools/METHYLKIT_FINAL.sif

# The Script used to perform the Analysis is: /storage/colddata/basesolve/scripts/METHYLSEQ/METHYLSEQ_ANALYSIS.R

=> The input required for this Analysis is the MethylKit output files in all the 3 context added into separate directories.

1. CpG_OUTPUT
2. CHH_OUTPUT
3. CHG_OUTPUT

The Entire Directory needs to be given as input

=> The Command line used to perform this Analysis is:

singularity exec /storage/colddata/basesolve/tools/METHYLKIT_FINAL.sif Rscript /storage/colddata/basesolve/scripts/METHYLSEQ/METHYLSEQ_ANALYSIS.R -p RESULTS -b /storage/colddata/basesolve/databases/methylseq/GRCh38/genes.bed -i CpG_OUTPUT -m SAMPLE_METADATA_TABLE.txt
