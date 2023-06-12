# RNAseq_pancancer
Run the RNAseq analysis on all cohorts with data preparation. The analysis includes differential expression analysis, gene set enrichment analysis and score creation with genes from pathways.


## Python

### dea_data_preparation
Merges the HRD-results with the sample sheet downloaded from GDC. Creates the countMatrix and the colData fro the differential expression analysis. Addtionally uses the cutoff from either GMM or K-means to add a column to the colData defining the type.

## R

### RNAseq_analysis_pan_cancer
Contains the pipeline which runs the RNAseq analysis on each cohort (or a single cohort). Creates dds object for the n2 and n3 case (if it is possible), creates volcano plots, heatmaps, runs GSEA with KEGG and Hallmark pathways. Creates new scores by first z-scoring the genes and summing up all genes from specific sets (Pathways) for each sample.
Folder structur needs to be: TODO

