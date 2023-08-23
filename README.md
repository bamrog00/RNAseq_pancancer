# RNAseq_pancancer
Run the RNAseq analysis on all cohorts with data preparation. The analysis includes differential expression analysis, gene set enrichment analysis and score creation with genes from pathways.


## Python

### dea_data_preparation
Merges the HRD results with the sample sheet downloaded from GDC. Creates the countMatrix and the colData for the differential expression analysis. Additionally uses the cutoff from either GMM or K-means to add a column to the colData defining the type.

## R

### RNAseq_analysis_pan_cancer
Contains the pipeline that runs the RNAseq analysis on each cohort (or a single cohort). Creates dds object for the n2 and n3 case (if it is possible), creates volcano plots, heatmaps, runs GSEA with KEGG and Hallmark pathways. Creates new scores by first z-scoring the genes and summing up all genes from specific sets (Pathways) for each sample.

## How to use
Data needs to prepared so that per cancer types there are two dataframes countMatrix and colData using the dea_data_preparation notebook.
Using the functions in RNAseq_analysis_pan_cancer, RNAseq analysis is conducted and figures and tables are saved. Main functions to use are runMultipleCohorts, runMultipleCoorts_SaveScores (saves the gene set scores in addition), and for Enrichment plots ( multipleDDS, multipleGSEA, rewriteCINSARC and plotEnrichment).
