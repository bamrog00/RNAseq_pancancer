library(DESeq2)
library(ggplot2)
library(kohonen)
library(tsne)
library(umap)
library(pheatmap)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(ComplexHeatmap)
library(ggstance)
library(forcats)
library(esquisse)
library(dbscan)
library(factoextra)
library(colorRamp2)
library(gplots)
library(corrplot)
library(Rtsne)
library(ggpubr)
library(tidyr)

## Does not work as there are some gene names duplicates
#count_matrix = read.csv('../data/counts_matrix.csv', sep = ',', row.names = 1)

# Load data and handle duplicate gene names



hrdgenes<-read.gmt("../data/HRDgenes.gmt")
kegg_symbols = read.gmt('../data/c2.cp.kegg.v7.5.1.symbols.gmt')
hallmark = read.gmt('../data/h.all.v7.5.1.symbols.gmt')
CINSARC = readxl::read_xlsx('../data/Genes-CINSARC.xlsx')

MMEJ = read.gmt('../data/REACTOME_HDR_THROUGH_MMEJ_ALT_NHEJ.v2023.1.Hs.gmt')
SSA = read.gmt('../data/REACTOME_HDR_THROUGH_SINGLE_STRAND_ANNEALING_SSA.v2023.1.Hs.gmt')
HRR = read.gmt('../data/REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR.v2023.1.Hs.gmt')




##### CHANGE pipline to exclude the 200 genes score ####### Done
##### Check nbew methods again for parameters


#### Functions ####

# Creates Volcano plots 
volcano_plotter = function(res, hrd_genes, reference_type, compare_type, mode, pathway){
  ## Creates a volcane plot of the DEA results highlighting the significant HRD genes
  ## Input:
  ## res: results from the DEA, optained by using the function result(dds)
  ## hrd_genes (dataframe): contains the pathway name in one column and in the column gene the genes
  ## reference_type (string): Name of the refernce type
  ## compare_type (string): Name of the 'other' type
  ## mode (string; default is 'n2'): Indicates if only HRD-High and HRD-Low (mode = n2) or HRD-High,-Medium and -Low (mode ='n3) are compared
  ## pathway (string): pathway to the folder where the plots should be saved
  
  res_df = as.data.frame(res)
  res_df$symbol = rownames(res_df)
  
  
  ### Extract the significant genes
  hrdgenes_all <- subset(res_df, rownames(res_df) %in% hrdgenes$gene)
  hrdgenes_all <- subset(hrdgenes_all, hrdgenes_all$log2FoldChange<0.58 | hrdgenes_all$log2FoldChange> -0.58 )
  hrdgenes_significant <- subset(res_df, rownames(res_df) %in% hrdgenes$gene & res_df$padj<=0.05 )
  hrdgenes_significant <- subset(hrdgenes_significant, hrdgenes_significant$log2FoldChange>=0.58 | hrdgenes_significant$log2FoldChange<= -0.58 )
  sig_genes_label<-subset(hrdgenes_significant)
  
  ### Create the plot
  volcano_plot = ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(color="grey87") +
    ggtitle("Genes belonging to the HR pathway") +
    theme_bw() +
    geom_point(data= hrdgenes_all, col="dodgerblue2") +
    geom_point(data= hrdgenes_significant, col="red") +
    geom_text_repel(data = sig_genes_label,
                    aes(x=log2FoldChange,
                        y=-log10(padj),label= symbol),
                    max.overlaps = 40) +
    theme(legend.position = "none") +
    scale_x_continuous(name = paste("log2(fold change) ",compare_type," vs ", reference_type, sep = '')) +
    scale_y_continuous(name = "-log10 p-value") +
    geom_hline(yintercept = -log10(0.05), linetype="dashed") +
    geom_vline(xintercept = 0.58, linetype="dashed") +
    geom_vline(xintercept = -0.58, linetype="dashed")
  
  #print(volcano_plot)
  
  if (mode == 'n2'){
    ggsave(paste(pathway,'volcanoplot_',compare_type,'_',reference_type,'_n2.png',sep = ''),volcano_plot, width = 8.63, height = 5.71)
  }else{
    ggsave(paste(pathway,'volcanoplot_',compare_type,'_',reference_type,'_n3.png',sep = ''),volcano_plot, width = 8.63, height = 5.71)
    
  }
}

## Kegg and Hallmark GSEA

kegg_hallmark_hrdgenes_gsea = function(res, hrdgenes, kegg_symbols, hallmark, reference_type, compare_type, mode, pathway){
  ##3 TODO: delete hrd GSEA
  ## Runs the GSEA with the KEGG and Hallmark pathways and creates plots of the results
  ## Input:
  ## res: results from the DEA, optained by using the function result(dds)
  ## hrd_genes (dataframe): contains the pathway name in one column and in the column gene the genes
  ## kegg_symbols (dataframe): contains the pathway name in one column and in the column gene the genes
  ## hallmark (dataframe): contains the pathway name in one column and in the column gene the genes
  ## reference_type (string): Name of the refernce type
  ## compare_type (string): Name of the 'other' type
  ## mode (string; default is 'n2'): Indicates if only HRD-High and HRD-Low (mode = n2) or HRD-High,-Medium and -Low (mode ='n3) are compared
  ## pathway (string): pathway to the folder where the plots should be saved
  ## Output:
  ## result[1] (dataframe): Results for the kegg pathway
  ## results[2] (dataframe): Results from the hallmark pathway
  
  gene_list = res$stat
  res$symbol = rownames(res)
  names(gene_list) = make.names(res$symbol, unique = T)
  gene_list = gene_list[order(gene_list, decreasing = T)]
  
  gene_names = make.names(res$symbol, unique = T)
  
  ### Running the GSEA for hrdgenes, kegg and hallmark pathway
  
  hrdgenes_gsea = GSEA(gene_list, TERM2GENE = hrdgenes,
                       pvalueCutoff = 0.1,
                       eps=1e-50,
                       seed=T)
  
  kegg_gsea = GSEA(gene_list, TERM2GENE = kegg_symbols,
                   pvalueCutoff = 0.1,
                   eps=1e-50,
                   seed=T)
  
  hallmark_gsea = GSEA(gene_list, TERM2GENE = hallmark,
                       pvalueCutoff = 0.1,
                       eps=1e-50,
                       seed=T)
  
  kegg_gsea_df = as.data.frame(kegg_gsea)
  hrdgenes_gsea_df = as.data.frame(hrdgenes_gsea)
  hallmark_gsea_df = as.data.frame(hallmark_gsea)
  
  
  ### Creating the barplots
  y <- mutate(kegg_gsea, ordering = abs(NES)) %>%
    arrange(desc(ordering))
  n <- 15
  y_bar <- group_by(y, sign(NES)) %>%
    slice(1:n)
  kegg_pathways = ggplot(y_bar, aes(NES, fct_reorder(Description, NES), fill = qvalue), showCategory=(n*2)) +
    geom_barh(stat='identity') +
    scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) +
    theme_minimal() + ylab(NULL) + ggtitle(paste("KEGG GSEA cell models, ", compare_type,' vs ', reference_type, sep = ''))+
    theme(plot.title = element_text(hjust = 1.0), panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"))
  
  y <- mutate(hallmark_gsea, ordering = abs(NES)) %>%
    arrange(desc(ordering))
  n <- 15
  y_bar <- group_by(y, sign(NES)) %>%
    slice(1:n)
  hallmark_pathways = ggplot(y_bar, aes(NES, fct_reorder(Description, NES), fill = qvalue), showCategory=(n*2)) +
    geom_barh(stat='identity') +
    scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) +
    theme_minimal() + ylab(NULL) + ggtitle(paste("Hallmark GSEA cell models, ", compare_type,' vs ', reference_type, sep = ''))+
    theme(plot.title = element_text(hjust = 1.0), panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"))
  
  #print(kegg_pathways)
  #print(hallmark_pathways)
  
  if (mode == 'n2'){
    ggsave(paste(pathway,'kegg_pathways_',compare_type,'_vs_',reference_type,'_n2.png',sep = ''),kegg_pathways, width = 8.63, height = 5.71)
    ggsave(paste(pathway,'hallmark_pathways_',compare_type,'_vs_',reference_type,'_n2.png',sep = ''),hallmark_pathways, width = 8.63, height = 5.71)
    
  }else{
    ggsave(paste(pathway,'kegg_pathways_',compare_type,'_vs_',reference_type,'_n3.png',sep = ''),kegg_pathways, width = 8.63, height = 5.71)
    ggsave(paste(pathway,'hallmark_pathways_',compare_type,'_vs_',reference_type,'_n3.png',sep = ''),hallmark_pathways, width = 8.63, height = 5.71)
    
  }
  
  return (list(hrdgenes_gsea = hrdgenes_gsea_df, kegg_gsea = kegg_gsea_df, hallmark_gsea = hallmark_gsea_df))
}

heatmap_DEA = function(dds, res, genes, mode, cancertype, save_heat = FALSE){
  ## creates a heatmap of the DEA results
  ## Input:
  ## dds (dds object)
  ## res: results from the DEA, optained by using the function result(dds)
  ## genes (list of strings): List of the genes that should be displayed on the heatmap
  ## mode (string; default is 'n2'): Indicates if only HRD-High and HRD-Low (mode = n2) or HRD-High,-Medium and -Low (mode ='n3) are compared
  ## cancertype (string): The cancertype has to be without the prefix like 'TCGA' or 'TARGET'
  ## save_heat (boolean): Indicates if the heatmap should be saved (TRUE) or not (FALSE). (default is not saving it)
  
  
  
  ## Constructing the dataframe for the heatmap
  res_df = as.data.frame(res)
  res_df = res_df[order(res_df$padj),]
  
  vst <- vst(dds, blind=FALSE)
  vst <- assay(vst)
  vst <- as.data.frame(vst)
  vst_sig <- vst
  tvst <- t(vst_sig)
  rownames(tvst) <- NULL
  
  heat <- t(scale(tvst))
  
  ### Creating the annotation, depeding on the mode
  if (mode == 'n2'){
    colnames(heat) = paste0(dds@colData$gmm_n2_cut)
    ha = HeatmapAnnotation(HRDtype = colnames(heat), col = list(HRDtype = c('High' = 'darkblue', 'Low' = 'lightblue')))
  }
  if (mode == 'n3'){
    colnames(heat) = paste0(dds@colData$gmm_n3_cut)
    ha = HeatmapAnnotation(HRDtype = colnames(heat), col = list(HRDtype = c('High' = 'darkblue', 'Medium' = 'green','Low' = 'orange')))
  }
  
  heat = heat[rownames(heat) %in% genes,]
  
  
  ### Get the HRD scores
  hrds_values = dds@colData$HRD_sum
  loh_values = dds@colData$LOH
  lst_values = dds@colData$LST
  tai_values = dds@colData$TAI
  
  hrd_annotation = HeatmapAnnotation(HRDs = hrds_values, col = list(HRDs = colorRamp2(c(min(hrds_values),max(hrds_values)),c('white','red'))))
  loh_annotation = HeatmapAnnotation(LOH = loh_values, col = list(LOH = colorRamp2(c(min(loh_values),max(loh_values)),c('lightgreen', 'darkgreen'))))
  lst_annotation = HeatmapAnnotation(LST = lst_values, col = list(LST = colorRamp2(c(min(lst_values),max(lst_values)),c('lightblue','darkblue'))))
  tai_annotation = HeatmapAnnotation(TAI = tai_values, col = list(TAI = colorRamp2(c(min(tai_values),max(tai_values)),c('lightcyan', 'darkcyan'))))
  
  if (save_heat){
    png(file=paste('../RNAseq_',cancertype,'/figures/heatmap_',mode, '.png', sep = ''))
  }
  
  hm = Heatmap(heat,show_column_names = FALSE, row_names_gp = gpar(fontsize = 4), bottom_annotation = ha, column_km = 4, name = 'z-score') %v% 
    hrd_annotation %v% loh_annotation %v% lst_annotation %v% tai_annotation
  
  draw(hm)
  
  if (save_heat){
    dev.off()
  }
  
}



runDeSeq = function(countData, colData, mode = 'n2'){
  ## Creates a DeSEQ2 object and runs the DEA
  ## Input:
  ## countData (dataframe): Count matrix where the rows are genes and the columns are samples. It contains the raw counts
  ## colData (dataframe): Dataframe containing information about the samples, the rows have to match the columns of countData
  ## mode (string; default is 'n2'): Indicates if only HRD-High and HRD-Low (mode = n2) or HRD-High,-Medium and -Low (mode ='n3) are compared
  ## Output:
  ## Returns the dds object
  
  if (mode == 'n2'){
    dds = DESeqDataSetFromMatrix(countData = countData,
                                     colData = colData,
                                     design= ~ gmm_n2_cut)
    
    ### Prefiltering
    keep = rowSums(counts(dds)) >= 10
    dds = dds[keep,]
    dds$gmm_n2_cut = relevel(dds$gmm_n2_cut, ref = 'Low')
  }
  else{
    dds = DESeqDataSetFromMatrix(countData = countData,
                                 colData = colData,
                                 design= ~ gmm_n3_cut)
    
    ### Prefiltering
    keep = rowSums(counts(dds)) >= 10
    dds = dds[keep,]
    #dds$gmm_n3_cut = relevel(dds$gmm_n3_cut, ref = 'Low')
  }
  
  dds = DESeq(dds)
  
  return(dds)

}

loadData = function(cancertype){
  ## Function loads the count matrix and the colData for the given cancertype
  ## Input:
  ## cancertype (string): The cancertype has to be without the prefix like 'TCGA' or 'TARGET'
  ## Output:
  ## returns the countmatrix (results[[1]]) containing the raw counts (rows are genes, columns are samples) 
  ## and the colData (results[[2]]) containing information about the samples (rows are the samples, columns are different features/information), the rows have to match the columns of the countmatrix
  
  
  count_matrix = read.csv(paste('../RNAseq_',cancertype,'/preparedData/counts_matrix.csv', sep = ''), sep = ',', row.names = NULL, check.names = FALSE)
  
  ### Make duplicate genenames unqiue
  cleaned_row_names <- make.unique(count_matrix[,1])
  rownames(count_matrix) <- cleaned_row_names
  count_matrix <- subset(count_matrix, select=-1)
  
  colData = read.csv(paste('../RNAseq_',cancertype,'/preparedData/colData.csv', sep = ''), sep = ',', row.names = 1)
  
  ### Extracting only the primary sample types
  colData_primary = colData[colData$Type == 'Primary',]
  count_matrix_primary = count_matrix[, colnames(count_matrix) %in% rownames(colData_primary)]
  return (list(count_matrix_primary,colData_primary))
}

saveGSEA = function(results, cancertype, resultInfo){
  ## Saves the results from the KEGG and Hallmark GSEA
  ## Input:
  ## results (list of dataframes): Contains the results for the GSEA of hrdgenes [1] (not saved), KEGG [2] and Hallmark [3] as dataframes
  ## cancertype (string): The cancertype has to be without the prefix like 'TCGA' or 'TARGET'
  ## resultsInfo (string): Info of the compared types (high, medium, low) and the mode (n2 or n3)
  
  res_hallmark = results[[3]]
  res_kegg = results[[2]]
  
  write.csv(res_hallmark, paste('../RNAseq_',cancertype,'/results/hallmark_',resultInfo,'.csv', sep = ''), row.names = TRUE)
  write.csv(res_kegg,paste('../RNAseq_',cancertype,'/results/kegg_',resultInfo,'.csv', sep = ''), row.names = TRUE)
}



runAnalysis = function(dds, cancertype, mode){
  ## Function that takes a dds object and creates volcano plots, runs GSEA, plots the results and saves them
  ## and creates a heatmap with the 20 most significant HRDgenes for the given cancertype and mode
  ## Input:
  ## dds (dds Object)
  ## cancertype (string): The cancertype has to be without the prefix like 'TCGA' or 'TARGET'
  ## mode (string; default is 'n2'): Indicates if only HRD-High and HRD-Low (mode = n2) or HRD-High,-Medium and -Low (mode ='n3) are compared
  
  
  if (mode == 'n2'){
    res = results(dds, contrast = c('gmm_n2_cut', 'High', 'Low'), alpha = 0.05)
    
    ### Create volcano plot
    volcano_plotter(res, hrd_genes = hrd_genes, paste(cancertype,'-Low',sep=''), paste(cancertype,'-High',sep=''), mode = mode, paste('../RNAseq_',cancertype,'/figures/', sep = ''))
    
    
    ### Run GSEA and save it
    res_LowHigh = kegg_hallmark_hrdgenes_gsea(res, hrdgenes, kegg_symbols, hallmark, paste(cancertype,'-Low',sep=''), paste(cancertype,'-High',sep=''), mode = mode, paste('../RNAseq_',cancertype,'/figures/', sep = ''))
    
    saveGSEA(res_LowHigh, cancertype, 'low_high_n2')
    
    res_df = as.data.frame(res)
    res_df = res_df[rownames(res_df) %in% hrdgenes$gene,]
    res_df = res_df[order(res_df$padj),]
    
    genes_sig = rownames(res_df)[1:20]
    
    ### Create the heatmap
    heatmap_DEA(dds, res, genes_sig, mode = mode, cancertype, save_heat = TRUE)
    
    
  }else{
    ### Extract all three possible combination of comparision
    res_hm = results(dds, contrast = c('gmm_n3_cut', 'High', 'Medium'), alpha = 0.05)
    res_ml = results(dds, contrast = c('gmm_n3_cut', 'Medium', 'Low'), alpha = 0.05)
    res_hl = results(dds, contrast = c('gmm_n3_cut', 'High', 'Low'), alpha = 0.05)
    
    ### Create the volcano plots
    volcano_plotter(res_hl, hrd_genes = hrd_genes, paste(cancertype,'-Low',sep=''), paste(cancertype,'-High',sep=''), mode = mode, paste('../RNAseq_',cancertype,'/figures/', sep = ''))
    volcano_plotter(res_hm, hrd_genes = hrd_genes, paste(cancertype,'-Medium',sep=''), paste(cancertype,'-High',sep=''), mode = mode, paste('../RNAseq_',cancertype,'/figures/', sep = ''))
    volcano_plotter(res_ml, hrd_genes = hrd_genes, paste(cancertype,'-Low',sep=''), paste(cancertype,'-Medium',sep=''), mode = mode, paste('../RNAseq_',cancertype,'/figures/', sep = ''))
    
    
    ## Run and save the GSEA
    res_LowHigh = kegg_hallmark_hrdgenes_gsea(res_hl, hrdgenes, kegg_symbols, hallmark, paste(cancertype,'-Low',sep=''), paste(cancertype,'-High',sep=''),  mode = mode, paste('../RNAseq_',cancertype,'/figures/', sep = ''))
    res_MediumHigh = kegg_hallmark_hrdgenes_gsea(res_hm, hrdgenes, kegg_symbols, hallmark, paste(cancertype,'-Medium',sep=''), paste(cancertype,'-High',sep=''),  mode = mode, paste('../RNAseq_',cancertype,'/figures/', sep = ''))
    res_LowMedium = kegg_hallmark_hrdgenes_gsea(res_ml, hrdgenes, kegg_symbols, hallmark, paste(cancertype,'-Low',sep=''), paste(cancertype,'-Medium',sep=''),  mode = mode, paste('../RNAseq_',cancertype,'/figures/', sep = ''))
    
    
    saveGSEA(res_LowHigh, cancertype, 'low_high_n3')
    saveGSEA(res_MediumHigh, cancertype, 'medium_high_n3')
    saveGSEA(res_LowMedium, cancertype, 'low_medium_n3')
    
    res_df = as.data.frame(res_hl)
    res_df = res_df[rownames(res_df) %in% hrdgenes$gene,]
    res_df = res_df[order(res_df$padj),]
    
    genes_sig = rownames(res_df)[1:20]
    
    ### Run the heatmap only on the HRD High versus low
    heatmap_DEA(dds, res_hl, genes_sig, mode = mode, cancertype, save_heat = TRUE)
  }
  
}

getScores = function(dds, colData){
  ## Function which calculates the different scores given the different gene sets and adds them to the colData.
  ## The scores are calculated by first z-scoring each gene and then summing up all z-scores of one sample.
  ## Currently the gene sets are: HRR_70, CINSARC, MMEJ, SSA, HRRv2
  ## Input:
  ## dds (dds object): Contains the normalized counts
  ## colData (dataframe): Dataframe containg information per sample (rows are samples and the columns are features/information)
  ## Output:
  ## colData_copy (dataframe): Copy of the colData with the added scores as new columns
  
  colData_copy = as.data.frame(colData)
  
  norm_counts = counts(dds, normalized = TRUE)
  
  ### Extact only the genes of the genesets
  norm_counts_HRR = norm_counts[rownames(norm_counts) %in% hrdgenes$gene,]
  norm_counts_CINSARC = norm_counts[rownames(norm_counts) %in% CINSARC$Genes.CINSARC,]
  norm_counts_MMEJ = norm_counts[rownames(norm_counts) %in% MMEJ$gene,]
  norm_counts_SSA = norm_counts[rownames(norm_counts) %in% SSA$gene,]
  norm_counts_HRRv2 = norm_counts[rownames(norm_counts) %in% HRR$gene,]
  
  ### Z-score the genes
  z_scored_matrix_HRR = t(apply(norm_counts_HRR, 1, function(x) (x - mean(x)) / sd(x)))
  z_scored_matrix_CINSARC = t(apply(norm_counts_CINSARC, 1, function(x) (x - mean(x)) / sd(x)))
  z_scored_matrix_MMEJ = t(apply(norm_counts_MMEJ, 1, function(x) (x - mean(x)) / sd(x)))
  z_scored_matrix_SSA = t(apply(norm_counts_SSA, 1, function(x) (x - mean(x)) / sd(x)))
  z_scored_matrix_HRRv2 = t(apply(norm_counts_HRRv2, 1, function(x) (x - mean(x)) / sd(x)))
  
  ### Create a score per sample by summing up all genes
  sample_scores_HRR = colSums(z_scored_matrix_HRR)
  sample_scores_CINSARC = colSums(z_scored_matrix_CINSARC)
  sample_scores_MMEJ = colSums(z_scored_matrix_MMEJ)
  sample_scores_SSA = colSums(z_scored_matrix_SSA)
  sample_scores_HRRv2 = colSums(z_scored_matrix_HRRv2)

  ### Save the score in teh colData dataframe
  colData_copy$scoreHRR = sample_scores_HRR[rownames(colData_copy)]
  colData_copy$scoreCINSARC = sample_scores_CINSARC[rownames(colData_copy)]
  colData_copy$scoreMMEJ = sample_scores_MMEJ[rownames(colData_copy)]
  colData_copy$scoreSSA = sample_scores_SSA[rownames(colData_copy)]
  colData_copy$scoreHRRv2 = sample_scores_HRRv2[rownames(colData_copy)]

  return(colData_copy)
}



# Create two dataframes each with all the genes, adjust the rownames with an extanstion _HRD and add both together
# Somehow annotated the genes with a new column either HRD or other pathway and use row_split = df$pathway (try it out)

construct_subdf = function(data, geneset, genesetname){
  ## Constructs two dataframe, one contains the normalized counts only of the genes that are in the gene set
  ## and the second contains the original genenames, a new name which is composed of genename_genesetname and the genesetname.
  ## (Function is need for the function construct_df which constructs the dataframe for the heatmap score comparison)
  ## Input:
  ## data (dataframe): Contains the normalized counts of all genes (rows) for each sample (columns)
  ## geneset (list of strings): List containing the names of the genes in the geneset
  ## gensetname (string): Name of the geneset
  ## Output:
  ## results[1] (dataframe): Contains the normalized counts of the genes of the geneset
  ## results[2] (dataframe): Contains the genenames, new genenames and genesetname
  
  
  ### Create two new empty dataframes
  new_df = data.frame(matrix(nrow = 0, ncol = length(colnames(data))))
  colnames(new_df) = colnames(data)
  
  gene_df = data.frame(matrix(nrow = 0, ncol = 3))
  col_names_genes_df = c('original_genename','new_genenames','geneset')
  colnames(gene_df) = col_names_genes_df
  
  for (gene in geneset){
    if (gene %in% rownames(data)){
      
      ### Add the gene to the new dataframe if it is in the data and in the geneset
      new_df = rbind(new_df,data[gene, , drop=FALSE])
      
      ### AConstruct row for the gene dataframe
      new_gene_row = c(gene, paste(gene,'_',genesetname,sep = ''),genesetname)
      
      names(new_gene_row) = col_names_genes_df
      new_gene_row <- data.frame(t(new_gene_row), stringsAsFactors = FALSE)
      gene_df = rbind(gene_df, new_gene_row)
    }
  }
  
  rownames(new_df) = gene_df$new_genenames
  return(list(new_df, gene_df))
}

construct_df = function(data, geneset_1, geneset_1_name, geneset_2, geneset_2_name){
  ## Constructs two new dataframes containing the normalized counts for each gene in both geneset 
  ## and the second containing the original gene names, the new gene names and the geneset name.
  ## These two dataframe asre used for the heatmap which compares two genesets
  ## The first rows are the genes from geneset1 and the second part are the genes from geneset2
  ## Input:
  ## data (dataframe): Contains the normalized counts of all genes (rows) for each sample (columns)
  ## geneset_1 (list of strings): List containing the names of the genes in the geneset1
  ## gensetset_1_name (string): Name of the geneset1
  ## geneset_2 (list of strings): List containing the names of the genes in the geneset2
  ## gensetset_2_name (string): Name of the geneset2
  ## Output:
  ## results[1] (dataframe): Contains the normalized counts of the genes of the geneset1 and geneset2
  ## results[2] (dataframe): Contains the genenames, new genenames and genesetname of both geneset
  
  results_1 = construct_subdf(data, geneset_1, geneset_1_name)
  results_2 = construct_subdf(data, geneset_2, geneset_2_name)
  
  new_df = rbind(results_1[[1]],results_2[[1]])
  gene_df = rbind(results_1[[2]],results_2[[2]])
  
  return(list(new_df, gene_df))
}


heatmap_two_Scores = function(dds, colData, geneset1, geneset1_name, geneset2, geneset2_name, mode, annot_mode, cancertype, save_heat = FALSE){
  ## Fucntion to create a heatmap comparing two genesets
  ## Input:
  ## dds (dds object)
  ## colData (dataframe): Dataframe containg information per sample (rows are samples and the columns are features/information)
  ## geneset_1 (list of strings): List containing the names of the genes in the geneset1
  ## gensetset_1_name (string): Name of the geneset1
  ## geneset_2 (list of strings): List containing the names of the genes in the geneset2
  ## gensetset_2_name (string): Name of the geneset2
  ## mode (string; default is 'n2'): Indicates if only HRD-High and HRD-Low (mode = n2) or HRD-High,-Medium and -Low (mode ='n3) are compared
  ## annot_mode (string): Indicates if only high and low are annotated ('n2') or high, medium and low ('n3')
  ## cancertype (string): The cancertype has to be without the prefix like 'TCGA' or 'TARGET'
  ## save_heat (boolean): Indicates if the heatmap should be saved (TRUE) or not (FALSE). (default is not saving it)
  
  norm_counts = counts(dds, normalized = TRUE)
  
  ### Construct the dataframe for the heatmap
  results = construct_df(norm_counts, geneset1, geneset1_name, geneset2, geneset2_name)
  
  geneset_df = results[[2]]
  norm_counts_geneset = results[[1]]
  
  ### Z-score the genes
  z_scored_matrix_geneset = t(apply(norm_counts_geneset, 1, function(x) (x - mean(x)) / sd(x)))
  
  ### Get either type High and Low or High and Medium and Low
  if (mode == 'n2'){
    hrdtypes = colData[match(colnames(z_scored_matrix_geneset), rownames(colData)), "gmm_n2_cut"]
  }else{
    hrdtypes = colData[match(colnames(z_scored_matrix_geneset), rownames(colData)), "gmm_n3_cut"]
  }
  
  colnames(z_scored_matrix_geneset) = hrdtypes
  heat = z_scored_matrix_geneset
  
  ### Create the annotation depending on the annotation mode selected
  if (annot_mode == 'n2'){
    ha = HeatmapAnnotation(HRDtype = colnames(heat), col = list(HRDtype = c('High' = 'darkblue', 'Low' = 'lightblue')))
  }else{
    ha = HeatmapAnnotation(HRDtype = colnames(heat), col = list(HRDtype = c('High' = 'darkblue', 'Medium' = 'green','Low' = 'orange')))
  }
  
  if (save_heat){
    png(file=paste('../RNAseq_',cancertype,'/figures/heatmap_',geneset1_name,'_vs_',geneset2_name,'_',mode,'_annotated_',annot_mode, '.png', sep = ''))
  }
  
  hm = Heatmap(heat,row_split = factor(geneset_df$geneset, levels = c(geneset1_name, geneset2_name)),show_row_dend = FALSE,show_column_names = FALSE, row_names_gp = gpar(fontsize = 4), bottom_annotation = ha, column_km = 3, name = 'z-score')
  
  draw(hm)
  
  if (save_heat){
    dev.off()
    #ggsave(paste('../RNAseq_',cancertype,'/figures/heatmap_',geneset1_name,'_vs_',geneset2_name,'_',mode,'_annotated_',annot_mode, '.png', sep = ''),hm, width = 8.63, height = 5.71)
  }
  
  
}



createGeneSigScore = function(dds, colData, mode, min_genes, max_genes){
  ## Function which extracts the most significant genes and creates a new score.
  ## Returns the genes and the colData with the new score
  ## Input:
  ## dds (dds object)
  ## colData (dataframe): Dataframe containg information per sample (rows are samples and the columns are features/information)
  ## min_genes (int): Minimum number of genes needed for the score
  ## max_genes (int): Maximum number of genes used for the score
  ## Output
  ## results[1] (boolean): Indicates if a score was made
  ## results[2] (dataframe): Dataframe containg information per sample (rows are samples and the columns are features/information) with the new score
  ## results[3] (list of strings): List of the significant genes use for the score (only if results[1] is TRUE)
  
  enough_genes = FALSE
  
  if (mode == 'n2'){
    res = results(dds, contrast = c('gmm_n2_cut', 'High', 'Low'), alpha = 0.05)
  }else{
    res = results(dds, contrast = c('gmm_n3_cut', 'High', 'Low'), alpha = 0.05)
  }
  res_df = as.data.frame(res)
  
  ### Extract the most significant genes
  genes_significant = subset(res_df, res_df$padj<=0.05 )
  genes_significant = subset(genes_significant, genes_significant$log2FoldChange>=0.58 | genes_significant$log2FoldChange<= -0.58 )
  
  genes_significant = genes_significant[order(genes_significant$padj),]
  
  gene_sign = rownames(genes_significant)[1:max_genes] 
  
  if (length(gene_sign) < min_genes){
    return(list(enough_genes, colData))
  }
  
  enough_genes = TRUE
  
  ### Create the new score
  norm_counts = counts(dds, normalized = TRUE)
  
  norm_counts_geneSig = norm_counts[rownames(norm_counts) %in% gene_sign,]
  z_scored_matrix_geneSig = t(apply(norm_counts_geneSig, 1, function(x) (x - mean(x)) / sd(x)))
  sample_scores_geneSig = colSums(z_scored_matrix_geneSig)
  colData$scoregeneSig = sample_scores_geneSig[rownames(colData)]
  
  return(list(enough_genes, colData, gene_sign))
}


checkMode_n3 = function(colData){
  ## Function to check if there are at least one sample for each case in the n3 mode
  ## If not the n3 analyisis cannot be made
  ## Input:
  ## colData (dataframe): Dataframe containg information per sample (rows are samples and the columns are features/information)
  ## Output (boolean): Indicates if n3 mode analysis can be made
  
  n_high = nrow(colData[colData['gmm_n3_cut'] == 'High',])
  n_medium = nrow(colData[colData['gmm_n3_cut'] == 'Medium',])
  n_low = nrow(colData[colData['gmm_n3_cut'] == 'Low',])
  
  if (n_high == 0 | n_medium == 0 | n_low == 0){
    return (FALSE)
  }else
    return(TRUE)
}


runSingleCohort = function(projectID){
  ## Runs the full analyiss on the given cohort
  ## Creates dds object (n2 and n3 if possible), creates volcano plot, runs GSEA, creates heatmaps, creates scores and heatmpas for scores
  ## Input: 
  ## projectID (string): Cohortname with the project prefix like TCGA-LAUD
  
  cancertype = strsplit(projectID,'-')[[1]]
  if (length(cancertype) > 2){
    cancertype = paste(cancertype[2],cancertype[3],sep = '-')
  } else{
    cancertype = cancertype[2]
  }
  
  ### Get the data
  data = loadData(cancertype)
  print(paste('Loaded data from cohort ', projectID, sep = ''))
  count_matrix = data[[1]]
  colData = data[[2]]
  
  n3_possible = checkMode_n3(colData)
  
  ### create both dds object for n2 and n3
  dds_n2 = runDeSeq(count_matrix, colData, mode = 'n2')
  print(paste('Created DeSeq Object for mode n2 for cohort ', projectID, sep = ''))
  
  if (n3_possible){
    dds_n3 = runDeSeq(count_matrix, colData, mode = 'n3')
    print(paste('Created DeSeq Object for mode n3 for cohort ', projectID, sep = ''))
  }else{
    print(paste('One class (high,medium,low) had 0 number of samples. N3 analysis could not be done'))
  }

  
  
  
  ### Run the analysis on the n2 and n3 data
  runAnalysis(dds_n2, cancertype, mode = 'n2')
  print(paste('Finished analysis for mode n2 for cohort ', projectID, sep = ''))
  if (n3_possible){
    runAnalysis(dds_n3, cancertype, mode = 'n3')
    print(paste('Finished analysis for mode n3 for cohort ', projectID, sep = ''))
  }

  
  ### Calculate the scores
  colData_scores = getScores(dds_n2, colData)
  
  results = createHRDGeneSigScore(dds_n2, colData_scores, mode ='n2', cancertype, min_genes = 5, save_genes = TRUE)
  print(paste('Created scores for cohort ', projectID, sep = ''))
  if (results[[1]]){
    enough_geneSig = results[[1]]
    colData_scores = results[[2]]
    geneSig = results[[3]]
  }else{
    enough_geneSig = results[[1]]
    colData_scores = results[[2]]
    geneSig = list()
  }
  
  ### Plot the scores as boxplot, in heatmap compared to HRR_70 genes and all scores together
  plotScores_boxplot_v2(colData_scores, cancertype, enough_geneSig)
  print(paste('Plotted the score for cohort ', projectID, sep = ''))
  
  
  createAll_heatmap_twoScores(dds_n2, colData_scores, mode = 'n2', annot_mode ='n2', cancertype, enough_geneSig, geneSig)
  print(paste('Created heatmap for score comparison for cohort ', projectID, sep = ''))
  
  heatmap_allScores(colData_scores, cancertype, mode = 'n2', annot_mode = 'n2', enough_geneSig, save_heat = TRUE)
  print(paste('Created heatmap with all scores for cohort ', projectID, sep = ''))
  
  
  print(paste('Finished cohort ', projectID, sep = ''))
  #return(list(number_samples, number_samples_low_n2, number_samples_high_n2, number_samples_low_n3, number_samples_medium_n3, number_samples_high_n3))
}


runMultipleCohorts = function(cohortlist){
  ## Runs the analysis for all cohorts in the input list
  ## Input:
  ## cohortlist (list of strings): List containing the names of the cohorts with the project prefix like TCGA-LUAD
  
  for (cohort in cohortlist){
    runSingleCohort(cohort)
  }
}


getStatisitcs = function(cohortlist, gmm_results){
  ## Gets the numbers for the cohorts including: number of casees, number of HRD-low cases with the n2 cutoff, number of HRD-High cases with the n2 cutoff,
  ## the n2 cutoff, the same for HRD-low, HRD-Medium and HRD-High for the n3 cutoff and the two n3 cutoffs. Saves the results in a csv file
  ## Input:
  ## cohortlist (list of strings): List containing the names of the cohorts with the project prefix like TCGA-LUAD
  ## gmm_results (dataframe): Results from the GMM, containg the n2 and the the two n3 cutoffs for each cohort (rows)
  
  column_names = c('projectID', 'n_samples', 'n_HRD_low_n2', 'n_HRD_high_n2', 'n2_cutoff_lh','n_HRD_low_n3', 'n_HRD_medium_n3', 'n_HRD_high_n3','n3_cutoff_lm','n3_cutoff_mh')
  summary_df = data.frame(matrix(ncol = length(column_names), nrow = 0))
  colnames(summary_df) = column_names
  
  for (cohort in cohortlist){
    
    cancertype = strsplit(cohort,'-')[[1]]
    if (length(cancertype) > 2){
      cancertype = paste(cancertype[2],cancertype[3],sep = '-')
    } else{
      cancertype = cancertype[2]
    }
    
    
    data = loadData(cancertype)
    colData = data[[2]]
    
    number_samples = nrow(colData)
    number_samples_high_n2 = nrow(colData[colData['gmm_n2_cut'] == 'High',])
    number_samples_low_n2 = nrow(colData[colData['gmm_n2_cut'] == 'Low',])
    
    number_samples_high_n3 = nrow(colData[colData['gmm_n3_cut'] == 'High',])
    number_samples_low_n3 = nrow(colData[colData['gmm_n3_cut'] == 'Low',])
    number_samples_medium_n3 = nrow(colData[colData['gmm_n3_cut'] == 'Medium',])
    
    n2_cutoff_lh = gmm_results[gmm_results['Project.ID'] == cohort, 'n2_cutoff']
    n3_cutoff_lm = gmm_results[gmm_results['Project.ID'] == cohort, 'n3_cutoff_lm']
    n3_cutoff_mh = gmm_results[gmm_results['Project.ID'] == cohort, 'n3_cutoff_mh']
    
    new_row <- data.frame(projectID = cohort, n_samples = number_samples, n_HRD_low_n2 = number_samples_low_n2, n_HRD_high_n2 = number_samples_high_n2, n2_cutoff_lh = n2_cutoff_lh, n_HRD_low_n3 = number_samples_low_n3, n_HRD_medium_n3 = number_samples_medium_n3, n_HRD_high_n3 = number_samples_high_n3, n3_cutoff_lm = n3_cutoff_lm, n3_cutoff_mh = n3_cutoff_mh)    
    summary_df = rbind(summary_df, new_row)
    print(paste('Finished cohort ',cohort, sep = ''))
  }
  write.csv(summary_df, '../data/summary_number_samples.csv', row.names = FALSE)
}




createHRDGeneSigScore = function(dds, colData, mode, cancertype, min_genes, save_genes = FALSE){
  ## Function which extracts the most significant HRR genes and creates a new score.
  ## Returns the genes and the colData with the new score
  ## Input:
  ## dds (dds object)
  ## colData (dataframe): Dataframe containg information per sample (rows are samples and the columns are features/information)
  ## min_genes (int): Minimum number of genes needed for the score
  ## Output
  ## results[1] (boolean): Indicates if a score was made
  ## results[2] (dataframe): Dataframe containg information per sample (rows are samples and the columns are features/information) with the new score
  ## results[3] (list of strings): List of the significant HRR genes use for the score (only if results[1] is TRUE)
  
  enough_genes = FALSE
  
  if (mode == 'n2'){
    res = results(dds, contrast = c('gmm_n2_cut', 'High', 'Low'), alpha = 0.05)
  }else{
    res = results(dds, contrast = c('gmm_n3_cut', 'High', 'Low'), alpha = 0.05)
  }
  res_df = as.data.frame(res)
  res_df = subset(res_df, rownames(res_df) %in% hrdgenes$gene)
  
  ### Extract the most significant genes
  genes_significant = subset(res_df, res_df$padj<=0.05 )
  genes_significant = subset(genes_significant, genes_significant$log2FoldChange>=0.58 | genes_significant$log2FoldChange<= -0.58 )
  
  genes_significant = genes_significant[order(genes_significant$padj),]
  
  
  if (length(rownames(genes_significant)) < min_genes){
    return(list(enough_genes, colData))
  }
  
  enough_genes = TRUE
  
  ### Create the new score
  norm_counts = counts(dds, normalized = TRUE)
  
  norm_counts_geneSig = norm_counts[rownames(norm_counts) %in% rownames(genes_significant),]
  z_scored_matrix_geneSig = t(apply(norm_counts_geneSig, 1, function(x) (x - mean(x)) / sd(x)))
  sample_scores_geneSig = colSums(z_scored_matrix_geneSig)
  print(sample_scores_geneSig[rownames(colData)])
  colData$scoreHRDgeneSig = sample_scores_geneSig[rownames(colData)]
  
  if (save_genes){
    genes_df = as.data.frame(genes_significant)
    write.csv(genes_df, paste('../RNAseq_',cancertype,'/results/hrd_significant_genes_',mode,'.csv', sep = ''), row.names = FALSE)
  }
  
  return(list(enough_genes, colData, rownames(genes_significant)))
}








plotScores_boxplot = function(colData, cancertype, hrdgeneSigScore){
  ## Function that creates and saves boxplots for each score showing the distribution of the HRD-type (high, medium, low) for each score.
  ## Input:
  ## colData (dataframe): Dataframe containg information per sample (rows are samples and the columns are features/information)
  ## cancertype (string): The cancertype has to be without the prefix like 'TCGA' or 'TARGET'
  ## hrdgeneSigScore (boolean): If true it means that this score is available and a boxplot can be created
  
  
  comparisons = list( c('High', 'Low'), c('High', 'Medium'), c('Low', 'Medium'))
  
  
  ### create the different boxplots
  bp_corr_score_HRD_HRR = ggplot(data = colData, aes(x = gmm_n3_cut, y = scoreHRR, fill = gmm_n3_cut)) +
    geom_boxplot(show.legend = FALSE, notch = TRUE) +
    stat_boxplot(geom='errorbar', linetype=1, width=0.5)+
    xlab("HRD type") + stat_compare_means(comparisons = comparisons, method = 't.test')
  
  bp_corr_score_HRD_CINSARC = ggplot(data = colData, aes(x = gmm_n3_cut, y = scoreCINSARC, fill = gmm_n3_cut)) +
    geom_boxplot(show.legend = FALSE, notch = TRUE) +
    stat_boxplot(geom='errorbar', linetype=1, width=0.5)+
    xlab("HRD type") + stat_compare_means(comparisons = comparisons, method = 't.test')
  
  bp_corr_score_HRD_MMEJ = ggplot(data = colData, aes(x = gmm_n3_cut, y = scoreMMEJ, fill = gmm_n3_cut)) +
    geom_boxplot(show.legend = FALSE, notch = TRUE) +
    stat_boxplot(geom='errorbar', linetype=1, width=0.5)+
    xlab("HRD type") + stat_compare_means(comparisons = comparisons, method = 't.test')
  
  bp_corr_score_HRD_SSA = ggplot(data = colData, aes(x = gmm_n3_cut, y = scoreSSA, fill = gmm_n3_cut)) +
    geom_boxplot(show.legend = FALSE, notch = TRUE) +
    stat_boxplot(geom='errorbar', linetype=1, width=0.5)+
    xlab("HRD type") + stat_compare_means(comparisons = comparisons, method = 't.test')
  
  bp_corr_score_HRD_HRRv2 = ggplot(data = colData, aes(x = gmm_n3_cut, y = scoreHRRv2, fill = gmm_n3_cut)) +
    geom_boxplot(show.legend = FALSE, notch = TRUE) +
    stat_boxplot(geom='errorbar', linetype=1, width=0.5)+
    xlab("HRD type") + stat_compare_means(comparisons = comparisons, method = 't.test')
  
  
  
  #print(bp_corr_score_HRD_HRR)
  #print(bp_corr_score_HRD_CINSARC)
  #print(bp_corr_score_HRD_MMEJ)
  #print(bp_corr_score_HRD_SSA)
  #print(bp_corr_score_HRD_HRRv2)
  
  if (hrdgeneSigScore){
    bp_corr_score_HRD_geneSig = ggplot(data = colData, aes(x = gmm_n3_cut, y = scoreHRDgeneSig, fill = gmm_n3_cut)) +
      geom_boxplot(show.legend = FALSE, notch = TRUE) +
      stat_boxplot(geom='errorbar', linetype=1, width=0.5)+
      xlab("HRD type") + stat_compare_means(comparisons = comparisons, method = 't.test')
    #print(bp_corr_score_HRD_geneSig)
    ggsave(paste('../RNAseq_',cancertype,'/figures/',cancertype,'_Significant_genes_n3.png', sep = ''),bp_corr_score_HRD_geneSig,width = 8.63, height = 5.71)
    
  }
  
  ggsave(paste('../RNAseq_',cancertype,'/figures/',cancertype,'_MMEJ_genes_n3.png', sep = ''),bp_corr_score_HRD_MMEJ,width = 8.63, height = 5.71)
  ggsave(paste('../RNAseq_',cancertype,'/figures/',cancertype,'_SSA_genes_n3.png', sep = ''),bp_corr_score_HRD_SSA,width = 8.63, height = 5.71)
  ggsave(paste('../RNAseq_',cancertype,'/figures/',cancertype,'_HRRv2_genes_n3.png', sep = ''),bp_corr_score_HRD_HRRv2,width = 8.63, height = 5.71)
  ggsave(paste('../RNAseq_',cancertype,'/figures/',cancertype,'_HRR_genes_n3.png', sep = ''),bp_corr_score_HRD_HRR,width = 8.63, height = 5.71)
  ggsave(paste('../RNAseq_',cancertype,'/figures/',cancertype,'_CINSARC_genes_n3.png', sep = ''),bp_corr_score_HRD_CINSARC,width = 8.63, height = 5.71)
  
  
}


createAll_heatmap_twoScores = function(dds, colData, mode, annot_mode, cancertype, hrdgeneSigScore, hrdgeneSig = list()){
  ## Function to create heatmaps which compares two geneset. All the genesets are compared with HRR_70. If possible the geneSigScore as well.
  ## Input:
  ## dds (dds object)
  ## colData (dataframe): Dataframe containg information per sample (rows are samples and the columns are features/information)
  ## mode (string; default is 'n2'): Indicates if only HRD-High and HRD-Low (mode = n2) or HRD-High,-Medium and -Low (mode ='n3) are compared
  ## annot_mode (string): Indicates if only high and low are annotated ('n2') or high, medium and low ('n3')
  ## cancertype (string): The cancertype has to be without the prefix like 'TCGA' or 'TARGET'
  ## geneSigScore (boolean): If true it means that this score is available and a boxplot can be created
  ## geneSig (list of strings): The list of significant genes (only needed if geneSigScore is TRUE), (default is an empty list)
  
  
  heatmap_two_Scores(dds, colData, MMEJ$gene, 'MMEJ', hrdgenes$gene, 'HRR_70', mode, annot_mode, cancertype, save_heat = TRUE)
  heatmap_two_Scores(dds, colData, SSA$gene, 'SSA', hrdgenes$gene, 'HRR_70', mode, annot_mode, cancertype, save_heat = TRUE)
  heatmap_two_Scores(dds, colData, CINSARC$Genes.CINSARC, 'CINSARC', hrdgenes$gene, 'HRR_70', mode, annot_mode, cancertype, save_heat = TRUE)
  heatmap_two_Scores(dds, colData, HRR$gene, 'HRRv2', hrdgenes$gene, 'HRR_70', mode, annot_mode, cancertype, save_heat = TRUE)
  
  if (hrdgeneSigScore){
    heatmap_two_Scores(dds, colData, hrdgeneSig, 'HRDgeneSign', hrdgenes$gene, 'HRR_70', mode, annot_mode, cancertype, save_heat = TRUE)
  }
  
}

heatmap_allScores = function(colData, cancertype, mode, annot_mode, hrdgeneSigScore, save_heat = FALSE){
  ## Creates a heatmap comparing the different scores from the different genesets
  ## Input:
  ## colData (dataframe): Dataframe containg information per sample (rows are samples and the columns are features/information)
  ## cancertype (string): The cancertype has to be without the prefix like 'TCGA' or 'TARGET'
  ## mode (string; default is 'n2'): Indicates if only HRD-High and HRD-Low (mode = n2) or HRD-High,-Medium and -Low (mode ='n3) are compared
  ## annot_mode (string): Indicates if only high and low are annotated ('n2') or high, medium and low ('n3')
  ## hrdgeneSigScore (boolean): If true it means that this score is available and a boxplot can be created
  ## save_heat (boolean): Indicates if the heatmap should be saved (TRUE) or not (FALSE). (default is not saving it)
  
  ### Get the scores
  if (hrdgeneSigScore){
    colData_scores = colData[,c('scoreHRR','scoreCINSARC','scoreMMEJ','scoreSSA','scoreHRRv2','scoreHRDgeneSig')]
  }else{
    colData_scores = colData[,c('scoreHRR','scoreCINSARC','scoreMMEJ','scoreSSA','scoreHRRv2')]
  }
  
  
  colData_scores = t(colData_scores)
  
  ### Get teh types of the samples
  if (mode == 'n2'){
    hrdtypes = colData[match(colnames(colData_scores), rownames(colData)), "gmm_n2_cut"]
  }else{
    hrdtypes = colData[match(colnames(colData_scores), rownames(colData)), "gmm_n3_cut"]
  }
  colnames(colData_scores) = hrdtypes
  
  heat = colData_scores
  
  ### Create the annotation
  if (mode == 'n2'){
    ha = HeatmapAnnotation(HRDtype = colnames(heat), col = list(HRDtype = c('High' = 'darkblue', 'Low' = 'lightblue')))
  }
  else{
    ha = HeatmapAnnotation(HRDtype = colnames(heat), col = list(HRDtype = c('High' = 'darkblue', 'Medium' = 'green','Low' = 'orange')))
  }
  
  
  if (hrdgeneSigScore){
    hm = Heatmap(heat,row_title_rot = 0,show_row_names = FALSE,row_split = factor(rownames(heat), levels = c('scoreHRR','scoreCINSARC','scoreMMEJ','scoreSSA','scoreHRRv2','scoreHRDgeneSig')),show_row_dend = FALSE,show_column_names = FALSE, row_names_gp = gpar(fontsize = 4), bottom_annotation = ha, column_km = 3, name = 'Score')
    
  }else{
    hm = Heatmap(heat,row_title_rot = 0,show_row_names = FALSE,row_split = factor(rownames(heat), levels = c('scoreHRR','scoreCINSARC','scoreMMEJ','scoreSSA','scoreHRRv2')),show_row_dend = FALSE,show_column_names = FALSE, row_names_gp = gpar(fontsize = 4), bottom_annotation = ha, column_km = 3, name = 'Score')
    
  }
  
  if (save_heat){
    png(file=paste('../RNAseq_',cancertype,'/figures/heatmap_all_scores_',mode,'_annotated_',annot_mode, '.png', sep = ''))
  }
  
  
  draw(hm)
  
  if (save_heat){
    dev.off()
    #ggsave(paste('../RNAseq_',cancertype,'/figures/heatmap_scores_',mode,'_annotated_',annot_mode, '.png', sep = ''),hm, width = 8.63, height = 5.71)
  }
  
  
}


newScoreRunSingleCohort = function(projectID){
  cancertype = strsplit(projectID,'-')[[1]]
  if (length(cancertype) > 2){
    cancertype = paste(cancertype[2],cancertype[3],sep = '-')
  } else{
    cancertype = cancertype[2]
  }
  
  ### Get the data
  data = loadData(cancertype)
  print(paste('Loaded data from cohort ', projectID, sep = ''))
  count_matrix = data[[1]]
  colData = data[[2]]
  
  dds_n2 = runDeSeq(count_matrix, colData, mode = 'n2')
  print(paste('Created DeSeq Object for mode n2 for cohort ', projectID, sep = ''))
  
  colData_scores = getScores(dds_n2, colData)
  
  results = createHRDGeneSigScore(dds_n2, colData_scores, mode ='n2', cancertype, min_genes = 5, save_genes = TRUE)
  print(paste('Created scores for cohort ', projectID, sep = ''))
  
  if (results[[1]]){
    enough_geneSig = results[[1]]
    colData_scores = results[[2]]
    geneSig = results[[3]]
  }else{
    enough_geneSig = results[[1]]
    colData_scores = results[[2]]
    geneSig = list()
  }
  
  ### Plot the scores as boxplot, in heatmap compared to HRR_70 genes and all scores together
  plotScores_boxplot_v2(colData_scores, cancertype, enough_geneSig)
  print(paste('Plotted the score for cohort ', projectID, sep = ''))
  
  
  createAll_heatmap_twoScores_v2(dds_n2, colData_scores, mode = 'n2', annot_mode ='n2', cancertype, enough_geneSig, geneSig)
  print(paste('Created heatmap for score comparison for cohort ', projectID, sep = ''))
  
  heatmap_allScores_v2(colData_scores, cancertype, mode = 'n2', annot_mode = 'n2', enough_geneSig, save_heat = TRUE)
  print(paste('Created heatmap with all scores for cohort ', projectID, sep = ''))
  
  
  print(paste('Finished cohort ', projectID, sep = ''))
  
  
}

newScoreRunMultipleCohort = function(cohortList){
  for (cohort in cohortList){
    newScoreRunSingleCohort(cohort)
  }
}



#### functions end ####



#### PCA Plotter ####
#pca_plotter(t(norm_counts_n2_HRR), colData_primary,'HRR genes')
#pca_plotter(t(norm_counts_n2_CINSARC), colData_primary, 'CINSARC genes')
pca_plotter = function(data, colData, genes, save =FALSE, filename='default_pca'){
  pca = prcomp(data, scale. = TRUE)
  pca_df <- data.frame(Case = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2], HRDtype_n2 = colData$gmm_n2_cut, HRDtype_n3 = colData$gmm_n3_cut)
  pca_plot = ggplot(pca_df, aes(x = PC1, y = PC2, color = HRDtype_n3)) +
    geom_point() +
    #labs(title = "PCA normalized counts", x = "PC1", y = "PC2")
    labs(title = paste("PCA normalized counts ",genes,sep=''), x = "PC1", y = "PC2")
  print(pca_plot)
  if (save) {
    ggsave(paste('../data/figures/',filename,'.png',sep = ''),pca_plot,width = 8.63, height = 5.71)
  }
}







#### Test pipeline
projectID = 'TCGA-LUAD'
  
cancertype = strsplit(projectID,'-')[[1]]
if (length(cancertype) > 2){
  cancertype = paste(cancertype[2],cancertype[3],sep = '-')
} else{
  cancertype = cancertype[2]
}
data = loadData(cancertype)
count_matrix = data[[1]]
colData = data[[2]]
  
  
dds_n2 = runDeSeq(count_matrix, colData, mode = 'n2')
dds_n3 = runDeSeq(count_matrix, colData, mode = 'n3')


res = results(dds_n2, contrast = c('gmm_n2_cut', 'High', 'Low'), alpha = 0.05)
res_df = as.data.frame(res)
res_df = res_df[rownames(res_df) %in% hrdgenes$gene,]
res_df = res_df[order(res_df$padj),]
genes_sig = rownames(res_df)[1:20]
heatmap_DEA(dds_n2, res, genes_sig, mode = 'n2', cancertype, save_heat = TRUE)

  
runAnalysis(dds_n2, cancertype, mode = 'n2')
runAnalysis(dds_n3, cancertype, mode = 'n3')
  
colData_scores = getScores(dds_n2, colData)
  
results = createGeneSigScore(dds_n2, colData_scores, mode ='n2', min_genes = 5, max_genes = 200)
  
if (results[[1]]){
  enough_geneSig = results[[1]]
  colData_scores = results[[2]]
  geneSig = results[[3]]
  
}else{
  enough_geneSig = results[[1]]
  colData_scores = results[[2]]
  geneSig = list()
}
  
  
plotScores_boxplot(colData_scores, cancertype, enough_geneSig)
  
createAll_heatmap_twoScores(dds_n2, colData_scores, mode = 'n2', annot_mode ='n2', cancertype, enough_geneSig, geneSig)
  
heatmap_allScores(colData_scores, cancertype, mode = 'n2', annot_mode ='n2', enough_geneSig, save_heat = TRUE)
  
number_samples = nrow(colData_scores)
number_samples_high_n2 = nrow(colData_scores[colData_scores['gmm_n2_cut'] == 'High',])
number_samples_low_n2 = nrow(colData_scores[colData_scores['gmm_n2_cut'] == 'Low',])
  
number_samples_high_n3 = nrow(colData_scores[colData_scores['gmm_n3_cut'] == 'High',])
number_samples_low_n3 = nrow(colData_scores[colData_scores['gmm_n3_cut'] == 'Low',])
number_samples_medium_n3 = nrow(colData_scores[colData_scores['gmm_n3_cut'] == 'Medium',])


test = runSingleCohort('TARGET-AML')

#### Test new score

projectID = 'TCGA-LUAD'

cancertype = strsplit(projectID,'-')[[1]]
if (length(cancertype) > 2){
  cancertype = paste(cancertype[2],cancertype[3],sep = '-')
} else{
  cancertype = cancertype[2]
}

### Get the data
data = loadData(cancertype)
print(paste('Loaded data from cohort ', projectID, sep = ''))
count_matrix = data[[1]]
colData = data[[2]]

dds_n2 = runDeSeq(count_matrix, colData, mode = 'n2')
print(paste('Created DeSeq Object for mode n2 for cohort ', projectID, sep = ''))

colData_scores = getScores(dds_n2, colData)

results = createHRDGeneSigScore(dds_n2, colData_scores, mode ='n2', cancertype, min_genes = 5, save_genes = TRUE)
print(paste('Created scores for cohort ', projectID, sep = ''))

if (results[[1]]){
  enough_geneSig = results[[1]]
  colData_scores = results[[2]]
  geneSig = results[[3]]
}else{
  enough_geneSig = results[[1]]
  colData_scores = results[[2]]
  geneSig = list()
}

### Plot the scores as boxplot, in heatmap compared to HRR_70 genes and all scores together
plotScores_boxplot_v2(colData_scores, cancertype, enough_geneSig)
print(paste('Plotted the score for cohort ', projectID, sep = ''))


createAll_heatmap_twoScores_v2(dds_n2, colData_scores, mode = 'n2', annot_mode ='n2', cancertype, enough_geneSig, geneSig)
print(paste('Created heatmap for score comparison for cohort ', projectID, sep = ''))

heatmap_allScores_v2(colData_scores, cancertype, mode = 'n2', annot_mode = 'n2', enough_geneSig, save_heat = TRUE)
print(paste('Created heatmap with all scores for cohort ', projectID, sep = ''))


print(paste('Finished cohort ', projectID, sep = ''))




newScoreRunSingleCohort('TCGA-LUAD')
#### END Test




#### RUN Pipeline
gmm_results = read.csv('../../GMM_Cutoff/data/gmm_cutoffs_pancancer.csv', sep = ',', row.names = 1)

cohortList = gmm_results$Project.ID

getStatisitcs(cohortList,gmm_results)

cohortList = cohortList[cohortList != 'TARGET-CCSK']

done = c('TARGET-ALL-P2', 'TARGET-AML', 'TARGET-OS', 'TCGA-ACC', 'TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CESC', 'TCGA-CHOL', 'TCGA-COAD','TCGA-DLBC','TCGA-ESCA','TCGA-GBM', 'TCGA-HNSC','TCGA-KICH','TCGA-KIRC')
cohortList = cohortList[cohortList != done]
#KIRP is start thrusday 
runMultipleCohorts(cohortList)



newScoreRunMultipleCohort(cohortList)

