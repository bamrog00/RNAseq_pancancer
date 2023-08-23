library(ggplot2)
library(stats)
library(ggpubr)
library(corrplot)
library(forcats)
library(dplyr)

cohorts = c('RNAseq_LGG', 'RNAseq_LIHC', 'RNAseq_LUAD', 'RNAseq_PAAD', 'RNAseq_PRAD', 'RNAseq_SARC', 'RNAseq_ACC', 'RNAseq_KICH','RNAseq_KIRC', 'RNAseq_KIRP', 'RNAseq_UCEC', 'RNAseq_BRCA')
cancertypes = c('TCGA-LGG', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-PAAD', 'TCGA-PRAD', 'TCGA-SARC', 'TCGA-ACC', 'TCGA-KICH','TCGA-KIRC', 'TCGA-KIRP', 'TCGA-UCEC', 'TCGA-BRCA')
geneCounts_up <- data.frame(Gene = character(), Count = integer(), stringsAsFactors = FALSE)
geneCounts_down <- data.frame(Gene = character(), Count = integer(), stringsAsFactors = FALSE)

test_counts_up = data.frame(Gene = character(), Count = integer(), Cancertype = character(), stringsAsFactors = FALSE)
test_counts_down = data.frame(Gene = character(), Count = integer(), Cancertype = character(), stringsAsFactors = FALSE)

for (cancertype in cohorts){
  hrr_sign_genes = read.csv(paste('../',cancertype,'/results_Kmeans/hrd_significant_genes_n2.csv', sep = ''), sep = ',', row.names = 1)
  
  genes = rownames(hrr_sign_genes)
  
  split_string <- strsplit(cancertype, "_")
  
  # Extract the part after the underscore
  name <- split_string[[1]][2]
  
  for (gene in genes) {
    # Check if the gene is already present in the geneCounts dataframe
    lfc = hrr_sign_genes[gene,'log2FoldChange']
    if (lfc < 0){
      test_counts_down = rbind(test_counts_down, data.frame(Gene = gene, Count = 1, Cancertype = name, stringsAsFactors = FALSE))
      if (gene %in% geneCounts_down$Gene) {
        # If the gene is already in the dataframe, increment its count by 1
        geneCounts_down$Count[geneCounts_down$Gene == gene] = geneCounts_down$Count[geneCounts_down$Gene == gene] + 1
      } else {
        # If the gene is new, add it to the geneCounts dataframe with a count of 1
        geneCounts_down = rbind(geneCounts_down, data.frame(Gene = gene, Count = 1, stringsAsFactors = FALSE))
      }
    }else{
      test_counts_up = rbind(test_counts_up, data.frame(Gene = gene, Count = 1, Cancertype = name, stringsAsFactors = FALSE))
      if (gene %in% geneCounts_up$Gene) {
        # If the gene is already in the dataframe, increment its count by 1
        geneCounts_up$Count[geneCounts_up$Gene == gene] = geneCounts_up$Count[geneCounts_up$Gene == gene] + 1
      } else {
        # If the gene is new, add it to the geneCounts dataframe with a count of 1
        geneCounts_up = rbind(geneCounts_up, data.frame(Gene = gene, Count = 1, stringsAsFactors = FALSE))
      }
    }
    
    

  }
}



geneCounts_up = geneCounts_up[order(-geneCounts_up$Count), ]
geneCounts_up$Gene <- factor(geneCounts_up$Gene, levels = geneCounts_up$Gene)

geneCounts_down = geneCounts_down[order(-geneCounts_down$Count), ]
geneCounts_down$Gene <- factor(geneCounts_down$Gene, levels = geneCounts_down$Gene)

count_plot_up  = ggplot(geneCounts_up, aes(x = Gene, y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Genes") +
  ylab("Frequency") +
  ggtitle("HRR significant upregulated genes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

print(count_plot_up)
ggsave('../data/figures/hrr_signi_genes_counts_up.png',count_plot_up,width = 8.63, height = 5.71)


count_plot_down  = ggplot(geneCounts_down, aes(x = Gene, y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Genes") +
  ylab("Frequency") +
  ggtitle("HRR significant downregulated genes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

print(count_plot_down)
ggsave('../data/figures/hrr_signi_genes_counts_down.png',count_plot_down,width = 8.63, height = 5.71)



category_count = table(test_counts_up$Gene)

test_counts_up$Gene_Count <- category_count[test_counts_up$Gene]

test_counts_up <- test_counts_up[order(test_counts_up$Gene_Count, decreasing = TRUE), ]


stacked_plot_up  = ggplot(test_counts_up, aes(x  = reorder(Gene, -Gene_Count), y = Count, fill = Cancertype)) +
  geom_bar(position = 'stack', stat = "identity") +
  xlab("Genes") +
  ylab("Number of appearances") +
  theme_bw()+
  ggtitle("HRR significant upregulated genes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

print(stacked_plot_up)
ggsave('../data/figures/hrr_signi_genes_counts_upregulated_per_cancertype_new.png',stacked_plot_up,width = 8.63, height = 5.71)

category_count = table(test_counts_down$Gene)

test_counts_down$Gene_Count <- category_count[test_counts_down$Gene]

test_counts_down <- test_counts_down[order(test_counts_down$Gene_Count, decreasing = TRUE), ]

stacked_plot_down  = ggplot(test_counts_down, aes(x = reorder(Gene, -Gene_Count), y = Count, fill = Cancertype)) +
  geom_bar(position = 'stack', stat = "identity") +
  xlab("Genes") +
  ylab("Number of appearances") +
  theme_bw()+
  ggtitle("HRR significant downregulated genes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

print(stacked_plot_down)
ggsave('../data/figures/hrr_signi_genes_counts_downregulated_per_cancertype_new.png',stacked_plot_down,width = 8.63, height = 5.71)



for (i in 1:length(cohorts)){
  cohort = cohorts[[i]]
  cancertype = cancertypes[[i]]
  if (i == 1){
    hrr_sign_genes = read.csv(paste('../',cohort,'/results_Kmeans/hrd_significant_genes_n2.csv', sep = ''), sep = ',', row.names = 1)
    hrr_sign_genes$CancerType = cancertype
    hrr_sign_genes$Gene = rownames(hrr_sign_genes)
  }else{
    hrr_sign_genes_other = read.csv(paste('../',cohort,'/results_Kmeans/hrd_significant_genes_n2.csv', sep = ''), sep = ',', row.names = 1)
    hrr_sign_genes_other$CancerType = cancertype
    hrr_sign_genes_other$Gene = rownames(hrr_sign_genes_other)
    hrr_sign_genes = rbind(hrr_sign_genes,hrr_sign_genes_other)
  }
}
hrr_sign_genes = hrr_sign_genes[,c('CancerType', 'Gene', 'log2FoldChange','pvalue', 'padj')]
write.csv(hrr_sign_genes,'../data/hrd_signature.csv', row.names = FALSE)

for (cancertype in unique(hrr_sign_genes$CancerType)){
  subset = hrr_sign_genes[hrr_sign_genes$CancerType == cancertype,]
  print(cancertype)
  print(nrow(subset))
}
