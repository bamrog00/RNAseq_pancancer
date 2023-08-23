library(ggplot2)
library(stats)
library(ggpubr)
library(corrplot)
library(ComplexHeatmap)

no = c('TARGET-ALL-P2', 'TARGET-AML','TCGA-CESC', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-DLBC', 'TCGA-ESCA', 'TCGA-GBM', 'TCGA-LAML', 'TCGA-MESO', 'TCGA-PCPG', 'TCGA-READ', 'TCGA-SKCM', 'TCGA-TGCT', 'TCGA-THCA', 'TCGA-THYM', 'TCGA-UCS', 'TCGA-UVM')


getSubfolders = function(folder_path) {
  ## Gets the folder names
  subfolders = list.dirs(folder_path, recursive = FALSE)
  subfolder_names <- basename(subfolders)
  return(subfolder_names)
}

stackData = function(methode, folders){
  ## Get the data and stacks all cancer types togehter
  data = data.frame()
  
  for (folder in folders){
    
    if (file.exists(paste('../',folder,'/results_',methode,'/colData_scores.csv', sep = ''))) {
      print("File exists!")
    } else {
      print(paste('No normCounts for cohort ', folder, sep = ''))
      next
    }
    
    
    colData = read.csv(paste('../',folder,'/results_',methode,'/colData_scores.csv', sep = ''), sep = ',', row.names = 1)
    
    data = rbind(data, colData)
  }
  return (data)
}



folders = getSubfolders(paste('../', version, sep = ''))

folders = folders[grepl("^RNAseq_", folders)]
folders = folders[folders != 'RNAseq_CCSK']



### GMM

methode = 'GMM'

data_gmm = stackData(methode, folders)

data_gmm_filtered = data_gmm[!(data_gmm$Project.ID %in% no), ]



# Create the boxplot
cinsarc_bp_gmm = ggplot(data_gmm_filtered, aes(x = Project.ID, y = scoreCINSARC, fill = n2_cut)) +
  geom_boxplot() +
  labs(x = "Project ID", y = "ScoreCINSARC", fill = "HRD-Type")  +
  ggtitle("ScoreCINSARC, GMM") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  stat_compare_means(aes(group = n2_cut), method = 't.test', label ="p.signif") +
  scale_fill_manual(values = c("High" = "#1F77B4", "Low" = "#E377C2"), 
                    labels = c("HRD-high", "HRD-low"))
print(cinsarc_bp_gmm)
ggsave('../data/figures/cinsarc_gmm_bp.png',cinsarc_bp_gmm,width = 8.63, height = 5.71)
#stat_compare_means(aes(group = n2_cut), method = 't.test',label.y = c(200, 215, 230,245,260,275,290,305,320,290,300,310,320,330,340,350,360), size = 2.3) +
  
HRR70_bp_gmm = ggplot(data_gmm_filtered, aes(x = Project.ID, y = scoreHRR, fill = n2_cut)) +
  geom_boxplot() +
  labs(x = "Project ID", y = "ScoreHRR70", fill = "HRD-Type")  +
  ggtitle("scoreHRR70, GMM") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  stat_compare_means(aes(group = n2_cut), method = 't.test', label ="p.signif") +
  scale_fill_manual(values = c("High" = "#1F77B4", "Low" = "#E377C2"), 
                    labels = c("HRD-high", "HRD-low"))
print(HRR70_bp_gmm)
ggsave('../data/figures/hrr70_gmm_bp.png',HRR70_bp_gmm,width = 8.63, height = 5.71)


MMEJ_bp_gmm = ggplot(data_gmm_filtered, aes(x = Project.ID, y = scoreMMEJ, fill = n2_cut)) +
  geom_boxplot() +
  labs(x = "Project ID", y = "ScoreMMEJ", fill = "HRD-Type")  +
  ggtitle("scoreMMEJ, GMM") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  stat_compare_means(aes(group = n2_cut), method = 't.test', label ="p.signif") +
  scale_fill_manual(values = c("High" = "#1F77B4", "Low" = "#E377C2"), 
                    labels = c("HRD-high", "HRD-low"))
print(MMEJ_bp_gmm)
ggsave('../data/figures/mmek_gmm_bp.png',MMEJ_bp_gmm,width = 8.63, height = 5.71)


SSA_bp_gmm = ggplot(data_gmm_filtered, aes(x = Project.ID, y = scoreSSA, fill = n2_cut)) +
  geom_boxplot() +
  labs(x = "Project ID", y = "ScoreSSA", fill = "HRD-Type")  +
  ggtitle("scoreSSA, GMM") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  stat_compare_means(aes(group = n2_cut), method = 't.test', label ="p.signif") +
  scale_fill_manual(values = c("High" = "#1F77B4", "Low" = "#E377C2"), 
                    labels = c("HRD-high", "HRD-low"))
print(SSA_bp_gmm)
ggsave('../data/figures/ssa_gmm_bp.png',SSA_bp_gmm,width = 8.63, height = 5.71)

HRR_bp_gmm = ggplot(data_gmm_filtered, aes(x = Project.ID, y = scoreHRRv2, fill = n2_cut)) +
  geom_boxplot() +
  labs(x = "Project ID", y = "ScoreHRR", fill = "HRD-Type")  +
  ggtitle("scoreHRR, GMM") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  stat_compare_means(aes(group = n2_cut), method = 't.test', label ="p.signif") +
  scale_fill_manual(values = c("High" = "#1F77B4", "Low" = "#E377C2"), 
                    labels = c("HRD-high", "HRD-low"))
print(HRR_bp_gmm)
ggsave('../data/figures/hrr_gmm_bp.png',HRR_bp_gmm,width = 8.63, height = 5.71)



###### Kmeans

data_kmean = stackData('Kmeans', folders)

data_kmean_filtered = data_kmean[!(data_kmean$Project.ID %in% no), ]


cinsarc_bp_km = ggplot(data_kmean_filtered, aes(x = Project.ID, y = scoreCINSARC, fill = n2_cut)) +
  geom_boxplot() +
  labs(x = "Project ID", y = "ScoreCINSARC", fill = "HRD-Type")  +
  ggtitle("ScoreCINSARC, K-means") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  stat_compare_means(aes(group = n2_cut), method = 'wilcox.test', label ="p.signif") +
  scale_fill_manual(values = c("High" = "#1F77B4", "Low" = "#E377C2"), 
                    labels = c("HRD-high", "HRD-low"))
print(cinsarc_bp_km)
ggsave('../data/figures/cinsarc_km_bp_wilcox.png',cinsarc_bp_km,width = 8.63, height = 5.71)
#stat_compare_means(aes(group = n2_cut), method = 't.test',label.y = c(200, 215, 230,245,260,275,290,305,320,290,300,310,320,330,340,350,360), size = 2.3) +

HRR70_bp_km = ggplot(data_kmean_filtered, aes(x = Project.ID, y = scoreHRR, fill = n2_cut)) +
  geom_boxplot() +
  labs(x = "Project ID", y = "ScoreHRR70", fill = "HRD-Type")  +
  ggtitle("scoreHRR70, K-means") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  stat_compare_means(aes(group = n2_cut), method = 'wilcox.test', label ="p.signif") +
  scale_fill_manual(values = c("High" = "#1F77B4", "Low" = "#E377C2"), 
                    labels = c("HRD-high", "HRD-low"))
print(HRR70_bp_km)
ggsave('../data/figures/hrr70_km_bp_wilcox.png',HRR70_bp_km,width = 8.63, height = 5.71)


MMEJ_bp_km = ggplot(data_kmean_filtered, aes(x = Project.ID, y = scoreMMEJ, fill = n2_cut)) +
  geom_boxplot() +
  labs(x = "Project ID", y = "ScoreMMEJ", fill = "HRD-Type")  +
  ggtitle("scoreMMEJ, K-means") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  stat_compare_means(aes(group = n2_cut), method = 'wilcox.test', label ="p.signif") +
  scale_fill_manual(values = c("High" = "#1F77B4", "Low" = "#E377C2"), 
                    labels = c("HRD-high", "HRD-low"))
print(MMEJ_bp_km)
ggsave('../data/figures/mmek_km_bp_wilcox.png',MMEJ_bp_km,width = 8.63, height = 5.71)


SSA_bp_km = ggplot(data_kmean_filtered, aes(x = Project.ID, y = scoreSSA, fill = n2_cut)) +
  geom_boxplot() +
  labs(x = "Project ID", y = "ScoreSSA", fill = "HRD-Type")  +
  ggtitle("scoreSSA, K-means") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  stat_compare_means(aes(group = n2_cut), method = 'wilcox.test', label ="p.signif") +
  scale_fill_manual(values = c("High" = "#1F77B4", "Low" = "#E377C2"), 
                    labels = c("HRD-high", "HRD-low"))
print(SSA_bp_km)
ggsave('../data/figures/ssa_km_bp_wilcox.png',SSA_bp_km,width = 8.63, height = 5.71)

HRR_bp_km = ggplot(data_kmean_filtered, aes(x = Project.ID, y = scoreHRRv2, fill = n2_cut)) +
  geom_boxplot() +
  labs(x = "Project ID", y = "ScoreHRR", fill = "HRD-Type")  +
  ggtitle("scoreHRR, K-means") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  stat_compare_means(aes(group = n2_cut), method = 'wilcox.test', label ="p.signif") +
  scale_fill_manual(values = c("High" = "#1F77B4", "Low" = "#E377C2"), 
                    labels = c("HRD-high", "HRD-low"))
print(HRR_bp_km)
ggsave('../data/figures/hrr_km_bp_wilcox.png',HRR_bp_km,width = 8.63, height = 5.71)



scores = c('scoreCINSARC', 'scoreHRR', 'scoreHRRv2', 'scoreMMEJ', 'scoreSSA')

ids = unique(data_kmean_filtered$Project.ID)


empty_df <- data.frame(matrix(ncol = length(ids), nrow = length(scores)))
colnames(empty_df) = ids
rownames(empty_df) = scores

for (score in scores){
  for (cancer in ids){
    subdata = data_kmean_filtered[data_kmean_filtered['Project.ID'] == cancer,]
    correlation = cor(subdata$HRD_sum,subdata[,score], method = 'spearman')
    empty_df[score,cancer] = correlation
  }
}

corrmatrix = as.matrix(empty_df)


# Create heatmap
#heatmap(corrmatrix, Rowv = NA, Colv = NA, scale = "none",
        #labRow = row.names(empty_df), labCol = colnames(empty_df))

png(file='../data/figures/correlation_scores_hrd_hm_spearman.png')
Heatmap(corrmatrix, show_column_dend = FALSE, show_row_dend = FALSE, name = 'Spearman correlation', column_title = 'Correlation gene-scores and HRDsum')
dev.off()


empty_df <- data.frame(matrix(nrow = length(ids), ncol = length(scores)))
rownames(empty_df) = ids
colnames(empty_df) = c(scores)

for (score in scores){
  for (cancer in ids){
    subdata = data_kmean_filtered[data_kmean_filtered['Project.ID'] == cancer,]
    formula <- as.formula(paste(score, "~ n2_cut"))
    res = wilcox.test(formula, data = subdata,exact = FALSE, paired = FALSE)
    empty_df[cancer,score] = res$p.value
  }
}

write.csv(empty_df,'../data/p_values_wilcoxon_genesetscores_HRDtype.csv', row.names = TRUE)

