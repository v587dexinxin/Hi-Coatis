setwd('H:/work/niulongjian/HiRPC_processed_data/plots/HiRPC_classify/K562_VS_Hemin')

library(ggforce)
library(reshape2)
library(tidyr)
library('ggpubr')

data <- reshape2::melt(Titanic)
data <- gather_set_data(data, 1:4)

ggplot(data, aes(x, id = id, split = y, value = value)) +
  geom_parallel_sets(aes(fill = Class), alpha = 0.3, axis.width = 0.1, show.legend = F) +
  geom_parallel_sets_axes(axis.width = 0.3, colour = "black", fill = "white" ) +  #colour和fill分别修改边框色和填充艿
  geom_parallel_sets_labels(colour = 'black', angle = 0)+ # angle修改方框填充文字的角庿
  theme_classic()


df_K562_Merge_GeneGroup <- read.table('H:/work/niulongjian/HiRPC_processed_data/plots/HiRPC_classify/classify_K562/K562_union_gene_classify.txt', header = T)
df_Hemin_Merge_GeneGroup <- read.table('H:/work/niulongjian/HiRPC_processed_data/plots/HiRPC_classify/classify_Hemin/Hemin_union_gene_classify.txt', header = T)


K562_c1_Genes <- df_K562_Merge_GeneGroup[df_K562_Merge_GeneGroup$Level == 'c1',]$Gene_Name
K562_c2_Genes <- df_K562_Merge_GeneGroup[df_K562_Merge_GeneGroup$Level == 'c2',]$Gene_Name
K562_c3_Genes <- df_K562_Merge_GeneGroup[df_K562_Merge_GeneGroup$Level == 'c3',]$Gene_Name
K562_c4_Genes <- df_K562_Merge_GeneGroup[df_K562_Merge_GeneGroup$Level == 'c4',]$Gene_Name

Hemin_c1_Genes <- df_Hemin_Merge_GeneGroup[df_Hemin_Merge_GeneGroup$Level == 'c1',]$Gene_Name
Hemin_c2_Genes <- df_Hemin_Merge_GeneGroup[df_Hemin_Merge_GeneGroup$Level == 'c2',]$Gene_Name
Hemin_c3_Genes <- df_Hemin_Merge_GeneGroup[df_Hemin_Merge_GeneGroup$Level == 'c3',]$Gene_Name
Hemin_c4_Genes <- df_Hemin_Merge_GeneGroup[df_Hemin_Merge_GeneGroup$Level == 'c4',]$Gene_Name



df_K562_Merge_GeneGroup[, c('Gene_Name', 'Level')]
df_Hemin_Merge_GeneGroup[, c('Gene_Name', 'Level')]



## 2. PCI_EX genes sankey plot.
groups <- c('c1', 'c2', 'c3', 'c4')

df_group <- data.frame(G_group=character(),
                       K_group=character(),
                       GeneNum=numeric())

for (G_group in groups){
  for (K_group in groups){
    cat("K562", G_group, "VS" , "Hemin", K_group, '\n')
    # g_group_genes
    gene_num <- length(intersect(df_K562_Merge_GeneGroup[df_K562_Merge_GeneGroup$Level == G_group,]$Gene_Name,
                                 df_Hemin_Merge_GeneGroup[df_Hemin_Merge_GeneGroup$Level == K_group,]$Gene_Name))
    df_tmp <- data.frame(G_group=G_group,
                         K_group=K_group,
                         GeneNum=gene_num)
    df_group <- rbind(df_group, df_tmp)
  }
}

df_group <- gather_set_data(df_group, x = 1:2)
df_group$y <- factor(df_group$y, levels = c('c1', 'c2', 
                                            'c3', 'c4'))
df_group$G_group <- factor(df_group$G_group, levels = c('c1', 'c2', 
                                                        'c3', 'c4'))
p_group_sankey <- ggplot(df_group, aes(x, id = id, split = y, value = GeneNum)) +
  geom_parallel_sets(aes(fill = G_group), alpha = 0.5, axis.width = 0.1, show.legend = F) +
  geom_parallel_sets_axes(axis.width = 0.2, colour = "black", fill = "white" ) +  #colour和fill分别修改边框色和填充艿
  geom_parallel_sets_labels(colour = 'black', angle = 0)+ # angle修改方框填充文字的角庿
  scale_fill_manual(values = c("#EE0000", "#CD5C5C", "#EE7600", "#EEAD0E", "#00008B", "#4876FF", "#454545", "#c2c2c2"))+
  theme_void()


ggarrange(p_group_sankey, nrow=1, ncol=1)
ggsave('sankey_plot.pdf', width = 9, height = 16)

