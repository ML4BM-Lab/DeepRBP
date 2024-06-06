
rm(list = ls())
setwd("C:/Users/lvalcarcel/OneDrive - Tecnun/Postdoc/102 - Joseba DeepSF/2024-05-30_Figura2/data_paper_figure2")
# setwd("C:/Users/Luisvi/OneDrive - ---/Postdoc/102 - Joseba DeepSF/2024-05-30_Figura2/data_paper_figure2")

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(latex2exp)

plotlist <- list()
plotlist2 <- list()

df_combined1 <- data.table::fread('ordered_GxS_scores_of_k_genes_with_more_rbps_in_all_Postar.csv') %>% as.data.frame()
df_combined2 <- data.table::fread('ordered_GxS_scores_of_k_rbps_with_more_genes_in_all_Postar.csv') %>% as.data.frame()
df_scores <- data.table::fread('df_knockout_ref_DeepLIFT_z_scores_GxRBPs.csv') %>% as.data.frame()

df_combined1$Postar <- as.character(df_combined1$Postar)
df_combined1$Postar[is.na(df_combined1$Postar)] <- "NA"
df_combined1 <- df_combined1 %>% mutate(Postar = factor(Postar, levels = c("1","0", "NA"), labels = c("1","0","NA")))

df_combined2$Postar <- as.character(df_combined2$Postar)
df_combined2$Postar[is.na(df_combined2$Postar)] <- "NA"
df_combined2 <- df_combined2 %>% mutate(Postar = factor(Postar, levels = c("1","0","NA"), labels = c("1","0","NA")))


head(df_combined1)
head(df_combined2)

actual_rbps <- c( 'UTP3', 'SF3B1', 'AQR', 'UCHL5')
actual_rbps2 <- c( '                       UTP3', '          SF3B1', '           AQR', '           UCHL5')
actual_genes <- c('ENSG00000082898', 'ENSG00000126883', 'ENSG00000109685', 'ENSG00000079805')

actual_df_combined1 <- df_combined1 %>% filter(RBP_name %in% actual_rbps) %>% 
  mutate(RBP_name = factor(RBP_name, levels = actual_rbps, labels = actual_rbps2))

actual_df_combined2 <- df_combined2 %>% filter(Gene_ID %in% actual_genes) %>% 
  mutate(Gene_ID = factor(Gene_ID, levels = actual_genes))



# Add p-values onto the bar plots
stat.test <- actual_df_combined1 %>%
  group_by(RBP_name) %>%
  rstatix::wilcox_test(Scores ~ Postar, comparisons = list(c("1","0"))) %>%
  # %>%#, dodge = 0.8)# %>%
  # filter(group1 == "1") %>% filter(group2 == "0")  %>%
  add_xy_position(fun = "max", x = "RBP_name") %>% 
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")



plotlist[[1]] <- ggboxplot(actual_df_combined1, notch = FALSE,
                           x = "RBP_name", y = "Scores", fill = "Postar") + xlab("RBP name") + 
  scale_fill_manual(values = c("#7fc97f", "#beaed4", "#fdc086"),  breaks = c("1", "0", "NA")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_continuous(expand = c(0,0.25)) + 
  # theme(legend.position = "right", legend.spacing = unit(3, "cm")) + 
  # labs(fill="Postar     ") + 
  #           "Transcripts"
  stat_pvalue_manual( stat.test,  label = "p.adj.signif")#, tip.length = 0.01, bracket.nudge.y = -2)
# plotlist[[1]]

# Add p-values onto the bar plots
stat.test <- actual_df_combined2 %>%
  group_by(Gene_ID) %>%
  rstatix::wilcox_test(Scores ~ Postar, comparisons = list(c("1","0"))) %>%
  # %>%#, dodge = 0.8)# %>%
  # filter(group1 == "1") %>% filter(group2 == "0")  %>%
  add_xy_position(fun = "max", x = "Gene_ID") %>% 
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  filter(Gene_ID %in% actual_genes) 

stat.test <- merge(stat.test, 
                   actual_df_combined2 %>% filter(Postar %in% c("1", "0")) %>% 
                     group_by(Gene_ID) %>% summarise(y.position.2 = max(Scores)+2))

plotlist[[2]] <- ggboxplot(actual_df_combined2,
                           x = "Gene_ID", y = "Scores", fill = "Postar") + xlab("Gene ID") + 
  scale_fill_manual(values = c("#7fc97f", "#beaed4", "#fdc086"),  breaks = c("1", "0", "NA")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  # labs(fill="Postar     ") +
  scale_y_continuous(expand = c(0,0.25)) +
  stat_pvalue_manual( stat.test,  label = "p.adj.signif", y.position = "y.position.2")#, tip.length = 0.01, bracket.nudge.y = -2)
# plotlist[[2]]

ggarrange(plotlist = plotlist, nrow = 1, ncol = 2, labels = c("A", "B"), common.legend = T, legend = "right")


rbp_interest_list <- c('MBNL1', 'RBM47', 'FUS', 'TAF15', 'TARDBP', 'TARDBP')
experiment_list <- c('PRJEB39343', 'GSE75491', 'GSE77702-kd1', 'GSE77702-kd2', 'GSE77702-kd3', 'GSE136366')
limma_list <- c('limma_HFE_MBNL1_results.csv', 
                'limma_SRR296_RBM47_results.csv',
                'limma_GSE77702-kd1_FUS_results.csv',
                'limma_GSE77702-kd2_TAF15_results.csv',
                'limma_GSE77702-kd3_TARDBP_results.csv',
                'limma_GSE136366_TARDBP_results.csv')


for (i in 1:length(rbp_interest_list)){
  
  # i = 1
  
  trans_file='df_DeepLIFT_zscores_TxRBPs' # name of the transcript score matrix. ESTA LA NECESITAMOS SIEMPRE.
  genes_file='df_z_score_GxRBP_knockout'
  
  
  
  
  rbp_interest = rbp_interest_list[i]
  experiment = experiment_list[i]
  file_path = limma_list[i]
  df_limma = data.table::fread(limma_list[i]) %>% as.data.frame()
  
  # df_limma.set_index('Transcript_ID', inplace=True)
  df_score_trans = data.table::fread(paste0(trans_file, '_', experiment, '.csv')) %>% 
    as.data.frame() %>% tibble::column_to_rownames("V1") %>% abs()
  df_significant = df_limma %>% filter(Transcript_ID %in% rownames(df_score_trans))
  df_significant = df_significant %>% filter(adj.P.Val < 0.05)
  # df_significant['Transcript_ID'] = df_significant.index
  
  # df_significant.reset_index(inplace=True, drop=True)
  significant_samples = df_significant$Transcript_ID
  
  df_score_gns = data.table::fread(paste0(genes_file, '_', experiment, '.csv')) %>% 
    as.data.frame() %>% tibble::column_to_rownames("V1")
  
  significant_samples_gns = df_significant['Gene_ID'] %>% unique() %>% unlist() %>% as.character()
  
  # plot_boxplot_with_annotations(
  #   ax=ax,
  #   df=df_score_gns,
  #   significant_samples=significant_samples_gns,
  #   rbp_interest=rbp_interest,
  #   experiment=experiment,
  #   data_type=data_type
  # )
  
  # Identify the differential expressed transcripts
  
  df_scores = df_score_gns[, rbp_interest,drop = F] %>% abs() %>% as.data.frame()
  df_scores['DE_Limma'] = 'No'
  df_scores[significant_samples_gns, 'DE_Limma'] = 'Yes'
  
  # Hacer una copia del DataFrame original
  df_original = df_scores
  df_melted = reshape2::melt(df_scores) %>% 
    mutate(Transcripts = factor(DE_Limma, levels = c("Yes", "No"), labels = c("DE", "non-DE"))) %>%
    mutate(Scores = log10(value+1))
  
  
  # Statistical test
  stat.test <- df_melted %>% rstatix::wilcox_test(Scores ~ Transcripts) %>%   
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj") %>% add_xy_position(x = "Transcripts")
  
  
  plotlist2[[i]] <- ggboxplot(df_melted, 
                              add = "jitter", add.params = list(shape = 21, alpha = 0.5, color = NULL), 
                              x = "Transcripts", y = "Scores", fill = "Transcripts") + 
    ylab(TeX(r'(Scores in $log_{10}-scale$)')) + xlab(rbp_interest) + 
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
    stat_pvalue_manual(stat.test, label = "p.adj.signif") + scale_y_continuous(expand = c(0,0.25)) + 
    ggtitle(experiment) + scale_fill_manual(values = c("#d95f02", "#7570b3"),  breaks = c("DE", "non-DE"))
  
}


p <- ggarrange(ggarrange(ggplot()+theme_void()),
               ggarrange(ggarrange(plotlist = plotlist, nrow = 1, ncol = 2, labels = c("A", "B"), common.legend = T, legend = "right"),
                    ggplot()+theme_void(), nrow = 1, ncol = 2, widths = c(100,5)),
          ggarrange(ggplot()+theme_void()),
          ggarrange(plotlist = plotlist2, nrow = 3, ncol = 2, labels = c("C", "D", "E", "F", "G", "H"), common.legend = T, legend = "right"),
          nrow = 4, ncol = 1, heights = c(0.5,30,1,55))

ggsave("zscores_plot_lvv.pdf", plot = p, height = 140, width = 100, units = "mm",scale = 2)

