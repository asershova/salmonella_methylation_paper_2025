#install.packages('skimr', 'here', 'kableExtra')
#install.packages("patchwork")
#install.packages("pROC")
library("hrbrthemes")
library(patchwork)
library("pROC")
library(ggplot2)
library(dplyr)
library(tidyverse)
library(readxl)
library(here)
library(skimr) # install.packages('skimr')
library(kableExtra) # install.packages('kableExtra')
library(pROC) #version 1.18.5
library("ggpubr")
site = "GAATTC"
paper_pal <-paletteer::paletteer_d("ggsci::category20_d3")
status_pal <- paper_pal[1:4]
violin_color <- c("#C49C94FF","#F7B6D2FF","#C7C7C7FF")
chi_test_color <- c("#7F7F7FFF","#C7C7C7FF")
status_pal_2 <- c("#AEC7E8FF","#FFBB78FF","#98DF8AFF","#FF9896FF")
sample_pal <- c("#BCBD22FF","#DBDB8DFF","#17BECFFF","#9EDAE5FF")

col_names = c('chrom', 'm_start', 'm_stop', 'met_type', 'total_reads', 'met_dir', 'met_st_1', 'met_st_2', 'color','reads_all', 'methyl', 'methyl_reads', 'non_methyl_reads',
              'v1','v2','v3', 'v4', 'v5', 'chrom2', 'site_start', 'site_stop', 'site', 'v6', 'dir',
              'chrom3', 'f_start', 'f_stop', 'f_dir','f_type','locus_tag','id1','id2','locus_tag1','id11','id21','d_chrom','d_start', 'd_stop','d_old','d_id','domain')
#working_dir = "/home/anna/old/mnt/Documents/Methylome/results/experiments/nanopore_salmonella/sal_ent_barcodes_6mA/aligned/0_4_4_bed/mod_data" # the path to the working directory, you should put your own
working_dir = "/app/data" # the path to the working directory, you should put your own
results_dir = "/app/results"
data_08 <- read.csv(paste(working_dir,"WT.GAATTC.bed", sep="/"), sep="", header = FALSE,
                    col.names = col_names)
data_14 <- read.csv(paste(working_dir,"M.EcoRI.GAATTC.bed", sep="/"), sep="", header = FALSE,
                    col.names = col_names)


paste("Loaded.")
data_08$sample<- 'LSP-1'
data_14$sample<- 'LSP-M.EcoRI'

data_08$group<- 'WT'
data_14$group<- 'M.EcoRI'

data_08$gaattc_stat <- 0
data_14$gaattc_stat <- 1

data_08 <- data_08 %>% 
  mutate(rel_pos = ifelse(met_dir=="+", 
                          data_08$m_start-data_08$site_start, 
                          data_08$site_stop-1-data_08$m_start))  

data_14 <- data_14 %>% 
  mutate(rel_pos = ifelse(met_dir=="+", 
                          data_14$m_start-data_14$site_start, 
                          data_14$site_stop-1-data_14$m_start)) 

data_all_ecori <-rbind(data_08,data_14)
#data_14_gaattc Salmonella low coverage not included
#M.EcoRI GAATTC
data_all_methyl_a_ecori_plus <- data_all_ecori[data_all_ecori$met_dir=="+"&data_all_ecori$dir=="+",]
data_all_methyl_a_ecori_minus <- data_all_ecori[data_all_ecori$met_dir=="-"&data_all_ecori$dir=="-",]
data_all_methyl_a_ecori <- rbind(data_all_methyl_a_ecori_plus, data_all_methyl_a_ecori_minus)
data_all_methyl_a_ecori_2 <- data_all_methyl_a_ecori[data_all_methyl_a_ecori$rel_pos==2,]
#data_all_methyl_a_chrom_ecori <- data_all_ecori[data_all_ecori$rel_pos==7 & data_all_ecori$chrom=='NZ_CP008706.1',]

data_all_methyl_a_ecori_2_chrom <- data_all_methyl_a_ecori_2[data_all_methyl_a_ecori_2$chrom=="ST4-74.fa",]
data_all_methyl_a_ecori_2_p <- data_all_methyl_a_ecori_2[data_all_methyl_a_ecori_2$chrom!="ST4-74.fa",]
#violin plot
reordered_data<-data_all_methyl_a_ecori_2 %>%
  mutate(group = fct_relevel(group, 
                             "M.EcoRI","WT"))
summary<- reordered_data %>%
  group_by(chrom, group) %>%
  summarize(n = n(), .groups="keep")
summary
reordered_data <- reordered_data %>% 
  mutate(dna_type = ifelse(chrom=='ST4-74.fa', 
                          "chromosome", 
                          "plasmid"))  

ecori_col <- c("#1F77B4FF", "#AEC7E8FF")
#ROC curve and thresholds

roc_gaattc <- roc(data_all_methyl_a_ecori_2$gaattc_stat, data_all_methyl_a_ecori_2$methyl)

coords(roc_gaattc, "best", ret=c("threshold", "specificity", "1-npv"))
coords(roc=roc_gaattc, x=35, input="threshold")
#threshold specificity       1-npv
#threshold    34.535   0.9924433 0.003162555


p1_roc<-ggroc(roc_gaattc,colour = "#1F77B4FF", linetype = 1, size = 1) +
  theme_ipsum(base_family = 'Arial', base_size = 12) +
  theme(
    axis.title.x = element_text(size = 12, vjust = 0.5, hjust = 0.5),
    axis.title.y = element_text(size = 12, vjust = 0.5, hjust = 0.5))


#Jitter plot
p1j <- ggplot(data_all_methyl_a_ecori[data_all_methyl_a_ecori$total_reads>1,], 
       aes(x = factor(rel_pos), y = methyl)) + 
  geom_jitter(aes(color = group), position=position_jitterdodge(dodge.width = 0.8), alpha = 0.4) +
  stat_summary(
    aes(group = group),
    fun.data="mean_sdl",  fun.args = list(mult=1), 
    geom = "pointrange",  size = 0.2,
    position = position_dodge(0.8),
    color = "black"
  )+
  theme_ipsum(base_family = 'Arial', base_size = 12) +
  scale_color_manual(values=ecori_col, labels = c("M.EcoRI-treated", "Untreated"))+
  scale_x_discrete(labels=c("0" = "G", "1" = "A", "2" = "A", "3" = "T", "4" = "T", "5" = "C")) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, vjust = 0.5, hjust = 0.5),
    axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, color= c("black","black","#D62728FF","black", "black", "black")),
    legend.position = "bottom",
    legend.title = element_blank()
                               
  ) +
  labs(x = "position in GAATTC", 
       y = "methylation (%)")
 

figure11 <- ggarrange(p1j, p1_roc,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1, vjust=1)

ggsave(filename = paste(site,"_fig_1_15_02_2026.tiff", sep=""), plot =figure11, path = results_dir,
       scale = 1, width = 180,
       height = 78, units = c("mm"),
       dpi = 300, bg = "white")
ggsave(filename = paste(site,"_fig_1_15_02_2026.png", sep=""), plot =figure11, path = results_dir,
       scale = 1, width = 180,
       height = 78, units = c("mm"),
       dpi = 300, bg = "white")

