#install.packages('skimr', 'here', 'kableExtra')
#install.packages("patchwork")
#install.packages("pROC")
require(extrafont)
extrafont::loadfonts(device="all")

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
working_dir = "/app/data" # the path to the working directory, you should put your own
results_dir = "/app/results"
site = "GAATTC"
col_names = c('chrom', 'm_start', 'm_stop', 'met_type', 'total_reads', 'met_dir', 'met_st_1', 'met_st_2', 'color','reads_all', 'methyl', 'methyl_reads', 'non_methyl_reads',
              'v1','v2','v3', 'v4', 'v5', 'chrom2', 'site_start', 'site_stop', 'site', 'v6', 'dir',
              'chrom3', 'f_start', 'f_stop', 'f_dir','f_type','locus_tag','id1','id2','locus_tag1','id11','id21','d_chrom','d_start', 'd_stop','d_old','d_id','domain')
#working_dir = "/home/anna/Documents/Methylome/results/experiments/nanopore_salmonella/sal_ent_barcodes_6mA/aligned/0_4_4_bed/mod_data" # the path to the working directory, you should put your own
data_08 <- read.csv(paste(working_dir,"WT.GAATTC.bed", sep="/"), sep="", header = FALSE,
                    col.names = col_names)
data_14 <- read.csv(paste(working_dir,"M.EcoRI.GAATTC.bed", sep="/"), sep="", header = FALSE,
                    col.names = col_names)


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
p1<- ggplot(reordered_data, 
            aes(x=group, y=methyl, fill=group)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin()+
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",
               colour = "black")+
  theme_ipsum(base_family = 'Arial', base_size = 12) +
  theme(
    axis.title.x = element_text(size = 14, vjust = 0.5, hjust = 0.5),
    axis.title.y = element_text(size = 14, vjust = 0.5, hjust = 0.5)
  )+labs(x = "", y = "methylation, %")
p2<- ggplot(reordered_data[reordered_data$chrom=='ST4-74.fa',], 
            aes(x=group, y=methyl, fill=group)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin()+ 
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",
               colour = "black") +
  theme_ipsum(base_family = 'Arial', base_size = 12) +
  theme(
    axis.title.x = element_text(size = 14, vjust = 0.5, hjust = 0.5),
    axis.title.y = element_text(size = 14, vjust = 0.5, hjust = 0.5)
  )+labs(x = "", y = "methylation, %")
p3<-ggplot(reordered_data[reordered_data$chrom!='ST4-74.fa',], 
           aes(x=group, y=methyl, fill=group)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin()+ 
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",
               colour = "black") +
  theme_ipsum(base_family = 'Arial', base_size = 12) +
  theme(
    axis.title.x = element_text(size = 14, vjust = 0.5, hjust = 0.5),
    axis.title.y = element_text(size = 14, vjust = 0.5, hjust = 0.5)
  )+labs(x = "", y = "methylation, %")


figure12 <- ggarrange(p2, p3, 
                      labels = c("A", "B"),
                      ncol = 2, nrow = 1, vjust=1, legend = "none")
figure12

ggsave(filename = paste(site,"_fig_S1.2.tiff", sep=""), plot =figure12, path = results_dir,
       scale = 1, width = 180,
       height = 156, units = c("mm"),
       dpi = 600, bg = "white")
ggsave(filename = paste(site,"_fig_S1.2.png", sep=""), plot =figure12, path = results_dir,
       scale = 1, width = 180,
       height = 156, units = c("mm"),
       dpi = 300, bg = "white")
#ROC curve and thresholds

roc_gaattc <- roc(data_all_methyl_a_ecori_2$gaattc_stat, data_all_methyl_a_ecori_2$methyl)

coords(roc_gaattc, "best", ret=c("threshold", "specificity", "1-npv"))
coords(roc=roc_gaattc, x=35, input="threshold")
#threshold specificity       1-npv
#threshold    34.535   0.9924433 0.003162555


p1_roc<-ggroc(roc_gaattc) +
  theme_ipsum(base_family = 'Arial', base_size = 12) +
  theme(
    axis.title.x = element_text(size = 12, vjust = 0.5, hjust = 0.5),
    axis.title.y = element_text(size = 12, vjust = 0.5, hjust = 0.5))
p1_roc
#Jitter plot
p1j <- ggplot(data_all_methyl_a_ecori[data_all_methyl_a_ecori$total_reads>1,], 
       aes(x = factor(rel_pos), y = methyl, color = group)) + 
  geom_jitter(position=position_jitterdodge(dodge.width = 1.0), alpha = 0.4) +
  theme_ipsum(base_family = 'Arial', base_size = 12) +
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
p1j 
p2j <- ggplot(data_all_methyl_a_ecori[data_all_methyl_a_ecori$total_reads>1,], 
              aes(x = rel_pos, y = total_reads, color = sample)) + 
  geom_jitter(position=position_jitterdodge(dodge.width = 1.0), alpha = 0.4) +
  theme_ipsum(base_family = 'Arial', base_size = 12) +
  theme(
    axis.title.x = element_text(size = 14, vjust = 0.5, hjust = 0.5),
    axis.title.y = element_text(size = 14, vjust = 0.5, hjust = 0.5)
  ) +
  labs(x = "position in GAATTC", 
       y = "methylated reads")
p2j 
figure11 <- ggarrange(p1j, p1_roc, 
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1, vjust=1)
figure11

ggsave(filename = paste(site,"_fig_S1.1.tiff", sep=""), plot =figure11, path = results_dir,
       scale = 1, width = 180,
       height = 156, units = c("mm"),
       dpi = 600, bg = "white")
ggsave(filename = paste(site,"_fig_S1.1.png", sep=""), plot =figure11, path = results_dir,
       scale = 1, width = 180,
       height = 156, units = c("mm"),
       dpi = 300, bg = "white")


