require(extrafont)
extrafont::loadfonts(device="all")

library(ggplot2)
library(dplyr)
library(tidyverse)
library(readxl)
library(here)
library(skimr) # install.packages('skimr')
library(kableExtra) # install.packages('kableExtra')
library(ggsci)
library("GGally")
library(RColorBrewer)
library(hrbrthemes)
library("ggpubr")
working_dir = "/app/data" # the path to the working directory, you should put your own
results_dir = "/app/results"

paper_pal <-paletteer::paletteer_d("ggsci::category20_d3")

status_pal <- paper_pal[1:4]
violin_color <- c("#C49C94FF","#F7B6D2FF","#C7C7C7FF")
chi_test_color <- c("#7F7F7FFF","#C7C7C7FF")
status_pal_2 <- c("#AEC7E8FF","#FFBB78FF","#98DF8AFF","#FF9896FF")

site="CAGAG"
sites_genome = 6125
meth_type = 'a'
data_08 <- read.csv(paste(working_dir,"LSP-1.CAGAG.bed", sep="/"), sep="", header = FALSE,
                    col.names = c('chrom', 'm_start', 'm_stop', 'met_type', 'total_reads', 'met_dir', 'met_st_1', 'met_st_2', 'color','reads_all', 'methyl', 'methyl_reads', 'non_methyl_reads',
                                  'v1','v2','v3', 'v4', 'v5', 'chrom2', 'site_start', 'site_stop', 'site', 'v6', 'dir',
                                  'chrom3', 'f_start', 'f_stop', 'f_dir','f_type','locus_tag','id1','id2','locus_tag1','id11','id21','d_chrom','d_start', 'd_stop','d_old','d_id','domain'))
data_09 <- read.csv(paste(working_dir,"MEP-1.CAGAG.bed", sep="/"), sep="", header = FALSE,
                    col.names = c('chrom', 'm_start', 'm_stop', 'met_type', 'total_reads', 'met_dir', 'met_st_1', 'met_st_2', 'color','reads_all', 'methyl', 'methyl_reads', 'non_methyl_reads',
                                  'v1','v2','v3', 'v4', 'v5', 'chrom2', 'site_start', 'site_stop', 'site', 'v6', 'dir',
                                  'chrom3', 'f_start', 'f_stop', 'f_dir','f_type','locus_tag','id1','id2','locus_tag1','id11','id21','d_chrom','d_start', 'd_stop','d_old','d_id','domain'))
data_10 <- read.csv(paste(working_dir,"MEP-2.CAGAG.bed", sep="/"), sep="", header = FALSE,
                    col.names = c('chrom', 'm_start', 'm_stop', 'met_type', 'total_reads', 'met_dir', 'met_st_1', 'met_st_2', 'color','reads_all', 'methyl', 'methyl_reads', 'non_methyl_reads',
                                  'v1','v2','v3', 'v4', 'v5', 'chrom2', 'site_start', 'site_stop', 'site', 'v6', 'dir',
                                  'chrom3', 'f_start', 'f_stop', 'f_dir','f_type','locus_tag','id1','id2','locus_tag1','id11','id21','d_chrom','d_start', 'd_stop','d_old','d_id','domain'))
data_12 <- read.csv(paste(working_dir,"LSP-2.CAGAG.bed", sep="/"), sep="", header = FALSE,
                    col.names = c('chrom', 'm_start', 'm_stop', 'met_type', 'total_reads', 'met_dir', 'met_st_1', 'met_st_2', 'color','reads_all', 'methyl', 'methyl_reads', 'non_methyl_reads',
                                  'v1','v2','v3', 'v4', 'v5', 'chrom2', 'site_start', 'site_stop', 'site', 'v6', 'dir',
                                  'chrom3', 'f_start', 'f_stop', 'f_dir','f_type','locus_tag','id1','id2','locus_tag1','id11','id21','d_chrom','d_start', 'd_stop','d_old','d_id','domain'))

data_08$sample<- 'LSP-1'
data_09$sample<- 'MEP-1'
data_10$sample<- 'MEP-2'
data_12$sample<- 'LSP-2'
data_08$group<- 'LSP'
data_09$group<- 'MEP'
data_10$group<- 'MEP'
data_12$group<- 'LSP'
data_all <-rbind(data_08, data_09, data_10, data_12)
data_all <- data_all %>% 
  mutate(rel_pos = ifelse(met_dir=="+", 
                          data_all$m_start-data_all$site_start, 
                          data_all$site_stop-1-data_all$m_start))  
data_all$pos_id<- paste(data_all$chrom,data_all$m_start, data_all$met_dir, sep = '_')

data_summary_sample <- data_all[data_all$met_dir=="+"&data_all$dir=="+",] %>%
  group_by(sample,rel_pos) %>%
  summarize(mean_methyl = mean(methyl), reads = mean(total_reads), .groups="keep")
data_all_plus <- data_all[data_all$met_dir=="+"&data_all$dir=="+",]
data_all_minus <- data_all[data_all$met_dir=="-"&data_all$dir=="-",]
data_all_methyl <- rbind(data_all_plus, data_all_minus)

set.seed(42)
sample_pal <- c("#BCBD22FF","#DBDB8DFF","#17BECFFF","#9EDAE5FF")
reordered_data<-data_all_methyl %>%
  mutate(sample = fct_relevel(sample, 
                              "MEP-1", "MEP-2", "LSP-1", "LSP-2"))
site_jitter<-ggplot(reordered_data[reordered_data$total_reads>1&reordered_data$met_type==meth_type,], 
                      aes(x = factor(rel_pos), y = methyl, color = sample)) + 
  geom_jitter(position=position_jitterdodge(dodge.width = 0.7), alpha = 0.7) +
  labs(x = paste("position in ",site, sep = ""), y = "methylation (%)") +
  scale_color_manual(values=sample_pal)+
  scale_x_discrete(labels=c("0" = "C", "1" = "A", "2" = "G", "3" = "A", "4" = "G")) +
  theme_ipsum(base_family = 'Arial', base_size = 12) +
  theme(
    axis.title.x = element_blank(), #text(size = 12, vjust = 0.5, hjust = 0.5),
    axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, color= c("black","black","black","#D62728FF", "black")),
    axis.title.y = element_text(size = 12, vjust = 0.5, hjust = 0.5),
    plot.margin=unit(c(10,10,10,10), "pt"),
    legend.position = "right",
    legend.title = element_blank()
  )


ggsave(filename = paste(site,"_jitter.tiff", sep=""), plot =site_jitter, path = results_dir,
       scale = 1, width = 89,
       height = 78, units = c("mm"),
       dpi = 300, bg = "white")


meth_pos = 3
data_all_methyl_a <- reordered_data[reordered_data$rel_pos==meth_pos&reordered_data$met_type==meth_type,]
#select unique positions - remove duplications due to annotations
uniq_sites <- unique(data_all_methyl_a %>% select(c(pos_id, methyl, site_start, sample, met_dir)))

data_summary_uniq_methyl <- uniq_sites[uniq_sites$methyl>=35,] %>% dplyr::count(sample)
data_summary_uniq_all <- uniq_sites %>% dplyr::count(sample)
data_summary_uniq_methyl$count_all<-data_summary_uniq_all$n
data_summary_uniq_methyl$ratio <- data_summary_uniq_methyl$n/data_summary_uniq_all$n
data_summary_a <- reordered_data[reordered_data$chrom=="ST4-74.fa"&
                                   reordered_data$rel_pos==meth_pos&
                                   reordered_data$met_type==meth_type,] %>%
  group_by(sample) %>%
  summarize(mean_reads = mean(total_reads), .groups="keep")
data_summary_a_all_dna <- reordered_data[reordered_data$rel_pos==meth_pos&
                                   reordered_data$met_type==meth_type,] %>%
  group_by(sample, chrom) %>%
  summarize(mean_reads = mean(total_reads), .groups="keep")

data_all_methyl_a_chrom <- reordered_data[reordered_data$rel_pos==meth_pos & 
                                            reordered_data$chrom=='ST4-74.fa' &
                                            reordered_data$met_type==meth_type,]
data_all_methyl_a_chrom_pv <- data_all_methyl_a_chrom %>% group_by(sample) %>%  
  summarize(coverage = mean(total_reads), .groups = "keep")
data_all_methyl_a_plasmid <- reordered_data[reordered_data$rel_pos==meth_pos & 
                                              reordered_data$chrom!='ST4-74.fa' &
                                              reordered_data$met_type==meth_type,]

#weightened average

#https://datacornering.com/how-to-calculate-weighted-mean-in-r/
## weight is a read coverage from data_summary_a_all_dna
##calculate a pivot table with methylation data
data_all_methyl_a_samples_pv <- data_all_methyl_a %>%
  group_by(sample,chrom,pos_id,m_start,m_stop,met_dir,site_start, site_stop,f_type,locus_tag,id2,f_dir,locus_tag1,id21,f_start,f_stop,domain)%>%
  summarize(mean_methyl = mean(methyl), .groups = "keep")
data_all_methyl_a_samples_pv_wide <- data_all_methyl_a_samples_pv%>% 
  pivot_wider(names_from = sample, 
              values_from = mean_methyl)

#weight is total_reads
data_all_total_reads_a_samples_pv <- data_all_methyl_a %>%
  group_by(sample,pos_id)%>%
  summarize(mean_reads = mean(total_reads), .groups = "keep")
data_all_total_reads_a_samples_pv_wide <- data_all_total_reads_a_samples_pv%>% 
  pivot_wider(names_from = sample, 
              values_from = mean_reads)
#weightened average
data_all_methyl_a_samples_pv_wide_wa <-merge(data_all_methyl_a_samples_pv_wide,data_all_total_reads_a_samples_pv_wide, by="pos_id") 
data_all_methyl_a_samples_pv_wide_wa<-data_all_methyl_a_samples_pv_wide_wa %>%
  rowwise() %>%
  mutate("MEP_wt_mean" = weighted.mean(across(c(`MEP-1.x`,`MEP-2.x`)), across(c(`MEP-1.y`,`MEP-2.y`)),na.rm=TRUE),
         "LSP_wt_mean" = weighted.mean(across(c(`LSP-1.x`,`LSP-2.x`)), across(c(`LSP-1.y`,`LSP-2.y`)),na.rm=TRUE)) %>%
  as.data.frame()
rows_with_na <- data_all_methyl_a_samples_pv_wide_wa[!complete.cases(data_all_methyl_a_samples_pv_wide_wa), ]

#table for annotation
data_all_methyl_a_samples_pv_wide_wa <- data_all_methyl_a_samples_pv_wide_wa %>% 
  mutate(methyl_status = ifelse(MEP_wt_mean<35&LSP_wt_mean<35, 
                                "UM",ifelse(MEP_wt_mean>=35&LSP_wt_mean>=35, "Constitutive",
                                                      ifelse(MEP_wt_mean<35&LSP_wt_mean>=35, "LSP-specific",
                                                             "MEP-specific")))) 
data_all_methyl_a_samples_pv_wide_wa <- data_all_methyl_a_samples_pv_wide_wa %>% 
  mutate(dna = ifelse(chrom=="ST4-74.fa","chromosome","plasmid")) 
#annotated - change table
#chromosome
data_all_sample_pv_wide_anno_summary <-  data_all_methyl_a_samples_pv_wide_wa %>%
  group_by(methyl_status, f_type) %>%
  summarize(n = n(), .groups="keep")

#palette
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)
rf_2 <- colorRampPalette(c("#1F77B4FF","#98DF8AFF","#DBDB8DFF","#FFBB78FF","#D62728FF"))
r_2 <- rf_2(32)

heat_all <-ggplot(data_all_methyl_a_samples_pv_wide_wa, aes(MEP_wt_mean,LSP_wt_mean)) + 
  stat_bin2d(bins=32) + 
  ylim(0,110) +
  scale_fill_gradientn(colours=r_2, breaks = c(10,100,200,300),labels=c(10,100, 200,300)) +
  labs(x = "MEP, methylation (%)",
       y = "LSP, methylation (%)")+
  theme_ipsum(base_family = 'Arial', base_size = 12) +
  theme(
    axis.title.y = element_text(size = 12, vjust = +5, hjust = 0.5, margin=margin(t=0,r=10,b=0,l=10, unit="pt")),
    axis.title.x = element_text(size = 12, vjust = 0.8, hjust = 0.5),
    axis.text.x=element_text(size = 12),
    plot.margin=unit(c(15,15,15,25), "pt")
    
  ) +
  geom_vline(xintercept = 35) +
  geom_hline(yintercept = 35)


ggsave(filename = paste(site,"_heat_all.tiff",sep=""), plot =heat_all, path = results_dir,
       scale = 1, width = 89,
       height = 78, units = c("mm"),
       dpi = 300, bg = "white")

figure <- ggarrange(site_jitter, heat_all,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1, vjust=1)


ggsave(filename = paste(site,"_figure_S2.2.AB.tiff", sep=""), plot =figure, path = results_dir,
       scale = 1, width = 180,
       height = 78, units = c("mm"),
       dpi = 300, bg = "white")
ggsave(filename = paste(site,"_figure_S2.2.AB.png", sep=""), plot =figure, path = results_dir,
       scale = 1, width = 180,
       height = 78, units = c("mm"),
       dpi = 300, bg = "white")


