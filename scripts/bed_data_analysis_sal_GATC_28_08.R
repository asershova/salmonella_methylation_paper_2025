require(extrafont)
    # need only do this once!
font_import(pattern="[A/a]rial", prompt=FALSE)
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
site="GATC"
sites_genome = 39428
meth_type = 'a'
paper_pal <-paletteer::paletteer_d("ggsci::category20_d3")
paper_pal
status_pal <- paper_pal[1:4]
violin_color <- c("#C49C94FF","#F7B6D2FF","#C7C7C7FF")
chi_test_color <- c("#7F7F7FFF","#C7C7C7FF")
status_pal_2 <- c("#AEC7E8FF","#FFBB78FF","#98DF8AFF","#FF9896FF")
data_08 <- read.csv(paste(working_dir,"LSP-1.GATC.bed", sep="/"), sep="", header = FALSE,
                    col.names = c('chrom', 'm_start', 'm_stop', 'met_type', 'total_reads', 'met_dir', 'met_st_1', 'met_st_2', 'color','reads_all', 'methyl', 'methyl_reads', 'non_methyl_reads',
                                  'v1','v2','v3', 'v4', 'v5', 'chrom2', 'site_start', 'site_stop', 'site', 'v6', 'dir',
                                  'chrom3', 'f_start', 'f_stop', 'f_dir','f_type','locus_tag','id1','id2','locus_tag1','id11','id21','d_chrom','d_start', 'd_stop','d_old','d_id','domain'))
data_09 <- read.csv(paste(working_dir,"MEP-1.GATC.bed", sep="/"), sep="", header = FALSE,
                    col.names = c('chrom', 'm_start', 'm_stop', 'met_type', 'total_reads', 'met_dir', 'met_st_1', 'met_st_2', 'color','reads_all', 'methyl', 'methyl_reads', 'non_methyl_reads',
                                  'v1','v2','v3', 'v4', 'v5', 'chrom2', 'site_start', 'site_stop', 'site', 'v6', 'dir',
                                  'chrom3', 'f_start', 'f_stop', 'f_dir','f_type','locus_tag','id1','id2','locus_tag1','id11','id21','d_chrom','d_start', 'd_stop','d_old','d_id','domain'))
data_10 <- read.csv(paste(working_dir,"MEP-2.GATC.bed", sep="/"), sep="", header = FALSE,
                    col.names = c('chrom', 'm_start', 'm_stop', 'met_type', 'total_reads', 'met_dir', 'met_st_1', 'met_st_2', 'color','reads_all', 'methyl', 'methyl_reads', 'non_methyl_reads',
                                  'v1','v2','v3', 'v4', 'v5', 'chrom2', 'site_start', 'site_stop', 'site', 'v6', 'dir',
                                  'chrom3', 'f_start', 'f_stop', 'f_dir','f_type','locus_tag','id1','id2','locus_tag1','id11','id21','d_chrom','d_start', 'd_stop','d_old','d_id','domain'))

data_12 <- read.csv(paste(working_dir,"LSP-2.GATC.bed", sep="/"), sep="", header = FALSE,
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
site_jitter<-ggplot(reordered_data[reordered_data$total_reads>2&reordered_data$met_type==meth_type,], 
                      aes(x = factor(rel_pos), y = methyl, color = sample)) + 
  geom_jitter(position=position_jitterdodge(dodge.width = 0.7), alpha = 0.4) +
  labs(x = paste("position in ",site, sep = ""), y = "methylation (%)") +
  scale_color_manual(values=sample_pal)+
  scale_x_discrete(labels=c("0" = "G", "1" = "A", "2" = "T", "3" = "C")) +
  theme_ipsum(base_family = 'Arial', base_size = 12) +
  theme(
    axis.title.x = element_blank(), #text(size = 16, vjust = 0.5, hjust = 0.5),
    axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, color= c("black", "#D62728FF", "black", "black")),
    axis.title.y = element_text(size = 12, vjust = 0.5, hjust = 0.5),
    plot.margin=unit(c(10,10,10,10), "pt"),
    legend.position = "right",
    legend.title = element_blank()
  )

site_jitter
ggsave(filename = paste(site,"_jitter.tiff", sep=""), plot =site_jitter, path = results_dir,
       scale = 1, width = 89,
       height = 78, units = c("mm"),
       dpi = 300, bg = "white")

site_jitter_total_reads<-ggplot(reordered_data[reordered_data$total_reads>1&reordered_data$met_type==meth_type,], 
                    aes(x = factor(rel_pos), y = total_reads, color = sample)) + 
  geom_jitter(position=position_jitterdodge(dodge.width = 0.7), alpha = 0.4) +
  labs(x = paste("position in ",site, sep = ""), y = "total reads") +
  scale_x_discrete(labels=c("0" = "G", "1" = "A", "2" = "T", "3" = "C")) +
  theme_ipsum(base_family = 'Arial', base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, color= c("black", "#D62728FF", "black", "black")),
    axis.title.y = element_text(size = 12, vjust = 0.5, hjust = 0.5),
    plot.margin=unit(c(10,10,10,10), "pt"),
  )

site_jitter_total_reads

ggsave(filename = paste(site,"_jitter_total_reads.tiff", sep=""), plot =site_jitter, path = results_dir,
       scale = 1, width = 89,
       height = 78, units = c("mm"),
       dpi = 300, bg = "white")

meth_pos = 1
data_all_methyl_a <- reordered_data[reordered_data$rel_pos==meth_pos&reordered_data$met_type==meth_type,]
#select unique positions - remove duplications due to annotations
uniq_sites <- unique(data_all_methyl_a %>% select(c(pos_id, chrom, methyl, site_start, sample, met_dir)))

data_summary_uniq_methyl <- uniq_sites[uniq_sites$methyl>=35,] %>% dplyr::count(sample)
data_summary_uniq_all_chrom <- uniq_sites %>% dplyr::count(sample,chrom)
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
#calculation of the weightened averages
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
r_2
heat_all <-ggplot(data_all_methyl_a_samples_pv_wide_wa, aes(MEP_wt_mean,LSP_wt_mean)) + 
  stat_bin2d(bins=32) + 
  ylim(0,100) +
  scale_fill_gradientn(colours=r_2, breaks = c(1,1000, 2000),labels=c(1,1000, 2000)) +
  labs(x = "MEP, methylation (%)",
       y = "LSP, methylation (%)")+
  theme_ipsum(base_family = 'Arial', base_size = 12) +
  theme(
    axis.title.y = element_text(size = 12, vjust = 0.8, hjust = 0.5),
    axis.title.x = element_text(size = 12, vjust = 0.8, hjust = 0.5),
    axis.text.x=element_text(size = 12),
    plot.margin=unit(c(10,10,10,10), "pt")
    
  ) +
  geom_vline(xintercept = 35) +
  geom_hline(yintercept = 35)

heat_all
ggsave(filename = paste(site,"_heat_all_ppt.tiff",sep=""), plot =heat_all, path = results_dir,
       scale = 1, width = 89,
       height = 78, units = c("mm"),
       dpi = 300, bg = "white")

data_all_methyl_a_samples_pv_wide_wa <- data_all_methyl_a_samples_pv_wide_wa %>% 
  mutate(feature = ifelse(f_type=="intergenic", 
                          "intergenic","intragenic")) 

#file for transcriptomic analysis
write.table(data_all_methyl_a_samples_pv_wide_wa, file = paste(results_dir,"/",site,"_data_all_flat_table_weight_av.csv",sep=""),
            sep = "\t", row.names = F)




#select unique positions for methyl status/feature barchart
#there are positions belonging to both, intergenic and gene bodies. In this case, we count them as intergenic.

unique_methyl_status_feature_domain_inter <- unique(select(data_all_methyl_a_samples_pv_wide_wa[data_all_methyl_a_samples_pv_wide_wa$feature=="intergenic",], 
                                                           c(pos_id, methyl_status, feature, dna, domain)))
unique_methyl_status_feature_domain_gb <- unique(select(data_all_methyl_a_samples_pv_wide_wa[!(data_all_methyl_a_samples_pv_wide_wa$pos_id%in%
                                                                                               unique_methyl_status_feature_domain_inter$pos_id),], 
                                                           c(pos_id, methyl_status, feature, dna, domain)))
unique_methyl_status_feature_domain<-rbind(unique_methyl_status_feature_domain_inter,unique_methyl_status_feature_domain_gb)
summary_test <- unique_methyl_status_feature_domain %>%
  group_by(pos_id) %>%
  summarize(n=n(),.groups="keep")
duplicates <- summary_test[summary_test$n>1,]
data_all_methyl_a_samples_pv_wide_wa[data_all_methyl_a_samples_pv_wide_wa$pos_id=="ST4-74.fa_1046346_+",]
summary_1 <- unique_methyl_status_feature_domain %>%
  group_by(methyl_status) %>%
  summarize(n = n(), .groups="keep")
summary_1_1 <-  unique_methyl_status_feature_domain %>%
  group_by(methyl_status, feature) %>%
  summarize(n = n(), .groups="keep")
summary_1_1
summary_1_2 <-  unique_methyl_status_feature_domain %>%
  group_by(methyl_status, dna, feature) %>%
  summarize(n = n(), .groups="keep")
summary_1_2
label = c()
for (x in 1:length(summary_1$methyl_status)){
  label<-append(label,(paste(summary_1$methyl_status[x],"\nn=",summary_1$n[x], sep="")))
}

summary_1_1
x_cu<-matrix((summary_1_1[summary_1_1$methyl_status=="Constitutive"|
                     summary_1_1$methyl_status=="UM",]$n), ncol=2)
x_cl<-matrix((summary_1_1[summary_1_1$methyl_status=="Constitutive"|
                            summary_1_1$methyl_status=="LSP-specific",]$n), ncol=2)
x_cm<-matrix((summary_1_1[summary_1_1$methyl_status=="Constitutive"|
                            summary_1_1$methyl_status=="MEP-specific",]$n), ncol=2)
p_chi_cl<-signif(fisher.test(x=x_cl)$p.value, digits = 3)
p_chi_cm<-signif(fisher.test(x=x_cm)$p.value, digits=3)
p_chi_cu<-signif(fisher.test(x=x_cu)$p.value,digits=3)
my_comparisons = list(c("Constitutive", "LSP-specific"),
                      c("Constitutive", "MEP-specific"),
                      c("Constitutive", "UM"))
bar_color <- violin_color[1:2]

stacked_bar<-ggplot(summary_1_1, 
                    aes(fill=feature, y=n, x=methyl_status)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = bar_color) +
  labs(x = "methylation status", y = paste(site," sites (%)",sep=""), fill = "feature type") +
  ggsignif::geom_signif(annotations = c(p_chi_cl,p_chi_cm,p_chi_cu), 
                        comparisons=my_comparisons, tip_length=0.000001, y_position=c(-1835.7, -1835.5, -1835.3),
                        textsize=3, size=0.2)+
  scale_x_discrete(labels=label)+
  scale_y_continuous(labels = scales::label_percent(suffix="")) +#scales::percent) +
  coord_cartesian(ylim = c(0, 1.8)) +
  theme_ipsum(base_family = 'Arial', base_size = 12) +
  theme(
    axis.title.y = element_text(size = 12, vjust = 0.8, hjust = 0.5, margin = margin(0, 10, 0, 0)),
    axis.title.x = element_blank(), #text(size = 14, vjust = 0.8, hjust = 0.5),
    axis.text.x=element_text(angle = 90, size = 12, vjust = 0.4, hjust = 0.4),
    axis.text.y=element_text(size = 12, vjust = 0.4, hjust = 0.4),
    plot.margin=unit(c(20,10,10,10), "pt"),
    legend.spacing = unit(1, "mm"),
    legend.position = "right",
    legend.title = element_blank()
    
  )
stacked_bar
ggsave(filename = paste(site,"_stacked_bar_all.tiff",sep=""), plot =stacked_bar, path = results_dir,
       scale = 1, width = 89,
       height = 78, units = c("mm"), 
       dpi = 300, bg = "white")
 
unique_methyl_status_feature_domain[unique_methyl_status_feature_domain$methyl_status=="LSP-specific",]
#macrodomain



gatc_pal<-status_pal_2[2:4]
gatc_pal
unique_methyl_status_feature_domain_summary_2 <- unique_methyl_status_feature_domain[unique_methyl_status_feature_domain$dna!="plasmid",] %>%
  group_by(methyl_status,domain) %>%
  summarize(n = n(), .groups="keep")
summary_2<- unique_methyl_status_feature_domain[unique_methyl_status_feature_domain$dna=="chromosome"&
                                                  unique_methyl_status_feature_domain$methyl_status!="Constitutive",] %>%
  group_by(domain) %>%
  summarize(n = n(), .groups="keep")
summary_2<-summary_2%>%arrange(factor(domain,levels=c("Origin","NSR","Right", "Ter", "Left", "NSL")))
summary_3<- data_all_methyl_a_samples_pv_wide_wa[data_all_methyl_a_samples_pv_wide_wa$dna!="plasmid",] %>%
  group_by(domain) %>%
  summarize(n = n(), .groups="keep")
summary_3
label_2 = c()
for (x in 1:length(summary_2$domain)){
  label_2<-append(label_2,(paste(summary_2$domain[x],", n=",summary_2$n[x], sep="")))
}

stacked_domains <- ggplot(unique_methyl_status_feature_domain_summary_2[unique_methyl_status_feature_domain_summary_2$methyl_status!=
                                                                          "Constitutive",], 
                          aes(fill=methyl_status, y=n, x= factor(domain,levels=c("Origin","NSR","Right", "Ter", "Left", "NSL")))) + 
  geom_bar(position="fill", stat="identity")+
  labs(x = "macrodomain", y = paste(site," sites (%)", sep=""), fill=' methylation status') +
  scale_x_discrete(labels=label_2)+
  scale_y_continuous(labels = scales::label_percent(suffix="")) +#scales::percent) +
  scale_fill_manual(values=gatc_pal) +
  theme_ipsum(base_family = 'Arial', base_size = 12) +
  theme(
    axis.title.x = element_blank(), #text(size = 14, vjust = 0.8, hjust = 0.5),
    axis.title.y = element_text(size = 12, vjust = 0.8, hjust = 0.5, margin = margin(0, 10, 0, 0)),
    plot.margin=unit(c(15,20,10,10), "pt"),
    axis.text.x=element_text(angle = 90, size = 12, vjust = 0.4, hjust = 0.4),
    legend.position = "right",
    legend.spacing = unit(1, "mm"),
    legend.title = element_blank()
  )
stacked_domains
gatc_pal
ggsave(filename = paste(site,"_dom_stacked_bar_chrom.tiff", sep=""), plot =stacked_domains, path = results_dir,
       scale = 1, width = 89,
       height = 78, units = c("mm"),
       dpi = 300, bg = "white")
figure <- ggarrange(site_jitter, heat_all, stacked_bar,stacked_domains,
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2, vjust=1)
figure
dev.off()
ggsave(filename = paste(site,"_figure_d3col_2.tiff", sep=""), plot =figure, path = results_dir,
       scale = 1, width = 180,
       height = 156, units = c("mm"),
       dpi = 300, bg = "white")
ggsave(filename = paste(site,"_figure_2.png", sep=""), plot =figure, path = results_dir,
       scale = 1, width = 180,
       height = 156, units = c("mm"),
       dpi = 300, bg = "white")

#feature analysis
data_all_methyl_a_samples_pv_wide_wa$feature_id <- paste(data_all_methyl_a_samples_pv_wide_wa$chrom, 
                                                         data_all_methyl_a_samples_pv_wide_wa$f_start,
                                                         data_all_methyl_a_samples_pv_wide_wa$f_stop, sep = "_")
summary_inter<- data_all_methyl_a_samples_pv_wide_wa[data_all_methyl_a_samples_pv_wide_wa$feature=="intergenic",] %>%
  group_by(feature_id) %>%
  summarize(n = n(), .groups="keep")
ggplot(summary_inter, aes(x=n)) +
  geom_density(color="#C7C7C7FF") +
  xlim(0,12)
feature_id_underm <- unique(data_all_methyl_a_samples_pv_wide_wa[data_all_methyl_a_samples_pv_wide_wa$methyl_status!="Constitutive"&
                                                            data_all_methyl_a_samples_pv_wide_wa$feature=="intergenic",]$feature_id)
feature_id_underm
data_all_methyl_a_samples_pv_wide_wa_underm <- data_all_methyl_a_samples_pv_wide_wa[data_all_methyl_a_samples_pv_wide_wa$feature_id%in%feature_id_underm,]
summary_inter_under<- data_all_methyl_a_samples_pv_wide_wa_underm %>%
  group_by(feature_id, methyl_status) %>%
  summarize(n = n(), .groups="keep")
summary_inter_under_wide <- summary_inter_under %>% 
  pivot_wider(names_from = methyl_status, 
              values_from = n)
summary_inter_under_wide <- summary_inter_under_wide %>% replace(is.na(.), 0)
summary_inter_under_wide[summary_inter_under_wide$UM>0&
                           summary_inter_under_wide$Constitutive==0,]
ggplot(summary_inter_under, aes(x=n)) +
  geom_density(color="#C7C7C7FF")
data_all_methyl_a_samples_pv_wide_wa_underm_heatmap <- select(data_all_methyl_a_samples_pv_wide_wa_underm, 
                                                              c("feature_id","locus_tag","id2","f_dir", "locus_tag1","id21","m_start","MEP_wt_mean","LSP_wt_mean","methyl_status", "dna"))
write.table(data_all_methyl_a_samples_pv_wide_wa_underm_heatmap, file = paste(results_dir, "/",site,"_UnderM_intergenic.csv",sep=""),
            sep = "\t", row.names = F)
write.table(summary_inter_under_wide, file = paste(results_dir,"/",site,"_UnderM_intergenic_pv.csv",sep=""),
            sep = "\t", row.names = F)
