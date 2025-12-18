#install.packages("reshape")
require(extrafont)
extrafont::loadfonts(device="all")
library(ggplot2)
library(hrbrthemes)
library(tidyverse)
library(readxl)
library(here)
library(skimr) # install.packages('skimr')
library(kableExtra) # install.packages('kableExtra')
library(ggsci)
library(reshape)
library(dplyr)
library(ggsci)
library("GGally")
library(data.table)
library(stringr)
library("ggpubr")

working_dir = "/app/data" # the path to the working directory, you should put your own
results_dir = "/app/results"

#open an example file
gatc_120 <- read.csv(paste(working_dir,"LSP-1.GATC.only.1.bed", sep = "/"),sep="", header = FALSE,
                     col.names = c('chrom', 'm_start', 'm_stop', 'met_type', 'total_reads', 'met_dir', 'met_st_1', 'met_st_2', 'color','reads_all', 'methyl', 'methyl_reads', 'non_methyl_reads',
                                   'v1','v2','v3', 'v4', 'v5', 'chrom-2', 'site_start', 'site_stop', 'site', 'v6', 'dir'))
gatc_120$sample <- "LSP-1.GATC.only.bed"
gatc_120$group <- "control"
#create a list of the files from your target directory
final_data <-data.frame()
for (k in c("5","10","20","30","40","50","60","70","80","90")){
  #setwd(paste(working_dir,"/coverage_",k, sep=""))
  file_list <- list.files(path=paste(working_dir,"/coverage_",k, sep=""))
  print(paste(working_dir,"/coverage_",k, sep=""))
  
  #initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable
  dataset <- data.frame()
  
  #had to specify columns to get rid of the total column
  for (i in 1:length(file_list)){
    temp_data <- read.csv(paste(working_dir, "/coverage_", k, "/",file_list[i],sep = ""),sep="", header = FALSE,
                          col.names = c('chrom', 'm_start', 'm_stop', 'met_type', 'total_reads', 'met_dir', 'met_st_1', 'met_st_2', 'color','reads_all', 'methyl', 'methyl_reads', 'non_methyl_reads','v1','v2','v3', 'v4', 'v5', 'chrom-2', 'site_start', 'site_stop', 'site', 'v6', 'dir'))#read in files using the fread function from the data.table package
    temp_data$sample <- file_list[i]
    temp_data$group<-"test"
    dataset <- rbindlist(list(dataset, temp_data), use.names = T, fill=TRUE) #for each iteration, bind the new data to the building dataset
  }
  
  dataset_all <- rbind(dataset, gatc_120)
  dataset_all$rel_pos <- dataset_all$m_start-dataset_all$site_start
  dataset_all <- dataset_all %>% 
    mutate(rel_pos = ifelse(met_dir=="+", 
                            dataset_all$m_start-dataset_all$site_start, 
                            dataset_all$site_stop-1-dataset_all$m_start))
  dataset_all_plus <- dataset_all[dataset_all$rel_pos==1&dataset_all$met_dir=="+"&dataset_all$dir=="+",]
  dataset_all_minus <- dataset_all[dataset_all$rel_pos==1&dataset_all$met_dir=="-"&dataset_all$dir=="-",]
  dataset_all_one<-rbind(dataset_all_plus, dataset_all_minus)
  dataset_plus_gg <- dataset_all_one %>%
    group_by(sample,chrom,m_start) %>%
    summarize(mean_methyl = mean(methyl), .groups="keep")
  
  dataset_plus_summary_pv_wide <-dataset_plus_gg %>% 
    pivot_wider(names_from = sample, 
                values_from = mean_methyl)
  dataset_plus_summary_pv_wide[is.na(dataset_plus_summary_pv_wide)] <- 0
  dif_methyl <- (dataset_plus_summary_pv_wide[,3:12]-dataset_plus_summary_pv_wide$`LSP-1.GATC.only.bed`)^2
  
  list_files <- colnames(dif_methyl)
  dif_methyl$sd <- sqrt(rowSums(dif_methyl[, list_files])/10)
  dif_methyl$mean <- dataset_plus_summary_pv_wide$`LSP-1.GATC.only.bed`
  dif_methyl$bins <- cut(dif_methyl$mean, breaks = c(-0.1,10,20,30,40,50,60,70,80,90,100))
  
  dif_methyl_pv_bins <- dif_methyl %>%
    group_by(bins) %>%
    summarize(mean_methyl = mean(mean), mean_sd = mean(sd), .groups="keep")
  dif_methyl_pv_bins$group <- paste(k,"x",sep='')
  final_data<- rbind(final_data, dif_methyl_pv_bins)  
}
p1<-ggplot(final_data) +
  geom_bar( aes(x=bins, y=mean_methyl), stat="identity", fill="forestgreen", alpha=0.5) +
  geom_errorbar( aes(x=bins, ymin=mean_methyl-mean_sd, ymax=mean_methyl+mean_sd), width=0.2, colour="orange", alpha=0.9,linewidth=0.5)+
  scale_x_discrete(labels=c("0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-100"))+
  facet_wrap(~factor(group, c("5x", "10x","20x","30x","40x","50x","60x","70x","80x","90x")),nrow=2)+
  labs(x = "Methylation in a control sample (%)", y="Methylation in subsamples (%)")+
  theme(axis.line = element_line(colour = "black", linewidth = .5),
        axis.title.y = element_text(size = 12, vjust = 0.8, hjust = 0.5),
        axis.title.x = element_text(size = 12, vjust = 0.8, hjust = 0.5),
        legend.position = "none",
        plot.margin=unit(c(30,20,10,10), "pt"),
        axis.text.x=element_text(size = 10, angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.y=element_text(size = 10, vjust = 0.5, hjust = 0.5),
        axis.ticks = element_line(color="black", linewidth = 0.8, linetype = 1),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey90"),
        axis.ticks.length = unit(.25, "cm"))

final_data$depth <- factor(final_data$group, levels = c("5x", "10x","20x","30x","40x","50x","60x","70x","80x","90x"))
p2<-ggplot(final_data, aes(x=bins, y=mean_sd, color=depth)) +
  geom_line(aes(group = depth, color = depth)) +
  scale_x_discrete(labels=c("0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-100")) +
  labs(x = "Methylation in a control sample (%)", y = "Mean deviation")+
  geom_point() +
  theme(axis.line = element_line(colour = "black", linewidth = .5),
        axis.title.y = element_text(size = 12, vjust = 0.8, hjust = 0.5),
        axis.title.x = element_text(size = 12, vjust = 0.8, hjust = 0.5),
        axis.text.x=element_text(size = 10, vjust = 0.5, hjust = 0.5),
        axis.text.y=element_text(size = 10, vjust = 0.5, hjust = 0.5),
        axis.ticks = element_line(color="black", linewidth = 0.8, linetype = 1),
        plot.margin=unit(c(15,20,10,10), "pt"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey90"),
        axis.ticks.length = unit(.25, "cm"))
figure <- ggarrange(p1, p2, nrow=2, labels = c("A", "B"), vjust=1,
          font.label = list(size = 12, face = "bold", color ="black", family = 'Arial'))

figure
ggsave(filename = paste("seq_depth_figure_S1.3.tiff", sep=""), plot =figure, path = results_dir,
       scale = 1, width = 180,
       height = 156, units = c("mm"),
       dpi = 600, bg = "white")
ggsave(filename = paste("seq_depth_figure_S1.3.png", sep=""), plot =figure, path = results_dir,
       scale = 1, width = 180,
       height = 156, units = c("mm"),
       dpi = 600, bg = "white")
