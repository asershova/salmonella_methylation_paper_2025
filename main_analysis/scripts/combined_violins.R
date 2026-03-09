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
library("viridis")
library(scales)
library(rgl)
library(data.table)
library("ggpubr")
library(extrafont)

working_dir = "/app/data"
results_dir = "/app/results"

gatc_violin <- readRDS(paste(working_dir, "GATC_violin_ggplot.rds", sep = "/"))
atgcat_violin <- readRDS(paste(working_dir, "ATGCAT_violin_ggplot.rds", sep = "/"))
ccwgg_violin <- readRDS(paste(working_dir, "CCWGG_violin_ggplot.rds", sep = "/"))

figure <- ggarrange(gatc_violin, ccwgg_violin, atgcat_violin,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1, vjust=1)
figure
annotated_figure <- annotate_figure(
  figure,
  top = NULL,
  bottom = NULL,
  left = "Log2FC",
  right = NULL,
  fig.lab = NULL,
  fig.lab.pos = c("top.left", "top", "top.right", "bottom.left", "bottom",
                  "bottom.right")
)
ggsave(filename = paste("All_fig_6.tiff", sep=""), plot =annotated_figure, 
       path = results_dir,
       scale = 1, width = 180,
       height = 78, units = c("mm"),
       dpi = 300, bg = "white")
ggsave(filename = paste("All_fig_6.png", sep=""), plot =annotated_figure, 
       path = results_dir,
       scale = 1, width = 180,
       height = 78, units = c("mm"),
       dpi = 300, bg = "white")


