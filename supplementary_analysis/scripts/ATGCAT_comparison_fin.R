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
library("viridis")
library(scales)
library(rgl)
library(data.table)
library("ggpubr")
library("VennDiagram")
working_dir = "/app/data" # the path to the working directory, you should put your own
results_dir = "/app/results"

site='ATGCAT'
paper_pal <-paletteer::paletteer_d("ggsci::category20_d3")
status_pal <- paper_pal[1:4]
violin_color <- c("#C49C94FF","#F7B6D2FF","#C7C7C7FF")
chi_test_color <- c("#7F7F7FFF","#C7C7C7FF")
status_pal_2 <- c("#AEC7E8FF","#FFBB78FF","#98DF8AFF","#FF9896FF")

gatc_data <- read.csv(paste(working_dir,"ATGCAT_data_all_flat_table_weight_av.csv", sep="/"), sep="\t", header = TRUE)
my_data <- read.csv(paste(working_dir,"ST4-74_mydata_ATGCAT_liftover_strand.bed", sep="/"), sep="\t", header = FALSE)
GSE185578_data <- read.csv(paste(working_dir,"GSE185578_ATGCAT_liftover_strand.bed", sep="/"), sep="\t", header = FALSE)
gatc_um_intre <- unique(gatc_data[gatc_data$feature=="intergenic"&gatc_data$methyl_status=="UM",]$pos_id)
gatc_um_intra <- unique(gatc_data[gatc_data$feature=="intragenic"&gatc_data$methyl_status=="UM",]$pos_id)
#adding affected genes
my_data <- my_data %>% mutate(
  affected_gene = case_when(V38 == "+" ~ V41,
                             V38 == "-" ~ V41,
                             V38 == "++" ~ V43,
                             V38 == "--" ~ V41,
                             V38 == "+-" ~ paste(V41, V43, "3'end", sep=' '),
                             V38 == "-+" ~ paste(V41, V43, sep='/'))
)

my_data_sel <- unique(select(my_data, c("V1", "V2", "V3", "V6","V10", "V28","V29", "V30", "V32", "affected_gene")))
colnames(my_data_sel) <- c("chrom", "start", "stop", "dir","aln_id","MEP", "LSP", 
                           "methylation status", "feature type", "gene")
GSE_data_sel <- unique(select(GSE185578_data, c("V10", "V15")))
colnames(GSE_data_sel) <- c("aln_id", "PMC9239280")
combined_data <- merge(my_data_sel, GSE_data_sel, by = "aln_id", all.x = TRUE)
combined_data$PMC9239280p<-combined_data$PMC9239280*100
combined_data <- combined_data %>% replace(is.na(.), -1)
colnames(combined_data) <- c("aln_id","chrom", "start", "stop", "dir", "MEP", "LSP", 
                             "methylation status", "feature type", "gene", "PMC9239280","PMC9239280p" )

combined_data_pv <- combined_data%>%
  group_by(`feature type`)%>%
  summarise(n=n())

characterised <- combined_data[combined_data$PMC9239280p > -1,]
characterised_pv <- characterised%>%
  group_by(`feature type`)%>%
  summarise(n=n())
gg1<-ggplot(combined_data[combined_data$PMC9239280p > -1,], aes(x=LSP, y=PMC9239280p, 
                               colour = `methylation status`)) + 
  geom_point(alpha=0.7, size = 2.5) +
  labs(y = 'set2', x = 'LSP',
       color = "methylation status") +
  scale_colour_manual(values=status_pal_2) +
  theme_ipsum(base_family = 'Arial', base_size = 12) +
  guides(colour = guide_legend(nrow = 2)) +
  facet_wrap(~`feature type`) + 
  theme(
    axis.title.x = element_text(size = 12, vjust = 0.8, hjust = 0.5),
    axis.title.y = element_text(size = 12, vjust = 0.8, hjust = 0.5),
    plot.margin=unit(c(15,15,10,10), "pt"),
    axis.text.x=element_text(angle = 90, size = 12, vjust = 0.4, hjust = 0.4),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.spacing = unit(1, "mm")
  )

ggsave(filename = "PMC9239280_vs_my_data_LSP_ATGCAT_fig.S4.1_AB.png", plot =gg1, 
       path = results_dir,
       scale = 1, width = 178,
       height = 78, units = c("mm"),
       dpi = 600, bg = "white")
combined_data_um <- combined_data[combined_data$`PMC9239280,%`<25,]
combined_data_um_pv <- combined_data_um%>%
  group_by(`methylation status`)%>%
  summarise(n=n())
intragenic_data <- combined_data[combined_data$`feature type`=='intragenic',]
intragenic_data_MEP <- intragenic_data[intragenic_data$`methylation status`=="MEP-specific",]



###Venn diagram###

venn_list_inter <- list(combined_data[combined_data$PMC9239280p==0&
                                  combined_data$`feature type`=='intergenic',]$aln_id,
                  combined_data[combined_data$`methylation status`=='UM'&
                                  combined_data$PMC9239280p>-1&
                                combined_data$`feature type`=='intergenic',]$aln_id)
venn_list_intragenic <- list(combined_data[combined_data$PMC9239280p==0&
                                        combined_data$`feature type`=='intragenic',]$aln_id,
                        combined_data[combined_data$`methylation status`=='UM'&
                                        combined_data$PMC9239280p>-1&
                                        combined_data$`feature type`=='intragenic',]$aln_id)

venn.diagram(
  x = venn_list_inter,
  category.names = c("set2","This work"),
  filename = paste(results_dir,'PMC9239280_vs_my_data_UM_ATGCAT_venn_inter_figS4.1C.png',sep="/"),
  output = TRUE ,
  imagetype="png" ,
  height = 78,
  width = 89, 
  units=c("mm"),
  resolution = 600,
  compression = "lzw",
  lwd = 1,
  col=c('#7F7F7FFF', '#FF9896FF'),
  fill = c(alpha('#7F7F7FFF',0.5), alpha('#FF9896FF',0.5)),
  cex = .6,
  fontfamily = "Arial",
  cat.cex = .7,
  cat.default.pos = "outer",
  cat.pos = c(-27, 135),
  #cat.dist = c(0.055, 0.085),
  cat.fontfamily = "Arial"
  #cat.col = c("#2CA02CFF", '#7F7F7FFF', '#FF9896FF'),
  #rotation = 1
)
venn.diagram(
  x = venn_list_intragenic,
  category.names = c("set2","This work"),
  filename = paste(results_dir, 'PMC9239280_vs_my_data_UM_ATGCAT_venn_intra_fig_S4.1D.png', sep = "/"),
  output = TRUE ,
  imagetype="png" ,
  height = 78,
  width = 89, 
  units=c("mm"),
  resolution = 600,
  compression = "lzw",
  lwd = 1,
  col=c('#7F7F7FFF', '#FF9896FF'),
  fill = c(alpha('#7F7F7FFF',0.5), alpha('#FF9896FF',0.5)),
  cex = .6,
  fontfamily = "Arial",
  cat.cex = .7,
  cat.default.pos = "outer",
  cat.pos = c(-27, 135),
  cat.fontfamily = "Arial"
)
