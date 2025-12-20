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
library("viridis")
library(scales)
library(rgl)
library(data.table)
library("ggpubr")

working_dir = "/app/data" # the path to the working directory, you should put your own
results_dir = "/app/results"

paper_pal <-paletteer::paletteer_d("ggsci::category20_d3")
status_pal <- paper_pal[1:4]
violin_color <- c("#C49C94FF","#F7B6D2FF","#C7C7C7FF")
chi_test_color <- c("#7F7F7FFF","#C7C7C7FF")
status_pal_2 <- c("#AEC7E8FF","#FFBB78FF","#98DF8AFF","#FF9896FF")

gatc_data <- read.csv(paste(working_dir,"GATC_data_all_flat_table_weight_av.csv", sep="/"), sep="\t", header = TRUE)
my_data <- read.csv(paste(working_dir,"ST4-74_mydata_GATC_liftover_strand.bed", sep="/"), sep="\t", header = FALSE)
GSE185578_data <- read.csv(paste(working_dir,"GSE185578_GATC_19_10_25_liftover_strand.bed", sep="/"), sep="\t", header = FALSE)
data_casadeus <- read.csv(paste(working_dir,"casadeus_undermethyl_liftover_strand.bed", sep = "/"), sep="\t", header =FALSE)
colnames(data_casadeus) <- c("chrom", 'site_start', 'site_stop','site', 'rel_num', 'site_dir',
                             'rel_start', 'rel_stop', 'cp_id', 'aln_id', 'chrom2', 'm_start', 
                             'm_stop', 'pos_abs', 'gene', 'm_dir')
data_casadeus$value <- 100
#adding affected genes
my_data <- my_data %>% mutate(
  affected_gene = case_when(V19 == "+" ~ V18,
                             V19 == "-" ~ V18,
                             V19 == "++" ~ V21,
                             V19 == "--" ~ V18,
                             V19 == "+-" ~ paste(V18, V21, "3'end", sep=' '),
                             V19 == "-+" ~ paste(V18, V21, sep='/'))
)

my_data_sel <- unique(select(my_data, c("V1", "V2", "V3", "V6","V10", "V33","V34", "V35", "V37", "affected_gene")))
colnames(my_data_sel) <- c("chrom", "start", "stop", "dir","aln_id","MEP", "LSP", 
                           "methylation status", "feature type", "gene")
GSE_data_sel <- unique(select(GSE185578_data, c("V10", "V20")))
colnames(GSE_data_sel) <- c("aln_id", "PMC9239280")
combined_data <- merge(my_data_sel, GSE_data_sel, by = "aln_id")
combined_data$PMC9239280p<-combined_data$PMC9239280*100
colnames(combined_data) <- c("aln_id","chrom", "start", "stop", "dir", "MEP", "LSP", 
                             "methylation status", "feature type", "gene", "PMC9239280","PMC9239280,%" )
combined_data_casadeus <- merge(my_data_sel, data_casadeus, by = "aln_id")
my_data_id_um <-my_data_sel[my_data_sel$`methylation status`=="UM",]$aln_id
GSE185578_data_zero_id <- GSE_data_sel[GSE_data_sel$PMC9239280==0,]$aln_id
data_casadeus_id <- data_casadeus$aln_id
undermethyl_id <- unique(c(my_data_id_um, GSE185578_data_zero_id, data_casadeus_id))
my_data_all_under <- my_data_sel[my_data_sel$aln_id%in%undermethyl_id,]
all_under_table <- merge(my_data_all_under,GSE_data_sel, by = "aln_id", all.x = TRUE)
all_under_table <- all_under_table %>% mutate(casadeus = ifelse
                                              (all_under_table$aln_id%in%data_casadeus$aln_id, 0, NA))
#all_under_table_intergenic <- all_under_table[all_under_table$`feature type`=='intergenic',]

all_under_table$PMC9239280p <- all_under_table$PMC9239280*100
all_under_table_score <- all_under_table %>% replace(is.na(.), -1)
all_under_table_score<-all_under_table_score%>%
  mutate(venn_category = case_when(`methylation status`=='UM'& casadeus==0 & PMC9239280p==0 ~ 'g',
                                   `methylation status`!='UM'& casadeus==0 & PMC9239280p!=0 ~ 'a',
                                   `methylation status`!='UM'& casadeus!=0 & PMC9239280p==0 ~ 'b',
                                   `methylation status`=='UM'& casadeus!=0 & PMC9239280p!=0 ~ 'c',
                                   `methylation status`=='UM'& casadeus==0 & PMC9239280p!=0 ~ 'd',
                                   `methylation status`!='UM'& casadeus==0 & PMC9239280p==0 ~ 'e',
                                   `methylation status`=='UM'& casadeus!=0 & PMC9239280p==0 ~ 'f',
                                   ))
all_under_table_score_only <- select(all_under_table_score, c('aln_id', 'venn_category'))
all_under_table_methyl_data <- select(all_under_table, c('MEP', 'LSP', 'PMC9239280p', 'casadeus'))
all_under_table_methyl_data_bed <- select(all_under_table, c('chrom', 'start', 'stop', 'dir','MEP', 'LSP', 'PMC9239280p', 'casadeus'))
all_under_table_flat_mep <- select(all_under_table, c('chrom', 'start', 'dir','gene','MEP','feature type', 'aln_id'))
all_under_table_flat_mep$sample <- 'MEP'
colnames(all_under_table_flat_mep) <- c('chrom', 'start', 'dir','gene','methyl','feature type', 'aln_id', 'sample')

all_under_table_flat_lsp <- select(all_under_table, c('chrom', 'start', 'dir','gene','LSP','feature type', 'aln_id'))
all_under_table_flat_lsp$sample <- 'LSP'
colnames(all_under_table_flat_lsp) <- c('chrom', 'start', 'dir','gene','methyl','feature type', 'aln_id', 'sample')
all_under_table_flat_PMC9239280 <- select(all_under_table, c('chrom', 'start', 'dir','gene','PMC9239280p','feature type', 'aln_id'))
all_under_table_flat_PMC9239280$sample <- 'PMC9239280'
colnames(all_under_table_flat_PMC9239280) <- c('chrom', 'start', 'dir','gene','methyl','feature type', 'aln_id','sample')
all_under_table_flat_casadeus <- select(all_under_table, c('chrom', 'start', 'dir','gene','casadeus','feature type', 'aln_id'))
all_under_table_flat_casadeus$sample <- 'PMC7708049'
colnames(all_under_table_flat_casadeus) <- c('chrom', 'start', 'dir','gene','methyl','feature type', 'aln_id', 'sample')
all_under_table_flat <- rbind(all_under_table_flat_mep, all_under_table_flat_lsp, 
                              all_under_table_flat_PMC9239280, all_under_table_flat_casadeus)
all_under_table_flat$id <- paste(all_under_table_flat$gene, all_under_table_flat$start,
                                 all_under_table_flat$dir, sep=" ")

all_under_table_flat_intragenic <- all_under_table_flat[all_under_table_flat$`feature type`=="intragenic",]
all_under_table_flat_intragenic <- merge(all_under_table_flat_intragenic, all_under_table_score_only, by='aln_id')
all_under_table_flat_intragenic<-arrange(all_under_table_flat_intragenic, venn_category) 
id_factor2 <- unique(all_under_table_flat_intragenic$id)
ht2v <- ggplot(all_under_table_flat_intragenic, 
               aes(factor(sample, c("MEP", "LSP","PMC7708049", "PMC9239280")),
                   factor(id, id_factor2),
                   fill= methyl)) + 
  geom_tile(color = "white",lwd = 1, linetype = 1) +
  scale_fill_gradient(low = "#FF9896FF",
                      high = "#AEC7E8FF") +
  scale_x_discrete(position = "top", labels = c("MEP", "LSP","set1", "set2")) +
  labs(fill = "methylation (%)") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin=unit(c(15,15,10,10), "pt"),
    axis.text.x=element_text(size = 10, vjust = 0.4, hjust = 0.4, color="black"),
    axis.text.y=element_text(size = 8, vjust = 0.4, hjust = 0.4, color="black"),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.spacing = unit(1, "mm"),
    
  ) +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

intragenic_data <- combined_data[combined_data$`feature type`=='intragenic',]
intragenic_data_MEP <- intragenic_data[intragenic_data$`methylation status`=="MEP-specific",]
gg12<-ggplot(combined_data[combined_data$`feature type`=='intragenic',], aes(x=LSP, y=`PMC9239280,%`, 
                                                                            colour = `methylation status`)) + 
  geom_point(alpha=0.7, size = 2.5) +
  labs(y = 'set2', x = 'LSP',
       color = "methylation status") +
  scale_colour_manual(values=status_pal_2) +
  theme_ipsum(base_family = 'Arial', base_size = 10) +
  guides(colour = guide_legend(nrow = 2)) +
  theme(
    axis.title.x = element_text(size = 12, vjust = 0.8, hjust = 0.5),
    axis.title.y = element_text(size = 12, vjust = 0.8, hjust = 0.5),
    plot.margin=unit(c(15,15,10,10), "pt"),
    axis.text.x=element_text(angle = 90, size = 10, vjust = 0.4, hjust = 0.4),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.spacing = unit(1, "mm")
  )
gg12

#intragenic figures

figure_intra <- ggarrange(ht2v, gg12, common.legend=FALSE,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1, vjust=1, 
                    legend="bottom")
ggsave(filename = "PMC9239280_vs_my_data_MEP_LSP_intragenic_fig_S3.1.png", 
       plot =figure_intra, path = results_dir,
       scale = 1, width = 180,
       height = 120, units = c("mm"),
       dpi = 600, bg = "white")
