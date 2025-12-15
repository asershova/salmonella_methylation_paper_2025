#options(warn=2)
require(extrafont)
extrafont::loadfonts(device="all")
library(RColorBrewer)
library("ggsci")
#library(ggtreeExtra)
library(ggstar)
library(ggplot2)
#library(ggtree)
#library(treeio)
library(ggnewscale)
library(data.table)
library(lattice)
library(reshape2)
library(dplyr)
library(forcats)
library(ggsci)
library("scales")
library(ggpubr)

working_dir = "/app/data" # the path to the working directory, you should put your own
results_dir = "/app/results"

tpm <- read.csv(paste(working_dir,"salmonella_MEP_ESP_LSP_TPM.csv",sep="/"),sep=",")
info <-read.csv(paste(working_dir,"sample_info.txt", sep="/"),sep="\t")
rm_info <-read.csv(paste(working_dir,"rm_genes_info.txt", sep="/"),sep=",")
colnames(info)[1]<-"sample"
colnames(tpm)[1]<-"locus_tag"

tpm_rm <- subset(tpm, locus_tag %in% rm_info$locus_tag)

rm_reshape <- melt(tpm_rm, id.vars=c("locus_tag"))
colnames(rm_reshape) <-c("locus_tag", "sample", "TPM")
rm_annotated<- merge(rm_reshape, info, by = "sample")

rm_annotated_group <- rm_annotated %>%
  group_by(locus_tag, group) %>%
  summarize(TPM_mean = mean(TPM), .groups="keep")

rm_annotated_group_info <- merge(rm_annotated_group, rm_info, 
                                 by = "locus_tag")

#  guides(fill = guide_legend(title="genes"))

rm_annotated_group_info_order <- rm_annotated_group_info %>%
  mutate(REBASE_Name = fct_reorder(REBASE_Name, System.Name))

palette<-c("#393B79FF", "#637939FF", "#637939FF", "#8C6D31FF", "#843C39FF", "#7B4173FF", "#5254A3FF", "#8CA252FF",
           "#BD9E39FF", "#AD494AFF","#AD494AFF","#A55194FF", "#6B6ECFFF","#6B6ECFFF","#6B6ECFFF" )

p1_all <- ggplot(rm_annotated_group_info_order, 
                aes(x = factor(REBASE_Name), y = TPM_mean, fill=REBASE_Name))+
  geom_bar(position="dodge", stat="identity")+
  facet_grid(~factor(group, levels=c('MEP','ESP','LSP'))) +
  labs(y = "Transcripts per million")+
  #ylim(0,22) +
  scale_fill_manual(values = palette) +
  coord_flip()+
  theme(axis.line = element_line(colour = "white", linewidth = .5),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.x=element_text(size = 12, vjust = 0.5, hjust = 0.5),
        axis.text.y=element_text(size = 12, vjust = 0.5, hjust = 0.5, color=palette, 
                                 margin = margin(t = .4, unit = "cm")),
        axis.ticks = element_line(color="black", linewidth = 0.2, linetype = 1),
        legend.position = "none",
        panel.background = element_rect(fill = NA),
        strip.text = element_text(size = 12),
        panel.grid.major = element_line(colour = "grey90"),
        axis.ticks.length = unit(.05, "cm")) 

#  scale_fill_manual(values = c("#1f77b4", "#ff7f0e","#2ca02c", "#d62728"))
#DESeq2
esp_WT <- read.csv(paste(working_dir,"ESP_MEP.clipped.deseq_results.05.csv",sep="/"),sep=",")
lsp_WT <- read.csv(paste(working_dir,"LSP_MEP.clipped.deseq_results.05.csv",sep="/"),sep=",")
esp_WT$sample <- "ESP"
lsp_WT$sample <- "LSP"
deseq2_table <- rbind(esp_WT,lsp_WT)
colnames(deseq2_table)[1] <- "locus_tag"
deseq2_subset <- subset(deseq2_table, locus_tag %in% rm_info$locus_tag)

deseq2_subset$padj_new <- round(deseq2_subset$padj,2)
deseq2_subset <-merge(deseq2_subset, rm_info, by = "locus_tag")
sample_levels <-c("ESP", "LSP")
name_level_all <- c("M.Sen474DamP", "M.Sen474DcmP", "V.Sen474DcmP",
                "M.Sen474ORF1050P", "M.Sen474ORF2068P", "M.Sen474ORF221P",
                "M.Sen474ORF2736P", "M.Sen474ORF2855P", "M.Sen474ORF3548P",
                "M.Sen474ORF371P","Sen474ORF371P", "RM.Sen474ORF4693P","M.Sen474ORF4729P", 
                "S.Sen474ORF4729P","Sen474ORF4729P")
deseq2_subset$LFC_mod <- with( 
  deseq2_subset, ifelse(deseq2_subset$padj<0.01, log2FoldChange, 0)) 


SLP1_0049 <- deseq2_subset[deseq2_subset$locus_tag=="SLP1_0049",]
SLP2_0027<- deseq2_subset[deseq2_subset$locus_tag=="SLP2_0027",]
SLP1_0049$sample="LSP"
SLP2_0027$sample="LSP"

deseq2_subset <- rbind(deseq2_subset,SLP1_0049,SLP2_0027)
p2_all <- ggplot(deseq2_subset, aes(y=factor(sample, levels = sample_levels), 
                                                   x=factor(REBASE_Name, levels = name_level_all), fill= LFC_mod)) + 
  geom_tile(color = "white",
            lwd = 1.0,
            linetype = 1) +
  geom_text(aes(label = round(LFC_mod,1)), color = "white", size = 3) +
  coord_fixed()+
  labs(fill = "Log2FC") +
  theme(axis.line = element_line(colour = "white", linewidth = .5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        #axis.text.x=element_blank(),
        axis.text.x=element_text(size = 12,angle = 90, vjust = 0.5, hjust = 0.5, color=palette),
        axis.text.y=element_text(size = 12, vjust = 0.5, hjust = 0.5, color = "black"),
        axis.ticks = element_line(color="black", linewidth = 0.1, linetype = 1),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey90"),
        strip.text = element_text(size = 12),
        axis.ticks.length = unit(.25, "cm")) 

tpm_lfc_graph <- ggarrange(p1_all, p2_all, nrow=2, labels = c("A", "B"), heights = c(1, 0.9))

ggsave(filename = "Fig1.RM_expression.tiff", plot =tpm_lfc_graph, path = results_dir,
       scale = 1, width = 180,
       height = 170, units = c("mm"),
       dpi = 300, bg = "white")
ggsave(filename = "Fig1.RM_expression.png", plot =tpm_lfc_graph, path = results_dir,
       scale = 1, width = 180,
       height = 170, units = c("mm"),
       dpi = 300, bg = "white")
