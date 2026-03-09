library(extrafont)
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
colors = c("#1f77b4", "#ff7f0e","#2ca02c", "#d62728")
paper_pal <-paletteer::paletteer_d("ggsci::category20_d3")
status_pal <- paper_pal[1:4]
violin_color <- c("#C49C94FF","#F7B6D2FF","#C7C7C7FF")
chi_test_color <- c("#7F7F7FFF","#C7C7C7FF")
status_pal_2 <- c("#AEC7E8FF","#FFBB78FF","#98DF8AFF","#FF9896FF")
scale_fill_cvi_d = function(name) {
  ggplot2::scale_fill_manual(values = colors(name,
                                             type = "discrete"))
}
working_dir = "/app/data"
results_dir = "/app/results"
#folder = "/home/anna/old/mnt/Documents/Methylome/results/experiments/nanopore_salmonella/transcriptomic"
tpm <- read.csv(paste(working_dir,"salmonella_MEP_ESP_LSP_TPM.csv",sep="/"),sep=",")
info <-read.csv(paste(working_dir,"sample_info.txt", sep="/"),sep="\t")
rm_info <-read.csv(paste(working_dir,"rm_genes_info.txt", sep="/"),sep=",")
methyl <-read.csv(paste(working_dir,"methylation_level.csv", sep="/"),sep=",")
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

#p1 as heatmap
rm_annotated_group_info_order[rm_annotated_group_info_order$locus_tag%in%methyl$locus_tag,]$Name
rm_info[rm_info$locus_tag%in%methyl$locus_tag,]$Name
name_order = c("M.Sen371P, CAGAG","RM.Sen4693P, GATCAG",
               "M.Sen4729P, GAGN(6)RTAYG","YhdJ, ATGCAT","Dcm, CCWGG","Dam, GATC")
names_sites <-  c("M.Sen371P,\nCAGAG", "RM.Sen4693P,\nGATCAG","M.Sen4729P,\nGAGN(6)RTAYG","YhdJ,\nATGCAT","Dcm,\nCCWGG","Dam,\nGATC")

p1_tile <- ggplot(rm_annotated_group_info_order
                  [rm_annotated_group_info_order$group!="ESP"&
                      rm_annotated_group_info_order$locus_tag%in%methyl$locus_tag,], 
                aes(x=factor(group, levels = c("MEP", "LSP")), #factor(sample, levels = c("ESP", "LSP")
              y=factor(Name, levels = name_order), fill= TPM_mean)) + 
  geom_tile(color = "grey",
            lwd = 0.5,
            linetype = 1) +
  scale_fill_gradientn(colors = c("white","yellow", "red"), 
                       values = c(0,0.5,1))+
  scale_y_discrete(labels = names_sites) +
  guides(fill = guide_colourbar(title = "TPM", theme(legend.title.position = "bottom"), 
                                barwidth = 0.7)) +
  theme(axis.line = element_line(colour = "white", linewidth = .5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x=element_text(size = 10, vjust = 0.5, hjust = 0.5),
        axis.text.y=element_text(size = 9, vjust = 0.5, hjust = 0.5, color = "black"),
        axis.ticks = element_line(color="black", linewidth = 0.1, linetype = 1),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey90"),
        strip.text = element_text(size = 8),
        axis.ticks.length = unit(.25, "cm"),
        legend.position = 'right') 
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

deseq2_subset$LFC_mod <- with( 
  deseq2_subset, ifelse(deseq2_subset$padj<0.05, log2FoldChange, 0)) 


SLP1_0049 <- deseq2_subset[deseq2_subset$locus_tag=="SLP1_0049",]
SLP2_0027<- deseq2_subset[deseq2_subset$locus_tag=="SLP2_0027",]
SLP1_0049$sample="LSP"
SLP2_0027$sample="LSP"
deseq2_subset <- rbind(deseq2_subset,SLP1_0049,SLP2_0027)
deseq2_subset_mt <- deseq2_subset[(deseq2_subset$Mtase_activity=="yes") &
                                    deseq2_subset$sample!="ESP",]
print(deseq2_subset)
deseq2_subset_mt$comp <- "Log2FC"

p2_mt <- ggplot(deseq2_subset_mt[deseq2_subset_mt$sample=="LSP"&
                                   deseq2_subset_mt$locus_tag%in%methyl$locus_tag,], 
                aes(x=sample, #factor(sample, levels = c("ESP", "LSP")
                                         y=factor(Name, levels = name_order), 
                    fill= LFC_mod)) + 
  geom_tile(color = "grey",
            lwd = 0.5,
            linetype = 1) +
  scale_x_discrete(labels = c("LSP vs MEP"))+
  scale_fill_gradient(low = "royalblue", high = "grey")+
  guides(fill = guide_colourbar(title = "Log2FC",
                                theme(legend.title.position = "bottom"),
                                barwidth = 0.7)) +
  theme(axis.line = element_line(colour = "white", linewidth = .5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        #axis.text.x=element_blank(),
        axis.text.x=element_text(size = 10, vjust = 0.5, hjust = 0.5),
        axis.text.y=element_blank(), #text(size = 12, vjust = 0.5, hjust = 0.5),
        axis.ticks = element_line(color="black", linewidth = 0.1, linetype = 1),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey90"),
        strip.text = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        legend.position = 'right') 


#methylation level
methyl_annotated <-merge(methyl, rm_info, by = "locus_tag")
p3_methyl <- ggplot(methyl_annotated, 
                aes(x=factor(sample, levels = c("MEP", "LSP")),
                    y=factor(Name, levels = name_order), fill= X..methylated.sites)) + 
  geom_tile(color = "grey",
            lwd = 0.5,
            linetype = 1) +
  #geom_text(aes(label = round(X..methylated.sites,1)), color = "black", size = 3) +
  scale_fill_gradient(low = "#FFBB78FF", high = "#FF7F0EFF")+
  guides(fill = guide_colourbar(title = "% methylated\nsites",
                                theme(legend.title.position = "bottom"),
                                barwidth = 0.7)) +
  theme(axis.line = element_line(colour = "white", linewidth = .5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        #axis.text.x=element_blank(),
        axis.text.x=element_text(size = 10, vjust = 0.5, hjust = 0.5),
        axis.text.y=element_blank(), #text(size = 12, vjust = 0.5, hjust = 0.5),
        axis.ticks = element_line(color="black", linewidth = 0.1, linetype = 1),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey90"),
        strip.text = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        legend.position = 'right') 

tpm_lfc_methyl <- ggarrange(p1_tile, p2_mt, p3_methyl, ncol = 3, nrow = 1, 
          widths = c(1.2,0.6,0.7), align = "h", labels = c("A", "B", "C"))
tpm_lfc_methyl

p2_rm <- ggplot(deseq2_subset_mt[deseq2_subset_mt$sample=="LSP"&
                                   deseq2_subset_mt$Orphan=="no",], aes(y=sample, #factor(sample, levels = c("ESP", "LSP")
                                                                   x=factor(Name, levels = name_level), fill= LFC_mod)) + 
  geom_tile(color = "white",
            lwd = 1.0,
            linetype = 1) +
  geom_text(aes(label = round(LFC_mod,2)), color = "white", size = 6) +
  facet_grid(~comp)+
  coord_flip()+
  theme(axis.line = element_line(colour = "white", linewidth = .5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        #axis.text.x=element_blank(),
        axis.text.x=element_text(size = 12, vjust = 0.5, hjust = 0.5),
        axis.text.y=element_blank(),
        axis.ticks = element_line(color="black", linewidth = 0.1, linetype = 1),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey90"),
        strip.text = element_text(size = 16),
        #legend.position = "none",
        axis.ticks.length = unit(.25, "cm")) 

#mean methylation level
gatc_methyl <- read.csv(paste(working_dir,"GATC_data_all_flat_table_weight_av.csv", sep = "/"),sep="\t")

gatc_flat_MEP<-data.frame(methyl = gatc_methyl$MEP_wt_mean, sample = c("MEP"))
gatc_flat_LSP<-data.frame(methyl = gatc_methyl$LSP_wt_mean, sample = c("LSP"))
gatc_flat <- rbind(gatc_flat_MEP,gatc_flat_LSP)
gatc_flat$site <- "GATC"
ccwgg_methyl <- read.csv(paste(working_dir,"CCWGG_data_all_flat_table_weight_av.csv", sep = "/"),sep="\t")
ccwgg_flat_MEP<-data.frame(methyl = ccwgg_methyl$MEP_wt_mean, sample = c("MEP"))
ccwgg_flat_LSP<-data.frame(methyl = ccwgg_methyl$LSP_wt_mean, sample = c("LSP"))
ccwgg_flat <- rbind(ccwgg_flat_MEP,ccwgg_flat_LSP)
ccwgg_flat$site <- "CCWGG"
atgcat_methyl <- read.csv(paste(working_dir, "ATGCAT_data_all_flat_table_weight_av.csv", sep = "/"),sep="\t")
atgcat_flat_MEP<-data.frame(methyl = atgcat_methyl$MEP_wt_mean, sample = c("MEP"))
atgcat_flat_LSP<-data.frame(methyl = atgcat_methyl$LSP_wt_mean, sample = c("LSP"))
atgcat_flat <- rbind(atgcat_flat_MEP,atgcat_flat_LSP)
atgcat_flat$site <- "ATGCAT"
gagn6_methyl <- read.csv(paste(working_dir, "GAGN6RTAYG_data_all_flat_table_weight_av.csv", sep = "/"),sep="\t")
gagn6_flat_MEP<-data.frame(methyl = gagn6_methyl$MEP_wt_mean, sample = c("MEP"))
gagn6_flat_LSP<-data.frame(methyl = gagn6_methyl$LSP_wt_mean, sample = c("LSP"))
gagn6_flat <- rbind(gagn6_flat_MEP,gagn6_flat_LSP)
gagn6_flat$site <- "GAGN6RTAYG"
cagag_methyl <- read.csv(paste(working_dir,"CAGAG_data_all_flat_table_weight_av.csv", sep = "/"),sep="\t")
cagag_flat_MEP<-data.frame(methyl = cagag_methyl$MEP_wt_mean, sample = c("MEP"))
cagag_flat_LSP<-data.frame(methyl = cagag_methyl$LSP_wt_mean, sample = c("LSP"))
cagag_flat <- rbind(cagag_flat_MEP,cagag_flat_LSP)
cagag_flat$site <- "CAGAG"
gatcag_methyl <- read.csv(paste(working_dir,"GATCAG_data_all_flat_table_weight_av.csv", sep = "/"),sep="\t")
gatcag_flat_MEP<-data.frame(methyl = gatcag_methyl$MEP_wt_mean, sample = c("MEP"))
gatcag_flat_LSP<-data.frame(methyl = gatcag_methyl$LSP_wt_mean, sample = c("LSP"))
gatcag_flat <- rbind(gatcag_flat_MEP,gatcag_flat_LSP)
gatcag_flat$site <- "GATCAG"
table_flat <- rbind(gatc_flat, ccwgg_flat, atgcat_flat, gagn6_flat, cagag_flat, gatcag_flat)
p <- ggplot(table_flat, aes(x=factor(sample, level = c("MEP", "LSP")), y=methyl, fill = sample)) + 
  geom_violin() +
  geom_boxplot(width=0.1) +
  #coord_flip() +
  facet_grid(.~factor(site, level = c("GATC", "CCWGG", "ATGCAT", "GAGN6RTAYG", "GATCAG", "CAGAG"))) +
  scale_fill_manual(values=c("#BCBD22FF","#17BECFFF")) +
  labs(y = "methylation (%)") +
  #stat_summary(fun.data=mean_sdl, geom="point", color="grey30") +
  theme(axis.line = element_line(colour = "white", linewidth = .5),
        axis.title.y = element_text(size = 10, vjust = 0.5, hjust = 0.5),
        axis.title.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"), 
        #                                size=1.5, linetype="solid"),
        strip.text.x = element_text(size=9, color="black"),
        axis.text.x=element_text(size = 10, vjust = 0.5, hjust = 0.5),
        axis.text.y=element_text(size = 10, vjust = 0.5, hjust = 0.5),
        axis.ticks = element_line(color="black", linewidth = 0.1, linetype = 1),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey90"),
        #strip.text = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        legend.position = 'none') 

tpm_lfc_methyl_mean <-ggarrange(ggarrange(p1_tile, p2_mt, p3_methyl,ncol = 3,
                    labels = c("A", "B", "C"),widths = c(1.2,0.5,0.8)),                                                 
          p, nrow=2, labels = c("","D")) 
tpm_lfc_methyl_mean
ggsave(filename = "Fig_5_TPM_LFC_methyl_mean.png", plot =tpm_lfc_methyl_mean, 
       path = results_dir,
       scale = 1, width = 190,
       height = 120, units = c("mm"),
       dpi = 300, bg = "white")

ggsave(filename = "Fig_5_TPM_LFC_methyl_mean.tiff", plot =tpm_lfc_methyl_mean, 
       path = results_dir,
       scale = 1, width = 190,
       height = 120, units = c("mm"),
       dpi = 300, bg = "white")
