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
gatc_methyl <- read.csv(paste(working_dir,"CCWGG_data_all_flat_table_weight_av.csv", sep="/"), sep="", header = TRUE)
lsp_mep <- read.csv(paste(working_dir,"LSP_MEP.clipped.deseq_results.05.csv", sep = "/"), sep=",", header = FALSE,skip=1,col.names=c("locus_tag", "baseMean","log2FoldChange", "lfcSE","stat","pvalue", "padj"))

site = "CCWGG"

paper_pal <-paletteer::paletteer_d("ggsci::category20_d3")
paper_pal
#status_pal <- paper_pal[1:4]
status_pal_2 <- c("#AEC7E8FF","#FFBB78FF","#98DF8AFF","#FF9896FF")
status_pal_2_ccwgg <- c("#AEC7E8FF","#FFBB78FF","#FF9896FF")
###DEG:abs(Log2FC)>=1&padj<0.05
lsp_mep <- lsp_mep %>% mutate(
  DEG=case_when(
    abs(log2FoldChange) >=1 &
      padj <0.05 ~"DEG",
    TRUE ~ "non-DEG"
  )
)
###add log2FC methylation ###
gatc_methyl<- gatc_methyl %>%
  mutate(log2FCMethyl= case_when(
    MEP_wt_mean<35&LSP_wt_mean<35 ~ 0,
    TRUE~log2(gatc_methyl$LSP_wt_mean/(gatc_methyl$MEP_wt_mean+0.01))
  )
  )

#get affected genes from different intergenic regions:
#-+ - take both
#+- - control
#++ - take the second gene
#-- - take the first gene
all_intergenic <- gatc_methyl[gatc_methyl$feature=="intergenic",]

plus_minus_1_all <- all_intergenic[all_intergenic$f_dir=="-+",]
plus_minus_1_all$affected_lt<-plus_minus_1_all$locus_tag
plus_minus_1_all$rel_genomic_pos<-plus_minus_1_all$f_start-plus_minus_1_all$site_start+1

plus_minus_2_all <- all_intergenic[all_intergenic$f_dir=="-+",]
plus_minus_2_all$rel_genomic_pos<-plus_minus_2_all$site_start-plus_minus_2_all$f_stop+1
plus_minus_2_all$affected_lt<-plus_minus_2_all$locus_tag1
plus_2_all <- all_intergenic[all_intergenic$f_dir=="++",]
plus_2_all$rel_genomic_pos<-plus_2_all$site_start-plus_2_all$f_stop+1
plus_2_all$affected_lt<-plus_2_all$locus_tag1
minus_1_all <- all_intergenic[all_intergenic$f_dir=="--",]
minus_1_all$affected_lt<-minus_1_all$locus_tag
minus_1_all$rel_genomic_pos<-minus_1_all$f_start-minus_1_all$site_start+1

minus_plus_all <- all_intergenic[all_intergenic$f_dir=="+-",]
minus_plus_all$affected_lt<-"None"
minus_plus_all$rel_genomic_pos<-minus_plus_all$f_start-minus_plus_all$site_start+1

all_intergenic_combined <- rbind(plus_minus_1_all, plus_minus_2_all,plus_2_all, minus_1_all, minus_plus_all)

gene_body_all <- gatc_methyl[gatc_methyl$feature=="intragenic",]
gene_body_all$affected_lt<- gene_body_all$locus_tag
gene_body_all<- gene_body_all %>%
  mutate(rel_genomic_pos= case_when(
    f_dir=="+" ~ site_start-f_start+1,
    f_dir=="-" ~ f_stop-site_stop+1
  )
  )
gene_body_nu<-gene_body_all[!(gene_body_all$affected_lt%in%all_intergenic_combined$affected_lt),]
gene_body_all[gene_body_all$methyl_status=="LSP-specific",]$affected_lt    
gatc_gene_data_all <- rbind(all_intergenic_combined, gene_body_all) 
gatc_gene_data_nu <- rbind(all_intergenic_combined, gene_body_nu)#gene body only genes not affected by intergenic methyllation

all_pos_transcript_total <- merge(gatc_gene_data_all, lsp_mep, by.x = "affected_lt", by.y="locus_tag")
write.csv(all_pos_transcript_total, file = paste(results_dir, "Suppl_Table_CCWGG_all_methyl_gene_expression.txt", sep=""),
          sep= ",", row.names = FALSE)

all_pos_transcript_total_pv <- all_pos_transcript_total %>%
  group_by(affected_lt, methyl_status,feature) %>%
  summarize(n=n(), .groups="keep")

all_pos_transcript <- merge(gatc_gene_data_nu, lsp_mep, by.x = "affected_lt", by.y="locus_tag")
um_gatc_list <-unique(gatc_gene_data_nu[gatc_gene_data_nu$methyl_status=="UM",]$affected_lt)

write.csv(unique(gatc_methyl[gatc_methyl$locus_tag%in%um_gatc_list,]$id2), file = paste(results_dir, "CCWGG_UM_upstream.txt", sep=""),
          sep= "\n", row.names = FALSE)
#plot with log2FCexpression/logFCmethyl

#figA
# New facet label names for supp variable
feature.labs <- c("Gene", "Upstream")
names(feature.labs) <- c("intragenic", "intergenic")

p_log2_FC<-ggplot(all_pos_transcript%>%
                    arrange(methyl_status),
                  aes(y=log2FCMethyl, x=log2FoldChange, color=methyl_status)) + 
  labs(x = "Log2FC") +
  geom_point(alpha=0.9) +
  facet_wrap(~feature, labeller = labeller(feature = feature.labs)) +
  scale_color_manual(values=status_pal_2_ccwgg) +
#  ylim(-2,4) +
  geom_vline(xintercept = -1, linetype="dashed") +
  geom_vline(xintercept = 1, linetype="dashed")+
  theme_ipsum(base_family = 'Arial', base_size = 12) +
  theme(
    axis.title.y = element_text(size = 12, vjust = 0.8, hjust = 0.5),
    axis.title.x = element_text(size = 12, vjust = 0.8, hjust = 0.5),
    axis.text.x=element_text(size = 12, vjust = 0.4, hjust = 0.4),
    axis.text.y=element_text(size = 12, vjust = 0.4, hjust = 0.4),
    plot.margin=unit(c(15,0,10,10), "pt"),
    legend.position="none",
    legend.spacing = unit(1, "mm")
    
  )
p_log2_FC


##methylation thresholds FC
inter_05_methyl_const<-all_pos_transcript[(all_pos_transcript$log2FCMethyl>=0.5)&
                     all_pos_transcript$feature=="intergenic"&
                       all_pos_transcript$methyl_status=="Constitutive",]
inter_1_methyl<-all_pos_transcript[(all_pos_transcript$log2FCMethyl>=1)&
                                            all_pos_transcript$feature=="intergenic",]
#LSP-specific genes vs remaining
inter_05_methyl_lsp<-all_pos_transcript[all_pos_transcript$feature=="intergenic"&
                                            (all_pos_transcript$methyl_status=="LSP-specific"|
                                               all_pos_transcript$methyl_status=="MEP-specific"),]
gb_05_methyl_lsp<-all_pos_transcript[all_pos_transcript$feature=="intragenic"&
                                          (all_pos_transcript$methyl_status=="LSP-specific"|
                                             all_pos_transcript$methyl_status=="MEP-specific"),]
inter_05_methyl_const<-inter_05_methyl_const[!(inter_05_methyl_const$affected_lt %in%inter_05_methyl_lsp$affected_lt),]

inter_05_methyl_lsp_genes_total_gatc <- all_pos_transcript_total_pv[all_pos_transcript_total_pv$affected_lt%in%inter_05_methyl_lsp$affected_lt,]
inter_05_methyl_lsp_genes_total_gatc_wide <- inter_05_methyl_lsp_genes_total_gatc %>%
  pivot_wider(names_from = methyl_status, 
              values_from = n)
inter_05_methyl_lsp_genes_total_gatc_wide_inter <- inter_05_methyl_lsp_genes_total_gatc_wide[
  inter_05_methyl_lsp_genes_total_gatc_wide$feature=="intergenic",]  #change NA to 0 and save the table as suppls
write.csv(unique(gatc_methyl[gatc_methyl$locus_tag%in%
            inter_05_methyl_lsp$affected_lt,]$id2), file = paste(results_dir, "CCWGG_DM_upstream_id.txt", sep=""),
          sep= "\n", row.names = FALSE)
gb_05_methyl_lsp[gb_05_methyl_lsp$log2FCMethyl>10,]
lsp_mep<- lsp_mep %>%mutate(
  DM_sites = case_when(
    locus_tag%in%inter_05_methyl_lsp$affected_lt~"DM",
    TRUE~"other")
)
lsp_mep<- lsp_mep %>%mutate(
  DM_sites_i_gb = case_when(
    locus_tag%in%gb_05_methyl_lsp$affected_lt~"DM in gene body",
    locus_tag%in%inter_05_methyl_lsp$affected_lt~"DM in upstream",
    TRUE~"other")
)
lsp_mep%>%
  group_by(DM_sites)%>%
    summarise(n=n())
lsp_mep%>%
  group_by(DM_sites_i_gb)%>%
  summarise(n=n())

my_comparisons_dm <- list(c("DM", "other"))
my_comparisons_igb <- list(c("DM in gene body", "other"),
                          c("DM in upstream", "other"),
                          c("DM in gene body", "DM in upstream"))
#figD
wilcox.test(lsp_mep[lsp_mep$DM_sites_i_gb=="DM in gene body",]$log2FoldChange,
            lsp_mep[lsp_mep$DM_sites_i_gb=="other",]$log2FoldChange)
violins_dm<-ggplot(lsp_mep, aes(y=log2FoldChange, x=factor(DM_sites_i_gb, 
                                                        levels=c("DM in upstream", "DM in gene body", "other")),
                             fill=DM_sites_i_gb)) + 
  geom_violin() +
  stat_compare_means(comparisons = my_comparisons_igb)+
  geom_boxplot(width=0.1, color="black", alpha=0.2) +
  scale_x_discrete(labels= c("Upstream\nn=268", "Gene\nn=1560", "Other\nn=3204"))+
  scale_fill_manual(values = c("#C49C94FF","#F7B6D2FF","#C7C7C7FF"))+  
  ylim(NA, 22) +
  theme_ipsum(base_family = 'Arial', base_size = 12) +
  theme(
    axis.title.y = element_text(size = 12, vjust = 0.8, hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.x=element_text(size = 12, vjust = 0.4, hjust = 0.4),
    axis.text.y=element_text(size = 12, vjust = 0.4, hjust = 0.4),
    plot.margin=unit(c(25,0,10,10), "pt"),
    legend.position="none")

violins_dm
#chi-square test
lsp_mep_igb_pv<-lsp_mep%>%
  group_by(DM_sites_i_gb, DEG)%>%
  summarise(n=n(), .groups="keep")

chi_inter_wide<-lsp_mep_igb_pv%>% 
  pivot_wider(names_from = DEG, 
              values_from = n)
chi_inter_wide
m_i_o = matrix(c(169,99,2025,1179), ncol=2)
m_i_gb = matrix(c(169,99,959,601), ncol=2)
m_gb_o=matrix(c(959,601,2025,1179),ncol=2)
p_chi_mio<-chisq.test(x=m_i_o)$p.value
p_chi_migb<-chisq.test(x=m_i_gb)$p.value
p_chi_gbo<-chisq.test(x=m_gb_o)$p.value

#figC
stacked_bar<-ggplot(lsp_mep_igb_pv, 
                    aes(fill=DEG, y=n, x=factor(DM_sites_i_gb, 
                                                       levels=c("DM in upstream", "DM in gene body", "other")))) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c("#7F7F7FFF","#C7C7C7FF"))+
  labs(x = "group", y = "genes (%)", fill = "") +
  scale_x_discrete(labels= c("Upstream\nn=268", "Gene\nn=1560","Other\nn=3204"))+
  ggsignif::geom_signif(annotations = c(0.26,1,0.67), comparisons=my_comparisons_igb, y_position=c(-95.15, -95.0, -95.25),
                        map_signif_level = TRUE,tip_length=0.00002, size = 0.3)+
  coord_cartesian(ylim = c(0, 1.4)) +
  scale_y_continuous(labels = scales::label_percent(suffix="")) +
  theme_ipsum(base_family = 'Arial', base_size = 12) +
  theme(
    axis.title.y = element_text(size = 12, vjust = 0.8, hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.x=element_text(size = 12, angle=45, vjust = 0.4, hjust = 0.4),
    axis.text.y=element_text(size = 12, vjust = 0.4, hjust = 0.4),
    plot.margin=unit(c(30,0,10,10), "pt"),
    legend.spacing = unit(1, "mm"),
    legend.text = element_text(size = 10)
    
  )
stacked_bar
#fig B
p_genomic_pos<-ggplot(all_pos_transcript%>%
                        arrange(methyl_status), #[all_pos_transcript$methyl_status!="Constitutive",]
                      aes(y=log2FCMethyl, x=rel_genomic_pos, color=methyl_status)) + 
  geom_point(alpha=0.7) +
  xlim(-500,+500) +
  scale_color_manual(values=status_pal_2_ccwgg) + #[2:4]
  theme_ipsum(base_family = 'Arial', base_size = 12) +
  geom_vline(xintercept = -250, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed")+
  labs(color = "methylation status", x = "Relative position") +
  theme(
    axis.title.y = element_text(size = 12, vjust = 0.8, hjust = 0.5),
    axis.title.x = element_text(size = 12, vjust = 0.8, hjust = 0.5),
    axis.text.x=element_text(size = 12, angle = 45, vjust = 0.4, hjust = 0.4),
    axis.text.y=element_text(size = 12, vjust = 0.4, hjust = 0.4),
    plot.margin=unit(c(15,0,10,10), "pt"),
    legend.spacing = unit(1, "mm")
  )
#D
p_genomic_pos
figure <- ggarrange(p_log2_FC,p_genomic_pos, stacked_bar, violins_dm,
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2, vjust=1)

ggsave(filename = paste(site,"_transcriptomic_figure_6.tiff", sep=""), plot =figure, path = results_dir,
       scale = 1, width = 180,
       height = 156, units = c("mm"),
       dpi = 300, bg = "white")
ggsave(filename = paste(site,"_transcriptomic_figure_6.png", sep=""), plot =figure, path = results_dir,
       scale = 1, width = 180,
       height = 156, units = c("mm"),
       dpi = 300, bg = "white")
