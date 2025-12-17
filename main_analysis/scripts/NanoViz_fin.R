require(extrafont)
extrafont::loadfonts(device="all")

library(ggplot2)
library(dplyr)
library(tidyverse)
library(hrbrthemes)
library("ggpubr")
library("pheatmap")
library("janitor")

paper_pal <-paletteer::paletteer_d("ggsci::category20_d3")
working_dir = "/app/data" # the path to the working directory, you should put your own
results_dir = "/app/results"

#wd = "/home/anna/Documents/Methylome/results/experiments/nanopore_salmonella/sal_ent_barcodes_6mA/aligned/"

lsp_1<-"LSP-1_yibF_mtlA_3896107-3896582.bgz"
mep_2<-"MEP-2_yibF_mtlA_3896107-3896582.bgz"

region=c(3896107,3896582)
site_pos = c(3896381,3896270,3896227) #GATC site position

site="GATC"
gene="SL1344_3650_mtlA" 

methy_data_gene <- read.table(
  gzfile(paste(working_dir, lsp_1, sep="/")))
methy_data_gene <- methy_data_gene %>% row_to_names(row_number = 1)
#plus
gene_exact_LSP_plus <- methy_data_gene[methy_data_gene$ref_position>=region[1]&methy_data_gene$ref_position<=region[2]&
                                    methy_data_gene$ref_strand=="+",]
gene_exact_LSP_plus$mod_qual<-as.numeric(gene_exact_LSP_plus$mod_qual)
gene_exact_LSP_plus$base_qual<-as.numeric(gene_exact_LSP_plus$base_qual)



data_summary_gene_plus_LSP <- gene_exact_LSP_plus %>%
  group_by(read_id,ref_position) %>%
  summarize(qual = mean(mod_qual), .groups="keep")
data_summary_gene_plus_LSP_pw <- data_summary_gene_plus_LSP%>% 
  pivot_wider(names_from = ref_position, 
              values_from = qual)
data_summary_gene_plus_LSP_pw[is.na(data_summary_gene_plus_LSP_pw)] <--1

rownames(data_summary_gene_plus_LSP_pw)<-data_summary_gene_plus_LSP_pw$read_id

data_summary_gene_plus_LSP_gatc_sites <- data_summary_gene_plus_LSP[data_summary_gene_plus_LSP$ref_position %in%
                                                                      c(site_pos+1),]
#base quality vs meth score for A in GATC
#pos <- gene_exact_LSP_plus[gene_exact_LSP_plus$ref_position%in%data_summary_gene_plus_LSP_gatc_sites$ref_position,]
#ggplot(pos, aes(x=base_qual, y=mod_qual)) + 
#  geom_point()
##heatmap###

#preparing clustered heatmap
data_summary_gene_plus_LSP_gatc_sites_pv <- data_summary_gene_plus_LSP_gatc_sites %>%
  pivot_wider(names_from = ref_position, 
              values_from = qual)
data_summary_gene_plus_LSP_gatc_sites_pv <- select(data_summary_gene_plus_LSP_gatc_sites_pv, c("3896228", "3896271", "3896382"))
data_summary_gene_plus_LSP_gatc_sites_pv<-data_summary_gene_plus_LSP_gatc_sites_pv[complete.cases(data_summary_gene_plus_LSP_gatc_sites_pv), ]

matrix_plus_lsp <- data_summary_gene_plus_LSP_gatc_sites_pv[,-1]
rownames(matrix_plus_lsp) <- data_summary_gene_plus_LSP_gatc_sites_pv$read_id

#all qvalue heatmap
colors<-c("#1F77B4FF","#98DF8AFF","#DBDB8DFF","#FFBB78FF","#D62728FF")
rev(colors)
rf_2 <- colorRampPalette(rev(colors))
r_2 <- rf_2(10)
r_2
pheatmap(as.matrix(matrix_plus_lsp), cluster_cols = FALSE, color=r_2,
           show_rownames = FALSE, show_colnames=FALSE, legend=F, fontsize=12, 
         labels_col = c("GATC-1", "GATC-2", "GATC-3"),
           filename=paste(results_dir, "pheatmap_plus_lsp_qval_fig9_D.tiff", sep = "/"))

#minus
gene_exact_LSP_minus <- methy_data_gene[methy_data_gene$ref_position>=region[1]&methy_data_gene$ref_position<=region[2]&
                                         methy_data_gene$ref_strand=="-",]
gene_exact_LSP_minus$mod_qual<-as.numeric(gene_exact_LSP_minus$mod_qual)
data_summary_gene_minus_LSP <- gene_exact_LSP_minus %>%
  group_by(read_id,ref_position) %>%
  summarize(qual = mean(mod_qual), .groups="keep")
data_summary_gene_minus_LSP_pw <- data_summary_gene_minus_LSP%>% 
  pivot_wider(names_from = ref_position, 
              values_from = qual)
data_summary_gene_minus_LSP_pw[is.na(data_summary_gene_minus_LSP_pw)] <--1

rownames(data_summary_gene_minus_LSP_pw)<-data_summary_gene_minus_LSP_pw$read_id

data_summary_gene_minus_LSP_gatc_sites <- data_summary_gene_minus_LSP[data_summary_gene_minus_LSP$ref_position %in%
                                                                        c(site_pos+2),]
data_summary_gene_minus_LSP_gatc_sites_pv <- data_summary_gene_minus_LSP_gatc_sites %>%
  pivot_wider(names_from = ref_position, 
              values_from = qual)
data_summary_gene_minus_LSP_gatc_sites_pv <- select(data_summary_gene_minus_LSP_gatc_sites_pv, c("3896229", "3896272", "3896383"))
data_summary_gene_minus_LSP_gatc_sites_pv<-data_summary_gene_minus_LSP_gatc_sites_pv[complete.cases(data_summary_gene_minus_LSP_gatc_sites_pv), ]
matrix_minus_lsp <- data_summary_gene_minus_LSP_gatc_sites_pv[,-1]
rownames(matrix_minus_lsp) <- data_summary_gene_minus_LSP_gatc_sites_pv$read_id

#all qvalue heatmap

pheatmap(as.matrix(matrix_minus_lsp), cluster_cols = FALSE, color=r_2,
         show_rownames = FALSE, show_colnames=FALSE, legend=F, fontsize=12, 
         labels_col = c("GATC-1", "GATC-2", "GATC-3"),
         filename=paste(results_dir,"pheatmap_minus_lsp_qval_fig_9_E.tiff", sep = "/"))

###MEP###
methy_data_gene_MEP <- read.table(
  gzfile(paste(working_dir, mep_2, sep="/")))
methy_data_gene_MEP <- methy_data_gene_MEP %>% row_to_names(row_number = 1)
#plus
gene_exact_MEP_plus <- methy_data_gene_MEP[methy_data_gene_MEP$ref_position>=region[1]&
                                             methy_data_gene_MEP$ref_position<=region[2]&
                                         methy_data_gene_MEP$ref_strand=="+",]
gene_exact_MEP_plus$mod_qual<-as.numeric(gene_exact_MEP_plus$mod_qual)
data_summary_gene_plus_MEP <- gene_exact_MEP_plus %>%
  group_by(read_id,ref_position) %>%
  summarize(qual = mean(mod_qual), .groups="keep")
data_summary_gene_plus_MEP_pw <- data_summary_gene_plus_MEP%>% 
  pivot_wider(names_from = ref_position, 
              values_from = qual)
data_summary_gene_plus_MEP_pw[is.na(data_summary_gene_plus_MEP_pw)] <--1

rownames(data_summary_gene_plus_MEP_pw)<-data_summary_gene_plus_MEP_pw$read_id

data_summary_gene_plus_MEP_gatc_sites <- data_summary_gene_plus_MEP[data_summary_gene_plus_MEP$ref_position %in%
                                                                      c(site_pos+1),]

data_summary_gene_plus_MEP_gatc_sites_pv <- data_summary_gene_plus_MEP_gatc_sites %>%
  pivot_wider(names_from = ref_position, 
              values_from = qual)
data_summary_gene_plus_MEP_gatc_sites_pv <- select(data_summary_gene_plus_MEP_gatc_sites_pv, c("3896228", "3896271", "3896382"))
data_summary_gene_plus_MEP_gatc_sites_pv<-data_summary_gene_plus_MEP_gatc_sites_pv[complete.cases(data_summary_gene_plus_MEP_gatc_sites_pv), ]
matrix_plus_MEP <- data_summary_gene_plus_MEP_gatc_sites_pv[,-1]
rownames(matrix_plus_MEP) <- data_summary_gene_plus_MEP_gatc_sites_pv$read_id

pheatmap(as.matrix(matrix_plus_MEP), cluster_cols = FALSE, color=r_2,
         show_rownames = FALSE, show_colnames=FALSE, legend=F, fontsize=12, 
         labels_col = c("GATC-1", "GATC-2", "GATC-3"),
         filename=paste(results_dir,"pheatmap_plus_mep_qval_fig_9_B.tiff", sep = "/"))

#minus
gene_exact_MEP_minus <- methy_data_gene_MEP[methy_data_gene_MEP$ref_position>=region[1]&
                                              methy_data_gene_MEP$ref_position<=region[2]&
                                          methy_data_gene_MEP$ref_strand=="-",]
gene_exact_MEP_minus$mod_qual<-as.numeric(gene_exact_MEP_minus$mod_qual)
data_summary_gene_minus_MEP <- gene_exact_MEP_minus %>%
  group_by(read_id,ref_position) %>%
  summarize(qual = mean(mod_qual), .groups="keep")
data_summary_gene_minus_MEP_pw <- data_summary_gene_minus_MEP%>% 
  pivot_wider(names_from = ref_position, 
              values_from = qual)
data_summary_gene_minus_MEP_pw[is.na(data_summary_gene_minus_MEP_pw)] <--1

rownames(data_summary_gene_minus_MEP_pw)<-data_summary_gene_minus_MEP_pw$read_id
data_summary_gene_minus_MEP_gatc_sites <- data_summary_gene_minus_MEP[data_summary_gene_minus_MEP$ref_position %in%
                                                                        c(site_pos+2),]

data_summary_gene_minus_MEP_gatc_sites_pv <- data_summary_gene_minus_MEP_gatc_sites %>%
  pivot_wider(names_from = ref_position, 
              values_from = qual)
data_summary_gene_minus_MEP_gatc_sites_pv <- select(data_summary_gene_minus_MEP_gatc_sites_pv, c("3896229", "3896272", "3896383"))
data_summary_gene_minus_MEP_gatc_sites_pv<-data_summary_gene_minus_MEP_gatc_sites_pv[complete.cases(data_summary_gene_minus_MEP_gatc_sites_pv), ]
matrix_minus_MEP <- data_summary_gene_minus_MEP_gatc_sites_pv[,-1]
rownames(matrix_minus_MEP) <- data_summary_gene_minus_MEP_gatc_sites_pv$read_id

pheatmap(as.matrix(matrix_minus_MEP), cluster_cols = FALSE, color=r_2,
         show_rownames = FALSE, show_colnames=FALSE, legend=T, fontsize=12, 
         labels_col = c("GATC-1", "GATC-2", "GATC-3"),
         filename=paste(results_dir, "pheatmap_minus_mep_qval_fig_9_C.tiff", sep = "/"))

