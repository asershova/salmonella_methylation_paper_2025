#!/bin/bash
#run all scripts to produce all results
if [ scripts/bed_data_analysis_sal_GATC_28_08.R -nt results/GATC_jitter.tiff \
  -o dockerfiles/Dockerfile.bed_data_analysis_sal_GATC_28_08 -nt results/GATC_jitter.tiff \
   ]; then
  echo "======= Rebuilding bed_data_analysis_sal_GATC_28_08 ======="
  docker build -t bed_data_analysis_sal_gatc:latest -f dockerfiles/Dockerfile.bed_data_analysis_sal_GATC_28_08 .

  docker run -v $(pwd)/results:/app/results:rw bed_data_analysis_sal_gatc:latest
fi

if [ scripts/bed_data_analysis_sal_CCWGG_7_09_2.R -nt results/CCWGG_jitter.tiff \
  -o dockerfiles/Dockerfile.bed_data_analysis_sal_CCWGG_7_09_2 -nt results/CCWGG_jitter.tiff \
   ]; then
  echo "======= Rebuilding bed_data_analysis_sal_CCWGG_7_09_2 ======="
  docker build -t bed_data_analysis_sal_ccwgg:latest -f dockerfiles/Dockerfile.bed_data_analysis_sal_CCWGG_7_09_2 .

  docker run -v $(pwd)/results:/app/results:rw bed_data_analysis_sal_ccwgg:latest
fi

if [ scripts/bed_data_analysis_sal_ATGCAT_7_09_2.R -nt results/ATGCAT_jitter.tiff \
  -o dockerfiles/Dockerfile.bed_data_analysis_sal_ATGCAT_7_09_2 -nt results/ATGCAT_jitter.tiff \
   ]; then
  echo "======= Rebuilding bed_data_analysis_sal_ATGCAT_7_09_2 ======="
  docker build -t bed_data_analysis_sal_atgcat:latest -f dockerfiles/Dockerfile.bed_data_analysis_sal_ATGCAT_7_09_2 .

  docker run -v $(pwd)/results:/app/results:rw bed_data_analysis_sal_atgcat:latest
fi

if [ gene_expression_vis_salmonella_fin -nt results/Fig1.RM_expression.tiff \
  -o dockerfiles/Dockerfile.gene_expression_vis_salmonella_fin -nt results/Fig1.RM_expression.tiff \
   ]; then
  echo "======= Rebuilding gene_expression_vis_salmonella_fin ======="
  docker build -t gene_expression_vis_salmonella_fin:latest -f dockerfiles/Dockerfile.gene_expression_vis_salmonella_fin .

  docker run -v $(pwd)/results:/app/results:rw gene_expression_vis_salmonella_fin:latest
fi
if [ DEM_transcriptomic_annotation_GATC_fin -nt results/GATC_transcriptomic_figure_5.tiff \
  -o dockerfiles/Dockerfile.DEM_transcriptomic_annotation_GATC_fin -nt results/GATC_transcriptomic_figure_5.tiff \
   ]; then
  echo "======= Rebuilding DEM_transcriptomic_annotation_GATC_fin ======="
  docker build -t dem_transcriptomic_annotation_gatc_fin:latest -f dockerfiles/Dockerfile.DEM_transcriptomic_annotation_GATC_fin .

  docker run -v $(pwd)/results:/app/results:rw dem_transcriptomic_annotation_gatc_fin:latest
fi

if [ DEM_transcriptomic_annotation_CCWGG_fin -nt results/CCWGG_transcriptomic_figure_6.tiff \
  -o dockerfiles/Dockerfile.DEM_transcriptomic_annotation_CCWGG_fin -nt results/CCWGG_transcriptomic_figure_6.tiff \
   ]; then
  echo "======= Rebuilding DEM_transcriptomic_annotation_CCWGG_fin ======="
  docker build -t dem_transcriptomic_annotation_ccwgg_fin:latest -f dockerfiles/Dockerfile.DEM_transcriptomic_annotation_CCWGG_fin .

  docker run -v $(pwd)/results:/app/results:rw dem_transcriptomic_annotation_ccwgg_fin:latest
fi
if [ DEM_transcriptomic_annotation_ATGCAT_fin -nt results/ATGCAT_transcriptomic_figure_7.tiff \
  -o dockerfiles/Dockerfile.DEM_transcriptomic_annotation_ATGCAT_fin -nt results/ATGCAT_transcriptomic_figure_7.tiff \
   ]; then
  echo "======= Rebuilding DEM_transcriptomic_annotation_ATGCAT_fin ======="
  docker build -t dem_transcriptomic_annotation_atgcat_fin:latest -f dockerfiles/Dockerfile.DEM_transcriptomic_annotation_ATGCAT_fin .

  docker run -v $(pwd)/results:/app/results:rw dem_transcriptomic_annotation_atgcat_fin:latest
fi
if [ GATC_comparison_fin_intergenic -nt results/all_vs_all_fig_v_figure_8_AB.tiff \
  -o dockerfiles/Dockerfile.GATC_comparison_fin_intergenic -nt results/ATGCAT_transcriptomic_figure_7.tiff \
   ]; then
  echo "======= Rebuilding intergenic GATC literature data comparisons ======="
  docker build -t gatc_comparison_fin_intergenic:latest -f dockerfiles/Dockerfile.GATC_comparison_fin_intergenic .

  docker run -v $(pwd)/results:/app/results:rw gatc_comparison_fin_intergenic:latest
fi
if [ NanoViz_fin -nt results/pheatmap_plus_lsp_qval_fig9_D.tiff \
  -o dockerfiles/Dockerfile.NanoViz_fin -nt results/pheatmap_plus_lsp_qval_fig9_D.tiff \
   ]; then
  echo "======= Rebuilding heterogeneity example ======="
  docker build -t nanoviz_fin:latest -f dockerfiles/Dockerfile.NanoViz_fin .

  docker run -v $(pwd)/results:/app/results:rw nanoviz_fin:latest
fi
