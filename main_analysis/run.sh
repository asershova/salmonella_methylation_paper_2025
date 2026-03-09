#!/bin/bash
#run all scripts to produce all results
if [ scripts/bed_data_analysis_GAATTC_ROC_salmonella.R -nt results/GAATTC_fig_1.tiff \
  -o dockerfiles/Dockerfile.bed_data_analysis_GAATTC_ROC_salmonella -nt results/GAATTC_fig_1.tiff \
   ]; then
  echo "======= Rebuilding bed_data_analysis_GAATTC_ROC_salmonella ======="
  docker build -t bed_data_analysis_gaattc_roc_salmonella:latest -f dockerfiles/Dockerfile.bed_data_analysis_GAATTC_ROC_salmonella .
  docker run -v $(pwd)/results:/app/results:rw bed_data_analysis_gaattc_roc_salmonella:latest
fi

if [ scripts/bed_data_analysis_sal_GATC.R -nt results/GATC_fig_2.tiff \
  -o dockerfiles/Dockerfile.bed_data_analysis_sal_GATC -nt results/GATC_fig_2.tiff \
   ]; then
  echo "======= Rebuilding bed_data_analysis_sal_GATC ======="
  docker build -t bed_data_analysis_sal_gatc:latest -f dockerfiles/Dockerfile.bed_data_analysis_sal_GATC .

  docker run -v $(pwd)/results:/app/results:rw bed_data_analysis_sal_gatc:latest
fi

if [ scripts/bed_data_analysis_sal_CCWGG.R -nt results/CCWGG_fig_4.tiff \
  -o dockerfiles/Dockerfile.bed_data_analysis_sal_CCWGG -nt results/CCWGG_fig_4.tiff \
   ]; then
  echo "======= Rebuilding bed_data_analysis_sal_CCWGG ======="
  docker build -t bed_data_analysis_sal_ccwgg:latest -f dockerfiles/Dockerfile.bed_data_analysis_sal_CCWGG .

  docker run -v $(pwd)/results:/app/results:rw bed_data_analysis_sal_ccwgg:latest
fi

if [ scripts/bed_data_analysis_sal_ATGCAT.R -nt results/ATGCAT_fig_3.tiff \
  -o dockerfiles/Dockerfile.bed_data_analysis_sal_ATGCAT -nt results/ATGCAT_fig_3.tiff \
   ]; then
  echo "======= Rebuilding bed_data_analysis_sal_ATGCAT ======="
  docker build -t bed_data_analysis_sal_atgcat:latest -f dockerfiles/Dockerfile.bed_data_analysis_sal_ATGCAT .

  docker run -v $(pwd)/results:/app/results:rw bed_data_analysis_sal_atgcat:latest
fi

if [ gene_expression_vis_salmonella.R -nt results/Fig_5_TPM_LFC_methyl_mean.tiff \
  -o dockerfiles/Dockerfile.gene_expression_vis_salmonella -nt results/Fig_5_TPM_LFC_methyl_mean.tiff \
   ]; then
  echo "======= Rebuilding gene_expression_vis_salmonella ======="
  docker build -t gene_expression_vis_salmonella:latest -f dockerfiles/Dockerfile.gene_expression_vis_salmonella .

  docker run -v $(pwd)/results:/app/results:rw gene_expression_vis_salmonella:latest
fi
if [ DEM_transcriptomic_annotation_GATC_fin.R -nt results/GATC_DM_upstream_id.txt \
  -o dockerfiles/Dockerfile.DEM_transcriptomic_annotation_GATC_fin -nt results/GATC_DM_upstream_id.txt \
   ]; then
  echo "======= Rebuilding DEM_transcriptomic_annotation_GATC_fin ======="
  docker build -t dem_transcriptomic_annotation_gatc_fin:latest -f dockerfiles/Dockerfile.DEM_transcriptomic_annotation_GATC_fin .

  docker run -v $(pwd)/results:/app/results:rw dem_transcriptomic_annotation_gatc_fin:latest
fi

if [ DEM_transcriptomic_annotation_CCWGG_fin.R -nt results/CCWGG_DM_upstream_id.txt \
  -o dockerfiles/Dockerfile.DEM_transcriptomic_annotation_CCWGG_fin -nt results/CCWGG_DM_upstream_id.txt \
   ]; then
  echo "======= Rebuilding DEM_transcriptomic_annotation_CCWGG_fin ======="
  docker build -t dem_transcriptomic_annotation_ccwgg_fin:latest -f dockerfiles/Dockerfile.DEM_transcriptomic_annotation_CCWGG_fin .

  docker run -v $(pwd)/results:/app/results:rw dem_transcriptomic_annotation_ccwgg_fin:latest
fi
if [ DEM_transcriptomic_annotation_ATGCAT_fin.R -nt results/ATGCAT_DM_upstream_id.txt \
  -o dockerfiles/Dockerfile.DEM_transcriptomic_annotation_ATGCAT_fin -nt results/ATGCAT_DM_upstream_id.txt \
   ]; then
  echo "======= Rebuilding DEM_transcriptomic_annotation_ATGCAT_fin ======="
  docker build -t dem_transcriptomic_annotation_atgcat_fin:latest -f dockerfiles/Dockerfile.DEM_transcriptomic_annotation_ATGCAT_fin .

  docker run -v $(pwd)/results:/app/results:rw dem_transcriptomic_annotation_atgcat_fin:latest
fi

if [ combined_violins.R -nt results/All_fig_6.tiff \
  -o dockerfiles/Dockerfile.combined_violins -nt results/All_fig_6.tiff \
   ]; then
  echo "======= Rebuilding combined_violins ======="
  docker build -t combined_violins:latest -f dockerfiles/Dockerfile.combined_violins .

  docker run -v $(pwd)/results:/app/results:rw combined_violins:latest
fi

if [ GATC_comparison_fin_intergenic.R -nt results/all_vs_all_fig_v_figure_8_AB.tiff \
  -o dockerfiles/Dockerfile.GATC_comparison_fin_intergenic -nt results/ATGCAT_transcriptomic_figure_7.tiff \
   ]; then
  echo "======= Rebuilding intergenic GATC literature data comparisons ======="
  docker build -t gatc_comparison_fin_intergenic:latest -f dockerfiles/Dockerfile.GATC_comparison_fin_intergenic .

  docker run -v $(pwd)/results:/app/results:rw gatc_comparison_fin_intergenic:latest
fi
if [ NanoViz_fin.R -nt results/Fig_8_D.tiff \
  -o dockerfiles/Dockerfile.NanoViz_fin -nt results/Fig_8_D.tiff \
   ]; then
  echo "======= Rebuilding heterogeneity example ======="
  docker build -t nanoviz_fin:latest -f dockerfiles/Dockerfile.NanoViz_fin .

  docker run -v $(pwd)/results:/app/results:rw nanoviz_fin:latest
fi
