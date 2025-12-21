#!/bin/bash
#run all scripts to produce all results
if [ scripts/bed_data_analysis_GAATTC_ROC_salmonella.R -nt results/GAATTC_fig_S1.1.tiff \
  -o dockerfiles/Dockerfile.bed_data_analysis_gaattc_roc_salmonella -nt results/GAATTC_fig_S1.1.tiff \
   ]; then
  echo "======= Rebuilding ROC curve based on GAATTC data ======="
  docker build -t bed_data_analysis_gaattc_roc_salmonella:latest -f dockerfiles/Dockerfile.bed_data_analysis_GAATTC_ROC_salmonella .

  docker run -v $(pwd)/results:/app/results:rw bed_data_analysis_gaattc_roc_salmonella:latest
fi

if [ scripts/test_coverage_18_12_25.R -nt results/seq_depth_figure_S1.3.tiff \
  -o dockerfiles/Dockerfile.test_coverage_18_12_25 -nt results/seq_depth_figure_S1.3.tiff \
   ]; then
  echo "======= Rebuilding sequencing depth test ======="
  docker build -t test_coverage_18_12_25:latest -f dockerfiles/Dockerfile.test_coverage_18_12_25 .

  docker run -v $(pwd)/results:/app/results:rw test_coverage_18_12_25:latest
fi

if [ scripts/bed_data_analysis_sal_CAGAG_18_12_25_ab.R -nt results/CAGAG_figure_S2.2.AB.tiff \
  -o dockerfiles/Dockerfile.bed_data_analysis_sal_CAGAG_18_12_25_ab -nt results/CAGAG_figure_S2.2.AB.tiff \
   ]; then
  echo "======= Rebuilding CAGAG analysis ======="
  docker build -t bed_data_analysis_sal_cagag_18_12_25_ab:latest -f dockerfiles/Dockerfile.bed_data_analysis_sal_CAGAG_18_12_25_ab .

  docker run -v $(pwd)/results:/app/results:rw bed_data_analysis_sal_cagag_18_12_25_ab:latest
fi
if [ scripts/bed_data_analysis_sal_GATCAG_18_12_25_cd.R -nt results/GATCAG_figure_S2.2.CD.tiff \
  -o dockerfiles/Dockerfile.bed_data_analysis_sal_GATCAG_18_12_25_cd -nt results/GATCAG_figure_S2.2.CD.tiff \
   ]; then
  echo "======= Rebuilding GATCAG analysis ======="
  docker build -t bed_data_analysis_sal_gatcag_18_12_25_cd:latest -f dockerfiles/Dockerfile.bed_data_analysis_sal_GATCAG_18_12_25_cd .

  docker run -v $(pwd)/results:/app/results:rw bed_data_analysis_sal_gatcag_18_12_25_cd:latest
fi

if [ scripts/GATC_comparison_fin_intragenic.R -nt results/PMC9239280_vs_my_data_MEP_LSP_intragenic_fig_S3.1.png \
  -o dockerfiles/Dockerfile.GATC_comparison_fin_intragenic -nt results/PMC9239280_vs_my_data_MEP_LSP_intragenic_fig_S3.1.png \
   ]; then
   echo "======= Rebuilding GATC intragenic comparative analysis (three datasets) ======="
  docker build -t gatc_comparison_fin_intragenic:latest -f dockerfiles/Dockerfile.GATC_comparison_fin_intragenic .

  docker run -v $(pwd)/results:/app/results:rw gatc_comparison_fin_intragenic:latest
fi
if [ scripts/ATGCAT_comparison_fin.R -nt results/PMC9239280_vs_my_data_UM_ATGCAT_venn_intra_fig_S4.1D.png \
  -o dockerfiles/Dockerfile.ATGCAT_comparison_fin -nt results/PMC9239280_vs_my_data_UM_ATGCAT_venn_intra_fig_S4.1D.png \
   ]; then
   echo "======= Rebuilding ATGCAT comparative analysis (three datasets) ======="
  docker build -t atgcat_comparison_fin:latest -f dockerfiles/Dockerfile.ATGCAT_comparison_fin .

  docker run -v $(pwd)/results:/app/results:rw atgcat_comparison_fin:latest
fi

