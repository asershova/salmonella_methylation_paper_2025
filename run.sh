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

if [ bed_data_analysis_sal_ATGCAT_7_09_2.R -nt results/ATGCAT_jitter.tiff \
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

