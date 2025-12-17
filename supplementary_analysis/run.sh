#!/bin/bash
#run all scripts to produce all results
if [ scripts/bed_data_analysis_GAATTC_ROC_salmonella.R -nt results/GAATTC_fig_S1.1.tiff \
  -o dockerfiles/Dockerfile.bed_data_analysis_gaattc_roc_salmonella -nt results/GAATTC_fig_S1.1.tiff \
   ]; then
  echo "======= Rebuilding ROC curve based on GAATTC data ======="
  docker build -t bed_data_analysis_gaattc_roc_salmonella:latest -f dockerfiles/Dockerfile.bed_data_analysis_GAATTC_ROC_salmonella .

  docker run -v $(pwd)/results:/app/results:rw bed_data_analysis_gaattc_roc_salmonella:latest
fi
