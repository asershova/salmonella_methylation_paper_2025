#!/bin/bash
#run all scripts to produce all results
if [ NanoViz_fin.R -nt results/Fig_8_D.tiff \
  -o dockerfiles/Dockerfile.NanoViz_fin -nt results/Fig_8_D.tiff \
   ]; then
  echo "======= Rebuilding heterogeneity example ======="
  docker build -t nanoviz_fin:latest -f dockerfiles/Dockerfile.NanoViz_fin .

  docker run -v $(pwd)/results:/app/results:rw nanoviz_fin:latest
fi


