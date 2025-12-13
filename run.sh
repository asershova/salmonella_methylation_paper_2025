#!/bin/bash
#run all scripts to produce all results
docker build -t bed_data_analysis_sal_GATC:latest -f dockerfiles/Dockerfile.bed_data_analysis_sal_GATC_28_08 .

docker run -v $(pwd)/results:/app/results:rw bed_data_analysis_sal_GATC:latest
