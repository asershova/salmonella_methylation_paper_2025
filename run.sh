#!/bin/bash
#run all scripts to produce all results
#docker build -t bed_data_analysis_sal_gatc:latest -f dockerfiles/Dockerfile.bed_data_analysis_sal_GATC_28_08 .

#docker run -v $(pwd)/results:/app/results:rw bed_data_analysis_sal_gatc:latest

docker build -t bed_data_analysis_sal_ccwgg:latest -f dockerfiles/Dockerfile.bed_data_analysis_sal_CCWGG_7_09_2 .

docker run -v $(pwd)/results:/app/results:rw bed_data_analysis_sal_ccwgg:latest

