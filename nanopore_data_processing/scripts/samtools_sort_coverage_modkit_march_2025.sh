#!/bin/bash

#SBATCH --job-name=samtools
#SBATCH --mem=20G
#SBATCH --cpus-per-task=40
#SBATCH --time=0-24:00
#SBATCH --output=samtools_mar_25_sal.log
#SBATCH --mail-user=asershova@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=samtools_mar_25_sal.err



module load StdEnv/2023 samtools/1.18  

for i in /home/ershovaa/scratch/14_12_23_14_barcodes/14_lib_bac_meth_all_mod/salmonella/*.bam; do samtools sort -@ 40 $i > ${i%%bam}sorted.bam; samtools coverage ${i%%bam}sorted.bam > ${i%%bam}coverage;samtools index -@ 40 ${i%%bam}sorted.bam;done




