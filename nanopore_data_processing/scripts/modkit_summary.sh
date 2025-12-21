#!/bin/bash

#SBATCH --job-name=modkit
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --time=0-11:00
#SBATCH --output=modkit_sum.log
#SBATCH --mail-user=asershova@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=modkit_sum.err



module load StdEnv/2023 samtools/1.18  

for i in /home/ershovaa/scratch/14_12_23_14_barcodes/14_lib_bac_meth_all_mod/salmonella/aligned/*sorted.bam; do /home/ershovaa/software/modkit/modkit summary --no-sampling $i > ${i%%bam}.summary; done



