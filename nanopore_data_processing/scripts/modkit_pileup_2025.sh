#!/bin/bash

#SBATCH --job-name=modkit_pileup
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
#SBATCH --time=0-24:00
#SBATCH --output=sal_all_mod_modkit_pileup.log
#SBATCH --mail-user=asershova@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=sal_all_mod_modkit_pileup.err



module load StdEnv/2023 samtools/1.18  

for i in /home/ershovaa/scratch/14_12_23_14_barcodes/14_lib_bac_meth_all_mod/salmonella/aligned/*.bam; do samtools sort -@ 8 $i > ${i%%bam}sorted.bam; samtools coverage ${i%%bam}sorted.bam > ${i%%bam}coverage;samtools index -@ 8 ${i%%bam}sorted.bam; /home/ershovaa/software/dist_modkit_v0.4.4_251055f/modkit pileup --log-filepath ${i%%bam}sorted.log --threads 8 ${i%%bam}sorted.bam ${i%%bam}sorted.bed;done


