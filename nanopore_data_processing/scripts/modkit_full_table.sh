#!/bin/bash

#SBATCH --job-name=motif_0
#SBATCH --mem=20G
#SBATCH --cpus-per-task=12
#SBATCH --time=0-12:00
#SBATCH --output=m6A_full.log
#SBATCH --mail-user=asershova@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=m6A_full.err



module load StdEnv/2023 samtools/1.18  

modkit="/home/ershovaa/software/dist_modkit_v0.4.4_251055f/modkit"
sal_dir="/home/ershovaa/scratch/14_12_23_14_barcodes/14_lib_bac_meth_all_mod/salmonella/aligned/"
ref="/home/ershovaa/scratch/references/salmonella/ST4-74.fa"
for bam in ${sal_dir}*sorted.bam; do
 ${modkit} extract full ${bam} ${bam}.modkit.bgz --bgzf --threads 12 --out-threads 12
done

