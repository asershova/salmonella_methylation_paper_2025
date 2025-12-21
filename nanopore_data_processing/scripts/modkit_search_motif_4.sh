#!/bin/bash

#SBATCH --job-name=motif
#SBATCH --mem=20G
#SBATCH --cpus-per-task=40
#SBATCH --time=0-20:00
#SBATCH --output=motif-4.log
#SBATCH --mail-user=asershova@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=motif-4.err



module load StdEnv/2023 samtools/1.18  

for i in /home/ershovaa/scratch/14_12_23_14_barcodes/14_lib_bac_meth_all_mod/salmonella/*.bed; do /home/ershovaa/software/dist_modkit_v0.4.4_251055f/modkit motif search --log-filepath ${i}.motif.log -o ${i}.motifs.tsv -t 40 --min-sites 150 --in-bedmethyl $i --ref /home/ershovaa/scratch/references/salmonella/ST4-74.fa 




