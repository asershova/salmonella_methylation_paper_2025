#!/bin/bash

#SBATCH --job-name=dorado_aln
#SBATCH --mem=20G
#SBATCH --cpus-per-task=40
#SBATCH --time=0-24:00
#SBATCH --output=dorado_aln_14_sal.log
#SBATCH --mail-user=asershova@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=dorado_aln_14_sal.err



module load StdEnv/2023 dorado/0.8.3  

mkdir /home/ershovaa/scratch/14_12_23_14_barcodes/14_lib_bac_meth_all_mod/salmonella/aligned
dorado aligner /home/ershovaa/scratch/references/salmonella/ST4-74.fa /home/ershovaa/scratch/14_12_23_14_barcodes/14_lib_bac_meth_all_mod/salmonella --output-dir /home/ershovaa/scratch/14_12_23_14_barcodes/14_lib_bac_meth_all_mod/salmonella/aligned --emit-summary --threads 40





