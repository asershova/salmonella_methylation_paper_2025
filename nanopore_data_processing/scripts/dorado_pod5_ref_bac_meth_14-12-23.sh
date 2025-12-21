#!/bin/bash

#SBATCH --job-name=dorado2
#SBATCH --gpus-per-node=h100:2 
#SBATCH --mem=20G
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --time=0-11:00
#SBATCH --output=dorado_bac_meth_14_12_all_mod.log
#SBATCH --mail-user=asershova@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=dorado_bac_meth_14_12_all_mod.err



module load StdEnv/2023 dorado/0.8.3  

nvidia-smi


tar -xvzf /home/ershovaa/scratch/14_12_23_14_barcodes.tar.gz --directory $SLURM_TMPDIR
#demultiplex
mkdir $SLURM_TMPDIR/14_lib_bac_meth_all_mod
dorado basecaller --modified-bases-models /home/ershovaa/software/dorado_models_feb_2025/dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v2,/home/ershovaa/software/dorado_models_feb_2025/dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v2 /home/ershovaa/software/dorado_models_feb_2025/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 $SLURM_TMPDIR/14_12_23_14_barcodes/14-12-23-14-barcodes/no_sample/20231214_1515_MN39980_FAX99382_af705a4e/pod5/ | dorado demux --kit-name SQK-NBD114-24 --output-dir $SLURM_TMPDIR/14_lib_bac_meth_all_mod

cp -r $SLURM_TMPDIR/14_lib_bac_meth_all_mod /home/ershovaa/scratch/14_12_23_14_barcodes/


