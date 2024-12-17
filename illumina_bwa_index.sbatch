#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=bwa_index_%J_stdout.txt
#SBATCH --error=bwa_index%J_stderr.txt
#SBATCH --job-name=bwa_index
# 

bash /scratch/biol726311/pseudomonas_gm41/Assembly/bwa/illumina_bwa_index.sh
