#!/bin/bash
#
#SBATCH --partition=64gb_24core
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem 48G
#SBATCH --output=pseud_quast_%J_stdout.txt
#SBATCH --error=pseud_quast_%J_stderr.txt
#SBATCH --job-name=pseud_quast
# 

bash /scratch/biol726311/pseudomonas_gm41/Assembly/pseudo_quast.sh