#!/bin/bash
#
#SBATCH --partition=64gb_24core
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem 62G
#SBATCH --output=pseud_long_%J_stdout.txt
#SBATCH --error=pseud_long_%J_stderr.txt
#SBATCH --job-name=pseud_long
# 

bash /scratch/biol726311/pseudomonas_gm41/Assembly/pseud_long_assembly.sh