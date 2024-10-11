#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=4
#SBATCH --mem 16G
#SBATCH --output=quali_rep_%J_stdout.txt
#SBATCH --error=quali_rep_%J_stderr.txt
#SBATCH --job-name=quali_rep
# 

bash /home/biol726311/BIOL7263_Genomics/scripts/BWA/quali_rep.sh