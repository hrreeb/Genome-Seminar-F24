#!/bin/bash
#
#SBATCH --partition=debug
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=blastn_%J_stdout.txt
#SBATCH --error=blastn_%J_stderr.txt
#SBATCH --job-name=blastn
# 

bash /home/biol726311/BIOL7263_Genomics/scripts/unmapped/unmapped_blastn.sh