#!/bin/bash
#
#SBATCH --partition=debug
#SBATCH --ntasks=4
#SBATCH --mem 16G
#SBATCH --output=pfam_%J_stdout.txt
#SBATCH --error=pfam_%J_stderr.txt
#SBATCH --job-name=pfam
# 

bash /home/biol726311/BIOL7263_Genomics/scripts/unmapped/unmapped_pfam.sh