#!/bin/bash
#
#SBATCH --partition=debug
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=orf_%J_stdout.txt
#SBATCH --error=orf_%J_stderr.txt
#SBATCH --job-name=orf
# 

bash /home/biol726311/BIOL7263_Genomics/scripts/unmapped/unmapped_orf.sh