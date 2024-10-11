#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=ecoli_call_%J_stdout.txt
#SBATCH --error=ecoli_call_%J_stderr.txt
#SBATCH --job-name=ecoli_call
# 

bash /home/mbtoomey/BIOL7263_Genomics/scripts/BWA/ecoli_filt.sh