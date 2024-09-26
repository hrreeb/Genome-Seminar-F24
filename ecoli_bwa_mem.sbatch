#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=4
#SBATCH --mem 16G
#SBATCH --output=ecoli_bwa_%J_stdout.txt
#SBATCH --error=ecoli_bwa_%J_stderr.txt
#SBATCH --job-name=ecoli_bwa
# 

bash /home/biol726311/BIOL7263_Genomics/scripts/BWA/ecoli_bwa_mem.sh