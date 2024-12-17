#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=4
#SBATCH --mem 16G
#SBATCH --output=pseud_bwa_mem_%J_stdout.txt
#SBATCH --error=pseud_bwa_mem_%J_stderr.txt
#SBATCH --job-name=pseud_bwa_mem
# 

bash /home/biol726311/scripts/pseudomonas/bwa/illumina_bwa_mem.sh