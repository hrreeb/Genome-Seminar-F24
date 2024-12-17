#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=pseud_samtools_sort_%J_stdout.txt
#SBATCH --error=pseud_samtools_sort_%J_stderr.txt
#SBATCH --job-name=pseudo_samtools_sort
# 

bash /home/biol726311/BIOL7263_Genomics/scripts/pseudomonas/bwa/illumina_sam_sort.sh