#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=pseud_samtools_view_%J_stdout.txt
#SBATCH --error=pseud_samtools_view_%J_stderr.txt
#SBATCH --job-name=pseudo_samtools_view
# 

bash /home/biol726311/BIOL7263_Genomics/scripts/pseudomonas/bwa/illumina_sam_view.sh