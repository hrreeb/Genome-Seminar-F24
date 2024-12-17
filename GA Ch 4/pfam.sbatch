#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=2
#SBATCH --mem 16G
#SBATCH --output=pfam_%J_stdout.txt
#SBATCH --error=pfam_%J_stderr.txt
#SBATCH --job-name=pfam
# 

bash /scratch/biol726311/COPYBIOL7263_Genomics/sequencing_data/ecoli/assembly/pfam.sh