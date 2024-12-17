#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=qmap_de_novo_%J_stdout.txt
#SBATCH --error=qmap_de_novo_%J_stderr.txt
#SBATCH --job-name=qmap_de_novo
# 

bash /scratch/biol726311/COPYBIOL7263_Genomics/sequencing_data/ecoli/assembly/qmap_de_novo.sh