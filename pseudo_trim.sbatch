#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem 24G
#SBATCH --output=pseudoTrim_%J_stdout.txt
#SBATCH --error=pseudoTrim_%J_stderr.txt
#SBATCH --job-name=trimPseudo
# 

bash /scratch/biol726311/pseudomonas_gm41/pseudo_trim.sh