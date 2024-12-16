#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem 24G
#SBATCH --output=pseudoSpades_%J_stdout.txt
#SBATCH --error=pseudoSpades_%J_stderr.txt
#SBATCH --job-name=spadesPseudo
# 

bash /scratch/biol726311/pseudomonas_gm41/Assembly/pseudo_spades.sh