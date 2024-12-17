#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=pseudo_long_quast%J_stdout.txt
#SBATCH --error=pseudo_long_quast%J_stderr.txt
#SBATCH --job-name=pseudo_long_quast
# 

bash /scratch/biol726311/pseudomonas_gm41/Assembly/pseudo_quast_long.sh