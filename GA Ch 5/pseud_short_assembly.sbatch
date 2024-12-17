#!/bin/bash
#
#SBATCH --partition=64gb_24core
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem 48G
#SBATCH --output=pseud_short_%J_stdout.txt
#SBATCH --error=pseud_short_%J_stderr.txt
#SBATCH --job-name=pseud_short
#SBATCH --time=00:50:00
# 

bash /scratch/biol726311/pseudomonas_gm41/pseud_short_assembly.sh%20
