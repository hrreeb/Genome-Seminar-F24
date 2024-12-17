#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem 16G
#SBATCH --output=blastp_%J_stdout.txt
#SBATCH --error=blastp_%J_stderr.txt
#SBATCH --job-name=blastp
# 

bash /scratch/biol726311/gouldian/annotate/diamond_blastp.sh