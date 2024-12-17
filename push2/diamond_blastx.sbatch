#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem 16G
#SBATCH --output=blastx_%J_stdout.txt
#SBATCH --error=blastx_%J_stderr.txt
#SBATCH --job-name=blastx
# 

bash /scratch/biol726311/gouldian/annotate/diamond_blastx.sh