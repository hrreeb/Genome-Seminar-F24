#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem 16G
#SBATCH --output=mkdb_%J_stdout.txt
#SBATCH --error=mkdb_%J_stderr.txt
#SBATCH --job-name=mkdb
# 

bash /scratch/biol736211/gouldian/annotate/diamond_mkdb.sh