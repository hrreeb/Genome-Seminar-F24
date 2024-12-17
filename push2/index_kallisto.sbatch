#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem 24G
#SBATCH --output=kallisto_ind_%J_stdout.txt
#SBATCH --error=kallisto_ind_%J_stderr.txt
#SBATCH --job-name=kallisto_ind
# 

bash /scratch/biol726311/gouldian/kallisto/index_kallisto.sh