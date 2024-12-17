#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=4
#SBATCH --mem 8G
#SBATCH --output=trim_gould%J_stdout.txt
#SBATCH --error=trim_gould_%J_stderr.txt
#SBATCH --job-name=trim_gould
# 

bash /scratch/biol726311/gouldian/trim_output/trimming_gould3.sh