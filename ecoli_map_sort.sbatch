#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=4
#SBATCH --mem 16G
#SBATCH --output=ecoli_map_sort_%J_stdout.txt
#SBATCH --error=ecoli_map_sort_%J_stderr.txt
#SBATCH --job-name=ecoli_map_sort
#

bash /home/biol726311/BIOL7263_Genomics/scripts/BWA/ecoli_map_sort.sh