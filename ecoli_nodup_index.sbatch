#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=ecoli_nodup_index_%J_stdout.txt
#SBATCH --error=ecoli_nodup_index_%J_stderr.txt
#SBATCH --job-name=ecoli_nodup_index
#

bash /home/biol726311/BIOL7263_Genomics/scripts/BWA/ecoli_nodup_index.sh