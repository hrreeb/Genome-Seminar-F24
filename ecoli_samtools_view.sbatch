#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=ecoli_view_%J_stdout.txt
#SBATCH --error=ecoli_view_%J_stderr.txt
#SBATCH --job-name=ecoli_view
#

bash /home/biol726311/BIOL7263_Genomics/scripts/BWA/ecoli_samtools_view.sh