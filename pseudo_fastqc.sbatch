#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 8G
#SBATCH --output=pseudo_fastqc_%J_stdout.txt
#SBATCH --error=pseudo_fastqc_%J_stderr.txt
#SBATCH --job-name=ecoli_fastqc
#

bash /home/biol726311/BIOL_Genomics/scripts/fastqcPseudo/pseudo_fastqc.sh