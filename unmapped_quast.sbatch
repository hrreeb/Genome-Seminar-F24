#!/bin/bash
#
#SBATCH --partition=debug
#SBATCH --ntasks=1
#SBATCH --mem 16
#SBATCH --output=quast_%J_stdout.txt
#SBATCH --error=quast_%J_stderr.txt
#SBATCH --job-name=quast
#

bash /home/biol726311/BIOL7263_Genomics/scripts/unmapped/unmapped_quast.sh