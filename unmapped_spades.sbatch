#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=spades_%J_stdout.txt
#SBATCH --error= =spades_%J_stderr.txt
#SBATCH --job-name=spades
#

bash /home/biol726311/BIOL7263_Genomics/scripts/unmapped/unmapped_spades.sh