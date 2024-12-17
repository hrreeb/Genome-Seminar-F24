#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=minimap_pseud_%J_stdout.txt
#SBATCH --error=minimap_pseud_%J_stderr.txt
#SBATCH --job-name=minimap_pseud
#

bash /home/biol726311/BIOL7263_Genomics/scripts/pseudomonas/minimap/minimap_pseud.sh