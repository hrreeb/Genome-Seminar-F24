#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=pacbio_samtools_view_%J_stdout.txt
#SBATCH --error=pacbio_samtools_view_%J_stderr.txt
#SBATCH --job-name=pacbio_samtools_view
# 

bash /home/biol726311/BIOL7263_Genomics/scripts/pseudomonas/minimap/pacbio_sam_view.sh