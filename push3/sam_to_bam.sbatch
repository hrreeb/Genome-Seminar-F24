#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=samtobam_%J_stdout.txt
#SBATCH --error=samtobam_fullmap_%J_stderr.txt
#SBATCH --job-name=samtobam_fullmap
#SBATCH --array=1-8
#SBATCH --time=00:50:00
# 

bash sam_to_bam.sh $(sed -n "${SLURM_ARRAY_TASK_ID}p" sam_to_bam.args)