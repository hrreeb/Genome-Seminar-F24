#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=hisat2_fullmap_%J_stdout.txt
#SBATCH --error=hisat2_fullmap_%J_stderr.txt
#SBATCH --job-name=hisat2_fullmap
#SBATCH --array=1-8
#SBATCH --time=00:50:00
# 

bash hisat2_full_mapping.sh $(sed -n "${SLURM_ARRAY_TASK_ID}p" hisat2_full_mapping.args)