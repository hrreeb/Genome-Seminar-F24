#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=fixmate_%J_stdout.txt
#SBATCH --error=fixmate_%J_stderr.txt
#SBATCH --job-name=fixmate
#SBATCH --array=1-8
#SBATCH --time=00:50:00
# 

bash fixmate.sh $(sed -n "${SLURM_ARRAY_TASK_ID}p" fixmate.args)