#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=indexed_markdup_%J_stdout.txt
#SBATCH --error=indexed_markdup_%J_stderr.txt
#SBATCH --job-name=indexed_markdup
#SBATCH --array=1-8
#SBATCH --time=00:50:00
# 

bash indexed_markdup.sh $(sed -n "${SLURM_ARRAY_TASK_ID}p" indexed_markdup.args)