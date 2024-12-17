#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=markdup_%J_stdout.txt
#SBATCH --error=markdup_%J_stderr.txt
#SBATCH --job-name=markdup
#SBATCH --array=1-8
#SBATCH --time=00:50:00
# 

bash markdup.sh $(sed -n "${SLURM_ARRAY_TASK_ID}p" markdup.args)