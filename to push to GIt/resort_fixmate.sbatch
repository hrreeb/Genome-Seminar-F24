#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=resort_%J_stdout.txt
#SBATCH --error=resort_%J_stderr.txt
#SBATCH --job-name=resort
#SBATCH --array=1-8
#SBATCH --time=00:50:00
# 

bash resort_fixmate.sh $(sed -n "${SLURM_ARRAY_TASK_ID}p" resort_fixmate.args)