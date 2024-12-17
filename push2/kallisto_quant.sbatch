#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem 24G
#SBATCH --output=kallisto_%J_stdout.txt
#SBATCH --error=kallisto_%J_stderr.txt
#SBATCH --job-name=kallisto
#SBATCH --array=1-8
# 

bash kallisto_quant.sh $(sed -n "${SLURM_ARRAY_TASK_ID}p" kallisto_quant.args)