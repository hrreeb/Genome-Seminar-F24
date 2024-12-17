#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=dump_sort_%J_stdout.txt
#SBATCH --error=dump_sort_%J_stderr.txt
#SBATCH --job-name=dump_sort
#SBATCH --array=1-3
# 

bash fastqc_array.sh $(sed -n "${SLURM_ARRAY_TASK_ID}p" fastqc_array.args)
