#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=map_stats_%J_stdout.txt
#SBATCH --error=map_stats_%J_stderr.txt
#SBATCH --job-name=map_stats
#SBATCH --array=1-8
#SBATCH --time=00:50:00
# 

bash map_stats.sh $(sed -n "${SLURM_ARRAY_TASK_ID}p" map_stats.args)