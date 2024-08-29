#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 16G
#SBATCH --output=[job_name]_%J_stdout.txt
#SBATCH --error=[job_name]_%J_stderr.txt
#SBATCH --job-name=[job_name]
# 

bash test.sh