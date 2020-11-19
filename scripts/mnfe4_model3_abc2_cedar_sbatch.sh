#!/bin/bash
#SBATCH --time=27-23:15:00
#SBATCH --account=def-er461912
#SBATCH --job-name=abc_july12_array
#SBATCH --output=%x-%j.out
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=500M
#SBATCH --array=1-1000

python abc_sampling_mnfe_setup.py -n '_july12_'$SLURM_ARRAY_TASK_ID
