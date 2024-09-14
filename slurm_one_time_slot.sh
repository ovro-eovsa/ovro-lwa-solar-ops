#!/bin/bash
#SBATCH --job-name=solar_one_time_slot
#SBATCH --partition=general
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=480G
#SBATCH --time=336:00:00
#SBATCH --output=/data07/peijinz/slurmtest.out
#SBATCH --export=ALL
#SBATCH --mail-user=pz47@gmail.com

conda activate suncasa


# input should be a dir str in the format of 'yyyymmdd_HHMMSS'
python proc_one_time_slot.py ${SLURM_ARRAY_TASK_ID}