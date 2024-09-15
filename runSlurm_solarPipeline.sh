#!/bin/bash
#SBATCH --job-name=solar-pipeline
#SBATCH --partition=general
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=336:00:00
#SBATCH --output=/home/peijinz/log/%j.out
#SBATCH --error=/home/peijinz/log/%j.err
#SBATCH --mail-user=pz47@njit.edu

# py env
source /home/peijinz/.bashrc
conda activate /opt/devel/peijin/solarenv


cd /lustre/peijin/ovro-lwa-solar-ops

# Run the Python script using srun
srun python osp_task_handler.py 