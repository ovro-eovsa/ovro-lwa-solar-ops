#!/bin/bash
#SBATCH --job-name=refcaldaily
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/lustre/solarpipe/slurmlog/%j.out
#SBATCH --error=/lustre/solarpipe/slurmlog/%j.err
#SBATCH --mail-user=pz47@njit.edu


DIRSOFT=/opt/devel/solarpipe/envs/lwasolarpipe/
DIR_PY_ENV=/opt/devel/bin.chen/envs/suncasa/
source /home/solarpipe/.bashrc
conda activate $DIR_PY_ENV

# add DIRSOFT to the python path
export PYTHONPATH=$DIRSOFT:$PYTHONPATH

cd /lustre/solarpipe/

srun $DIR_PY_ENV/bin/python $DIRSOFT/daily_refra_corr.py 