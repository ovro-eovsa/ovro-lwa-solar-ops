#!/bin/bash
#SBATCH --job-name=solar-test
#SBATCH --partition=general
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=336:00:00
#SBATCH --output=/home/solarpipe/log/%j.out
#SBATCH --error=/home/solarpipe/log/%j.err
#SBATCH --mail-user=pz47@njit.edu


DIRSOFT=/lustre/peijin/ovro-lwa-solar-ops/
DIRRUN=/lustre/peijin/testslurm/
DIR_PY_ENV=/opt/devel/peijin/solarenv    
source /home/peijinz/.bashrc
conda activate $DIR_PY_ENV


# run according to the case:
case "$1" in
    testnodes)
        srun $DIR_PY_ENV/bin/python slurm_taskid_test.py
        ;;
    slowtest)
        srun $DIR_PY_ENV/bin/python $DIRSOFT/solar_realtime_pipeline.py \
            --briggs -1.0 --slowfast slow --interval 600 --delay 180 --save_allsky \
            --start_time 2024-09-15T20:00:00 --end_time 2024-09-15T20:30:00 \
            --save_dir /lustre/solarpipe/test_realtime/
        ;;
    fast)
        srun $DIR_PY_ENV/bin/python $DIRSOFT/solar_realtime_pipeline.py \
            --briggs 1.0 --slowfast fast --interval 100 --delay 180
        ;;
    testnorealtime)
        srun $DIR_PY_ENV/bin/python $DIRSOFT/solar_realtime_pipeline.py \
            --briggs 1.0 --slowfast slow --interval 300 --delay 180 \
            --proc_dir /fast/peijinz/slurmtest/ --save_dir $DIRRUN/save/ \
            --logger_dir $DIRRUN/log/ \
            --start_time 2024-09-15T20:00:00 --end_time 2024-09-15T20:30:00
        ;;
    *)
        echo "Usage: sbatch $0 {test|slow|fast}"
        exit 1
        ;;
esac