#!/bin/bash
#SBATCH --job-name=solarpipedaily
#SBATCH --partition=general
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --time=16:00:00
#SBATCH --output=/lustre/solarpipe/slurmlog/%j.out
#SBATCH --error=/lustre/solarpipe/slurmlog/%j.err
#SBATCH --mail-user=pz47@njit.edu


DIRSOFT=/lustre/peijin/ovro-lwa-solar-ops/
DIRRUN=/lustre/peijin/testslurm/ # for no realtime test
DIR_PY_ENV=/opt/devel/bin.chen/envs/suncasa/
source /home/solarpipe/.bashrc
conda activate $DIR_PY_ENV

# add DIRSOFT to the python path
export PYTHONPATH=$DIRSOFT:$PYTHONPATH

cd /lustre/solarpipe/

# run according to the case:
case "$1" in
    testnodes)
        srun $DIR_PY_ENV/bin/python slurm_taskid_test.py
        ;;
    testslowfixedtime)
        srun $DIR_PY_ENV/bin/python $DIRSOFT/solar_realtime_pipeline.py \
            --briggs -1.0 --slowfast slow --interval 600 --delay 180 --save_allsky \
            --start_time $2 --end_time $3 \
            --save_dir /lustre/solarpipe/test_realtime/
        ;;
    slow)
        srun $DIR_PY_ENV/bin/python $DIRSOFT/solar_realtime_pipeline.py \
        --briggs -1.0 --slowfast slow --interval 600 --delay 180 --save_allsky --no_refracorr --slurm_kill_after_sunset
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
            --start_time 2024-09-22T13:30:00 --end_time 2024-09-23T00:30:00
        ;;
    slownorealtime)
        srun $DIR_PY_ENV/bin/python $DIRSOFT/solar_realtime_pipeline.py \
            --briggs -1.0 --slowfast slow --interval 600 --delay 180  --save_allsky \
            --start_time 2024-09-22T16:00:00 --end_time 2024-09-23T00:30:00
        ;;
    fastnorealtime)
        srun $DIR_PY_ENV/bin/python $DIRSOFT/solar_realtime_pipeline.py \
            --briggs 1.0 --slowfast fast --interval 100 --delay 180 \
            --start_time 2024-09-22T16:00:00 --end_time 2024-09-23T00:30:00
        ;;
    *)
        echo "Usage: sbatch $0 {test|slow|fast}"
        exit 1
        ;;
esac