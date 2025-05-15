#!/bin/bash
#SBATCH --job-name=solarpipedaily
#SBATCH --partition=solar
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=4
#SBATCH --distribution=cyclic
#SBATCH --nodelist=lwacalim[05-09]
#SBATCH --mem=38G
#SBATCH --time=16:00:00
#SBATCH --output=/lustre/solarpipe/slurmlog/%j.out
#SBATCH --error=/lustre/solarpipe/slurmlog/%j.err
#SBATCH --mail-user=pz47@njit.edu


######## SBATCH --ntasks-per-node=2

DIRSOFT=/lustre/peijin/ovro-lwa-solar-ops/
DIRRUN=/lustre/peijin/testslurm/ # for no realtime test
DIR_PY_ENV=/opt/devel/peijin/solarenv    #/opt/devel/bin.chen/envs/suncasa/
CLEAR_CACHE_BEFORE_RUN=True

source /home/solarpipe/.bashrc
conda activate $DIR_PY_ENV

# add DIRSOFT to the python path
export PYTHONPATH=$DIRSOFT:$PYTHONPATH

cd /lustre/solarpipe/

# clear cache before run
if [ "$CLEAR_CACHE_BEFORE_RUN" = "True" ]; then
    echo "Clearing cache before run"
    rm -rf /dev/shm/srtmp/*
fi

# run according to the case:
case "$1" in
    slow)
        srun $DIR_PY_ENV/bin/python $DIRSOFT/solar_realtime_pipeline.py \
        --briggs -0.5 --slowfast slow --interval 200 --delay 180 --save_allsky \
        --no_refracorr --slurm_kill_after_sunset --keep_working_fits --save_selfcaltab --use_jpl_ephem
        ;;
    fast)
        srun $DIR_PY_ENV/bin/python $DIRSOFT/solar_realtime_pipeline.py \
            --briggs 1.0 --slowfast fast --interval 3000 --delay 180 --use_jpl_ephem
        ;;
    slownorealtime)
        srun $DIR_PY_ENV/bin/python $DIRSOFT/solar_realtime_pipeline.py \
            --bands 23MHz 27MHz 32MHz 36MHz 41MHz 46MHz 50MHz 55MHz 59MHz 64MHz 69MHz 73MHz 78MHz 82MHz \
            --briggs -0.0 --slowfast slow --interval 200 --delay 180 --no_refracorr --alt_limit 0 \
            --start_time 2025-02-21T15:40:00 --end_time 2025-02-21T17:20:00 --save_selfcaltab --use_jpl_ephem
        ;;   
    slownorealtimecostumizeddir)
        srun $DIR_PY_ENV/bin/python $DIRSOFT/solar_realtime_pipeline.py  \
            --briggs -0.5 --slowfast slow --interval 100 --delay 180 --no_refracorr\
            --nolustre --file_path /lustre/solarpipe/20230927/slow/ \
            --alt_limit 0 \
            --bands 23MHz 27MHz 32MHz 36MHz 41MHz 46MHz 50MHz 55MHz 59MHz 64MHz 69MHz 73MHz 78MHz 82MHz \
            --calib_file 20230925_052635 \
            --start_time 2023-09-27T16:45:05 --end_time 2023-09-27T21:39:50 --save_selfcaltab
        ;;
    fastnorealtime)
        srun $DIR_PY_ENV/bin/python $DIRSOFT/solar_realtime_pipeline.py \
            --briggs 1.0 --slowfast fast --interval 100 --delay 180 \
            --start_time 2024-09-22T20:00:00 --end_time 2024-09-23T00:30:00
        ;;
    testnodes)
        srun $DIR_PY_ENV/bin/python slurm_taskid_test.py
        ;;
    testnorealtime)
        srun $DIR_PY_ENV/bin/python $DIRSOFT/solar_realtime_pipeline.py \
            --briggs 0.0 --slowfast slow --interval 100 --delay 180 \
            --proc_dir /fast/peijinz/slurmtest/ --save_dir $DIRRUN/save/ \
            --logger_dir $DIRRUN/log/ \
            --start_time 2024-10-07T19:50:00 --end_time 2024-10-07T22:00:00
        ;;
    testslowfixedtime)
        srun $DIR_PY_ENV/bin/python $DIRSOFT/solar_realtime_pipeline.py \
            --briggs -1.0 --slowfast slow --interval 600 --delay 180 --save_allsky \
            --start_time $2 --end_time $3 \
            --save_dir /lustre/solarpipe/test_realtime/
        ;;
    *)
        echo "Usage: sbatch $0 {test|slow|fast}"
        exit 1
        ;;
esac