#!/bin/bash
#SBATCH --job-name=solarpipedaily
#SBATCH --partition=solar
#SBATCH --ntasks=15 #10
#SBATCH --cpus-per-task=13
#SBATCH --distribution=cyclic
#SBATCH --nodelist=lwacalim[05-09]
#SBATCH --mem=399999M
#SBATCH --time=16:00:00
#SBATCH --output=/lustre/solarpipe/slurmlog/%j.out
#SBATCH --error=/lustre/solarpipe/slurmlog/%j.err
#SBATCH --mail-user=pz47@njit.edu


######## SBATCH --ntasks-per-node=2

DIRSOFT=/opt/devel/solarpipe/operation/ovro-lwa-solar-ops/
DIR_PY_ENV=/opt/devel/solarpipe/envs/lwasolarpipe   #/opt/devel/bin.chen/envs/suncasa/
#DIR_PY_ENV=/opt/devel/msurajit/envs/my_suncasa_env    #/opt/devel/bin.chen/envs/suncasa/
CLEAR_CACHE_BEFORE_RUN=True

source /home/solarpipe/.bashrc
conda activate
conda activate $DIR_PY_ENV

# add DIRSOFT to the python path
export PYTHONPATH=$DIRSOFT:$PYTHONPATH

cd /fast/solarpipe/realtime_pipeline/proc/

# clear cache before run
if [ "$CLEAR_CACHE_BEFORE_RUN" = "True" ]; then
    echo "Clearing cache before run"
    rm -rf /dev/shm/srtmp/*
    rm -rf /fast/solarpipe/realtime_pipeline/proc/*
    #find /dev/shm -maxdepth 1 -name '__KMP_REGISTERED_LIB_*' -delete 2>/dev/null
fi

#/dev/shm/srtmp/
#save_selfcaltab

#--bands 18MHz 23MHz 27MHz 32MHz 36MHz 41MHz 46MHz 50MHz 55MHz 59MHz 64MHz 69MHz 73MHz 78MHz 82MHz \            
# run according to the case:
case "$1" in
    slow)
        srun $DIR_PY_ENV/bin/python $DIRSOFT/solar_realtime_pipeline.py \
        --briggs -0.5 --slowfast slow --interval 300 --delay 180  \
        --proc_dir  /fast/solarpipe/realtime_pipeline/proc/ \
        --proc_dir_mem  /fast/solarpipe/realtime_pipeline/procmem/  \
        --no_refracorr --slurm_kill_after_sunset --alt_limit 15
        ;;
    fast)
        srun $DIR_PY_ENV/bin/python $DIRSOFT/solar_realtime_pipeline.py \
            --briggs 1.0 --slowfast fast --interval 3000 --delay 180
        ;;
    slownorealtime)
        srun $DIR_PY_ENV/bin/python $DIRSOFT/solar_realtime_pipeline.py \
            --briggs -0.5 --slowfast slow --interval 150 --delay 180 --no_refracorr --alt_limit 0 \
            --bands 23MHz 27MHz 32MHz 36MHz 41MHz 46MHz 50MHz 55MHz 59MHz 64MHz 69MHz 73MHz 78MHz 82MHz \
	    --proc_dir  /fast/solarpipe/realtime_pipeline/proc/ \
            --proc_dir_mem  /dev/shm/srtmp/   \
            --start_time 2026-01-03T17:23:55 --end_time 2026-01-03T17:54:50 --save_selfcaltab
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
