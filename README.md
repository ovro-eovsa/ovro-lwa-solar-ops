# ovro-lwa-solar-ops

Scripts and codes related to operations of OVRO-LWA for solar studies

# Realtime Pipeline

## Slurm managed run


By default, realtime pipeline runs on 10 nodes from calim cluster, every node requires 12 cpus and 32G Mem.

```bash
#!/bin/bash
#SBATCH --job-name=solar-test
#SBATCH --partition=general
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=14:00:00
#SBATCH --output=/lustre/solarpipe/slurmlog/%j.out
#SBATCH --error=/lustre/solarpipe/slurmlog/%j.err
#SBATCH --mail-user=pz47@njit.edu
```

To start (restart). First check the existing running job with `sinfo`, then `scancel` if necessary. Then:

```bash
sbatch runSlurm_solarPipeline.sh slow # for slow pipeline
sbatch runSlurm_solarPipeline.sh fast # for fast pipeline
```


To run slow imaging for given period of time, do:

```bash
sbatch runSlurm_solarPipeline.sh testslowfixedtime \
 "2024-09-15T20:00:00"  "2024-09-15T21:00:00" 
```