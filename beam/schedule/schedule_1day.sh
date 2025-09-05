#!/bin/bash

# --- conda initialize ---
__conda_setup="$('/opt/devel/solarpipe/miniconda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)" || true
if [ $? -eq 0 ] && [ -n "${__conda_setup:-}" ]; then
    eval "$__conda_setup"
elif [ -f "/opt/devel/solarpipe/miniconda/etc/profile.d/conda.sh" ]; then
    . "/opt/devel/solarpipe/miniconda/etc/profile.d/conda.sh"
else
    export PATH="/opt/devel/solarpipe/miniconda/bin:$PATH"
fi
unset __conda_setup
# --- end conda initialize ---

cd /opt/devel/solarpipe/operation/ovro-lwa-solar-ops
conda activate deployment
python schedule_1day.py

