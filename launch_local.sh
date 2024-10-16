#!/usr/bin/env bash

set -uo pipefail

json=params.json

# Nextflow Version
module purge all
module load Nextflow/23.04.2

# check if script is started via SLURM or bash
# if with SLURM: there variable '$SLURM_JOB_ID' will exist
# `if [ -n $SLURM_JOB_ID ]` checks if $SLURM_JOB_ID is not an empty string
if [[ -n $SLURM_JOB_ID ]];  then
    # check the original location through scontrol and $SLURM_JOB_ID
    SCRIPT_PATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
else
    # otherwise: started with bash. Get the real location.
    SCRIPT_PATH=$(realpath $0)
fi
shdir=$(dirname $SCRIPT_PATH)

nextflow run $shdir/main.nf -c $shdir/nextflow.config -params-file $json $@