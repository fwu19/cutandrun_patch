#!/usr/bin/env bash

set -uo pipefail

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

json=params.json

## run nf-core nextflow workflow for trimming and alignments
cut -f 1-5 -d , $input >nfcore_input.csv

module purge all
module load Nextflow/23.04.2
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH


nextflow run nf-core/cutandrun -r 2.0 -resume \
-with-report nextflow.report.nf-core.html \
-params-file $json -c $shdir/nextflow.singularity.config \
--input nfcore_input.csv \
--igenomes_base /shared/biodata/reference/iGenomes/ --normalisation_mode None --use_control false --save_spikein_aligned --save_align_intermed --trim_nextseq 20 --skip_multiqc --skip_heatmaps --skip_igv --skip_frip $@
