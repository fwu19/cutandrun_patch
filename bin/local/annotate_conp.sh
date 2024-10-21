#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --mem=12G
#SBATCH -n 1
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=fwu@fredhutch.org 
#SBATCH --output=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A_%a.%x.out 
#SBATCH --error=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A_%a.%x.err

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
shDir=$(realpath $(dirname $SCRIPT_PATH))

## Source params and check whether all required variables are defined
json=$( realpath $1 ); shift
pubDir=$( cat $json | jq '.outdir' | sed "s/\"//g; s/\'//g" )
genome=$( cat $json | jq '.genome' | sed "s/\"//g; s/\'//g" )
gtf4homer=$( cat $json | jq '.gtf4homer' | sed "s/\"//g; s/\'//g" )
subdir=$1; shift # under $pubDir

# IFS='\n' read -r -a peaks <<< "$(find $(realpath $pubDir)/$subdir/07_consensus_peaks/ | egrep "*.bed$" )"
# npeak=${#peaks[@]}

slurm="--parsable --mem=24G -n 6 --output=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.out --error=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.err --nodes=1"
for peak in $(find $(realpath $pubDir)/$subdir/07_consensus_peaks/ | egrep "*.bed$"); do
    echo -e Annotate consensus peak,$(sbatch $slurm $shDir/annotate_bed.sh $peak $genome $gtf4homer),$(basename $peak)
done

