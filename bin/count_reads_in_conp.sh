#!/usr/bin/env bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fwu@fredhutch.org
#SBATCH --output=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.out
#SBATCH --error=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.err
#SBATCH --nodes=1
#SBATCH --mem=24G
#SBATCH -n 6

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
shDir=$(dirname $SCRIPT_PATH)

## Source params and check whether all required variables are defined
samples=$(realpath $1); shift
json=$( realpath $1 ); shift
pubDir=$( cat $json | jq '.outdir' | sed "s/\"//g; s/\'//g" )
subdir=$1; shift # under $pubDir

## Source params and check whether all required variables are defined
indir=$pubDir/$subdir/07_consensus_peaks/

for conp in $(ls $indir/*.bed); do
  outdir=$indir/featureCounts/$(basename "$conp" | sed 's/.bed$//')
  [[ -d "$outdir" ]] || mkdir -p "$outdir"

    cat "$conp" | awk 'BEGIN {FS=OFS="\t"} { $2=$2+1; print $1":"$2"-"$3, $1, $2, $3, ".", $4, $5 }' >"$outdir/$(basename $conp | sed 's/.bed$/.saf/' )"

done

sbatch --array=2-$(egrep -c "^" $samples) $shDir/featureCounts.sh $samples

