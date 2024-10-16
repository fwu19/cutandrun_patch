#!/usr/bin/env bash
# set -euo pipefail

scratch=/hpc/temp/_SR/Genomics/users/fwu/slurm
slurm="--mail-type=END,FAIL --mail-user=fwu@fredhutch.org --output=$scratch/%A_%a.%x.out --error=$scratch/%A_%a.%x.err"

json=params.json

input=$( cat $json | jq '.input' | sed "s/\"//g; s/\'//g" ) # nf-core result dir
input=${1:-$input}; shift

nfDir=$( cat $json | jq '.indir' | sed "s/\"//g; s/\'//g" ) # save nf-core results
nfDir=${1:-$nfDir}; shift
[[ -d $nfDir ]] || mkdir -p $nfDir
nfDir=$(realpath $nfDir)
ln -sf $nfDir .

pubDir=$( cat $json | jq '.outdir' | sed "s/\"//g; s/\'//g" ) # final output 
pubDir=${1:-$pubDir}; shift
[[ -d $pubDir ]] || mkdir -p $pubDir
pubDir=$(realpath $pubDir)
ln -sf $pubDir .

log=cutandrun.log
echo -e "###### $( date ) ######" | tee -a $log
echo -e "Run Cut&Run workflow on $ss" | tee -a $log


## local pipeline
local=$(sbatch $slurm --parsable --mem=5G -n 1 --job-name="cnrLocal" $(dirname $0)/launch_local.sh --input $input --indir $nfDir --outdir $pubDir -resume)
echo -e "local workflow", $local | tee -a $log

