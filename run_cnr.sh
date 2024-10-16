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

## run nf-core nextflow workflow for trimming and alignments
nfcore=$(sbatch $slurm --parsable --mem=5G -n 1 --job-name="cnrNFCo" $(dirname $0)/launch_nfcore.sh --input $input --outdir $nfDir )
echo -e "NF-core workflow", $nfcore | tee -a $log

## local pipeline
local=$(sbatch $slurm --parsable --dependency=afterok:$nfcore --mem=5G -n 1 --job-name="cnrLocal" $(dirname $0)/launch_local.sh --input $input --indir $nfDir --outdir $pubDir )
echo -e "local workflow", $local | tee -a $log

## copy some nf-core results to $pubDir
rsync=$(sbatch $slurm --parsable --dependency=afterok:$nfcore --mem=12G -n 1 --job-name=cnrRsync --wrap="cd $PWD; rsync -vAXEWHhrLog --chmod=755 --no-compress --progress $nfDir/01_prealign/ $pubDir/01_prealign/; rsync -vAXEWHhrLog --chmod=755 --no-compress --progress --include='*/' --include='*.target.dedup.*' --include='*.target.markdup.*' --include='*.spikein.sorted.*' --include='*.bowtie2.log' --exclude='*' $nfDir/02_alignment/ $pubDir/02_alignment/" )
echo -e "copy nf-core results", $rsync | tee -a $log

## run MultiQC
multiqc=$(sbatch $slurm --parsable --dependency=afterok:$nfcore --mem=36G -n 6 --job-name=cnrMQC --wrap="cd $PWD; module purge all; module load MultiQC/1.21-foss-2023a; [[ -d $pubDir/04_reporting/ ]] || mkdir -p $pubDir/04_reporting/; multiqc $nfDir/01_prealign/ $nfDir/02_alignment/ -o $pubDir/04_reporting/multiqc/ --ignore-symlinks -f" )
echo -e "multiqc", $multiqc | tee -a $log

## render html report
#bash bin/render_rmd.sh cutandrun.Rmd $(date '+%y%m%d')_cutandrun.html