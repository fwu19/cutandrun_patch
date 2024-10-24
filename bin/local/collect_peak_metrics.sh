#!/usr/bin/env bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fwu@fredhutch.org
#SBATCH --output=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.out
#SBATCH --error=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.err
#SBATCH --nodes=1
#SBATCH --mem=24G
#SBATCH -n 1

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

## Load required modules
module purge all
module load fhR/4.1.2-foss-2021b
export R_LIBS=/fh/fast/_SR/Genomics/user/fwu/R/x86_64-pc-linux-gnu-library/4.1 # 

## Source params and check whether all required variables are defined
json=$( realpath $1 ); shift
pubDir=$( cat $json | jq '.outdir' | sed "s/\"//g; s/\'//g" )
minReps=$( cat $json | jq '.minReplicates' | sed "s/\"//g; s/\'//g" )

subdir=$1; shift # under $pubDir

outdir=$(realpath $pubDir)/$subdir/04_reporting/qc
[[ -d $outdir ]] || mkdir -p $outdir
mkdir -p $outdir
cd $outdir
[[ -d _work ]] || mkdir -p _work
cd _work
ln -s $pubDir/$subdir/04_reporting/qc/read_metrics.csv .
ln -s $pubDir/$subdir/03_peak_calling/MACS2/*.*Peak .
ln -s $pubDir/$subdir/03_peak_calling/SEACR/*/*.bed .
ln -s $pubDir/$subdir/04_reporting/qc/reads_in_peak/*.readsInPeak.csv .

Rscript $shDir/original_peak_metrics.r # read_metrics.csv peak_files *.readsInPeak.csv
Rscript $shDir/replicated_peaks.r $minReps # original_peaks.rds original_peak_metrics.csv
Rscript $shDir/consensus_peaks.r # replicated_peaks.rds replicated_peak_metrics.csv


mv original_peaks.rds original_peak_metrics.csv original_peak_widths.rds ../
mv replicated_peaks.rds replicated_peak_metrics.csv ../
mv 06_replicated_peaks ../../../
mv consensus_peaks.rds consensus_peak_metrics.csv ../
mv 07_consensus_peaks ../../../
cd ..; rm -r _work
