#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH -n 1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fwu@fredhutch.org
#SBATCH --output=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.out
#SBATCH --error=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.err

###################################################################################
#### Cut & Run
###################################################################################

## Print out help message if no arguments are provided
if [[ $# < 1 ]]; then
        echo "This script runs Cut&Run analysis."
        echo -e "USAGE: bash `basename $0` <params> <sample.csv> <path/to/input/bam/dir> <path/to/output/peak/dir>"
        exit 1;
fi

## Load required modules
module purge all
module load SAMtools/1.11-GCC-10.2.0
module load BEDTools/2.29.2-GCC-9.3.0

## Source params and check whether all required variables are defined
json=$( realpath $1 ); shift
nfDir=$( cat $json | jq '.indir' | sed "s/\"//g; s/\'//g" )
pubDir=$( cat $json | jq '.outdir' | sed "s/\"//g; s/\'//g" )
samples="$1"; shift # a file to store sample IDs, one per line
if [[ -f $samples ]]; then
	IFS=$'\n\r' read -d '' -r -a samplelist < $samples
else
	IFS=' ' read -r -a samplelist <<< "$samples"
fi	

IFS=',' read -r -a arr <<< "${samplelist[$(( ${SLURM_ARRAY_TASK_ID:-2} - 1 ))]}"   

subdir=$1; shift # under $pubDir

## cd working directory
bam=$(realpath $pubDir)/$subdir/02_alignment/bowtie2/target/${arr[12]}.target.markdup.bam

peakdir=$(realpath $pubDir)/$subdir/03_peak_calling/

outdir=$(realpath $pubDir)/$subdir/04_reporting/qc/reads_in_peak
[[ -d $outdir ]] || mkdir -p $outdir

find $peakdir/ | egrep "/${arr[12]}" | egrep "narrow_peaks.narrowPeak$|broad_peaks.broadPeak$|stringent.bed$|filtered" 2>/dev/null | while IFS='' read -r peak || [[ -n "$peak" ]]; do
  
		if [[ -f $peak && -f $bam ]]; then
			echo -e $(basename $peak ),$( bedtools intersect -a $bam -b $peak -u | samtools view -F 256 -f 64 -c - )
		fi
  
done >$outdir/${arr[12]}.readsInPeak.csv
