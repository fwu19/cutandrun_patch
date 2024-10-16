#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --mem=24G
#SBATCH -n 6
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=fwu@fredhutch.org 
#SBATCH --output=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A_%a.%x.out 
#SBATCH --error=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A_%a.%x.err

###################################################################################
#### Cut & Run
###################################################################################

## Print out help message if no arguments are provided
if [[ $# < 1 ]]; then
        echo "This script runs Cut&Run analysis."
        echo -e "USAGE: bash `basename $0` <params> <samples.csv> <path/to/inpit/bam/dir> <path/to/output/peak/dir>"
        exit 1;
fi

## Load required modules
module purge all
module load MACS2/2.2.6-foss-2019b-Python-3.7.4
module load BEDTools/2.29.2-GCC-9.3.0

## Source params and check whether all required variables are defined
json=$( realpath $1 ); shift
nfDir=$( cat $json | jq '.indir' | sed "s/\"//g; s/\'//g" )
pubDir=$( cat $json | jq '.outdir' | sed "s/\"//g; s/\'//g" )
genomeSize=$( cat $json | jq '.genomeSize' | sed "s/\"//g; s/\'//g" )

samples=$( realpath $1 ); shift
if [[ -f $samples ]]; then
	IFS=$'\n\r' read -d '' -r -a samplelist < $samples
else
	IFS=' ' read -r -a samplelist <<< "$samples"
fi	
IFS=',' read -r -a arr <<< "${samplelist[$(( ${SLURM_ARRAY_TASK_ID:-2} - 1 ))]}"        

subdir=$1; shift # under $pubDir

chbam=$(realpath $pubDir)/$subdir/02_alignment/bowtie2/target/${arr[12]}.target.markdup.bam
inbam=$(realpath $pubDir)/mergedIgG/02_alignment/bowtie2/target/${arr[13]}.target.markdup.bam
type=${arr[6]}
prefix=${arr[12]}
outdir=$pubDir/$subdir/03_peak_calling/MACS2
[[ -d $outdir ]] || mkdir -p $outdir
cd $outdir

###### call narrow peaks using keep-dup all ######
if [[ $type =~ 'narrow' || $type =~ 'both' ]]; then
	log=$prefix.narrow_peaks.log

	## call peaks with no control
	echo -e "############\n####Call narrow peaks without control\n" >$log
	macs2 callpeak --name $prefix.noControl.narrow  --treatment $chbam --outdir ./  --format BAMPE --gsize $genomeSize --keep-dup all -q 0.05 2>&1 | tee $log
	
	if [[ -f $inbam ]]; then
		## call peaks with respect to IgG ######
		echo -e "############\n####Call narrow peaks with control\n" >>$log
		macs2 callpeak --name $prefix.toIgG.narrow  --treatment $chbam --control $inbam --outdir ./  --format BAMPE --gsize $genomeSize --keep-dup all -q 0.05 2>&1 | tee -a $log

		echo -e "############\n####Intersect peaks with and without control\n" >>$log
		bedtools intersect -a $prefix.toIgG.narrow_peaks.narrowPeak -b $prefix.noControl.narrow_peaks.narrowPeak -u | awk -F "\t" '{if ($7 > 2 && $9 > 2) {print} }' >$prefix.toIgG.narrow_peaks.filtered.narrowPeak

	fi

fi


###### call broad peaks using keep-dup all ######
if [[ $type =~ 'broad' || $type =~ 'both' ]]; then
	log=$prefix.broad_peaks.log

	echo -e "############\n####Call broad peaks without control\n" >$clog
	macs2 callpeak --name $prefix.noControl.broad  --treatment $chbam --outdir ./  --format BAMPE --gsize $genomeSize --keep-dup all --broad --broad-cutoff 0.1 -q 0.05 --nomodel --nolambda 2>&1 | tee $log

	if [[ -f $inbam ]]; then

		echo -e "############\n####Call broad peaks with control\n" >$clog
		macs2 callpeak --name $prefix.toIgG.broad  --treatment $chbam --control $inbam --outdir ./  --format BAMPE --gsize $genomeSize --keep-dup all --broad --broad-cutoff 0.1 -q 0.05 2>&1 | tee -a $log

		echo -e "############\n####Intersect peaks with and without control\n" >>$log
		bedtools intersect -a $prefix.noControl.broad_peaks.broadPeak -b $prefix.toIgG.broad_peaks.broadPeak -u | awk -F "\t" '{if ($7 > 2 && $9 > 2) {print} }'  >$prefix.noControl.broad_peaks.filtered.broadPeak
	
	fi

fi