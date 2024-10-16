#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH -n 1

###################################################################################
#### Cut & Run
###################################################################################

## Print out help message if no arguments are provided
if [[ $# < 1 ]]; then
        echo "This script runs Cut&Run analysis."
        echo -e "USAGE: bash `basename $0` <params> <samples.csv> <path/to/input/bam/dir>"
        exit 1;
fi

## Load required modules
module load SAMtools/1.11-GCC-10.2.0

## Source params and check whether all required variables are defined
json=$( realpath $1 ); shift
nfDir=$( cat $json | jq '.indir' | sed "s/\"//g; s/\'//g" )
pubDir=$( cat $json | jq '.outdir' | sed "s/\"//g; s/\'//g" )

samples=$( realpath $1 ); shift
if [[ -f $samples ]]; then
	IFS=$'\n\r' read -d '' -r -a samplelist < $samples
else
	IFS=' ' read -r -a samplelist <<< "$samples"
fi	

subdir=$1; shift # under $pubDir


## cd working directory
indir=$pubDir/$subdir/02_alignment/bowtie2/target
outdir=$pubDir/$subdir/04_reporting/qc
[[ -d $outdir ]] || mkdir -p $outdir
cd $outdir
cp $pubDir/04_reporting/qc/read_metrics.csv .

## Process one sample
echo "id,bt2_aligned_short_fragments" >read_metrics.short_fragments.csv
for str in "${samplelist[@]}"; do
    IFS=',' read -r -a arr <<< "$str"   
    prefix=${arr[12]}
    tgt=${prefix}.target.markdup
    if [[ -f $indir/$tgt.bam ]]; then
        echo -e "${prefix}",$(samtools view -cf 64 -F 260 $indir/$tgt.bam)
    fi
done >>read_metrics.short_fragments.csv
    
