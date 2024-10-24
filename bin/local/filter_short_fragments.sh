#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH -n 1
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
        echo -e "USAGE: bash `basename $0` <params> <samples.csv> <path/to/input/bam/dir>"
        exit 1;
fi

## Load required modules
module load SAMtools/1.11-GCC-10.2.0

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
prefix=${arr[12]}
indir=$nfDir/02_alignment/bowtie2/target
outdir=$pubDir/$subdir/02_alignment/bowtie2/target
[[ -d $outdir ]] || mkdir -p $outdir
cd $outdir

## Process one sample
tgt=${prefix}.target.markdup
(samtools view -H $indir/$tgt.bam; samtools view -F 0x04 $indir/$tgt.bam | awk '$9 < 121 && $9 > -121' | egrep -v chrEBV ) | samtools view -bh -@ 6 >${tgt}.bam 

if [[ -f $indir/${prefix}.target.dedup.bam  ]]; then     
	tgt=${prefix}.target.dedup   
    (samtools view -H $indir/${tgt}.bam; samtools view -F 0x04 $indir/${tgt}.bam | awk '$9 < 121 && $9 > -121' | egrep -v chrEBV ) | samtools view -bh -@ 6 >${tgt}.bam 
    samtools index ${tgt}.bam
fi
