#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --mem=32G
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
        echo -e "USAGE: bash `basename $0` <params> <samples.csv> <path/to/input/bam/dir>"
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
genomeSize=$( cat $json | jq '.genomeSize' | sed "s/\"//g; s/\'//g" )

samples=$( realpath $1 ); shift
if [[ -f $samples ]]; then
	IFS=$'\n\r' read -d '' -r -a samplelist < $samples
else
	IFS=' ' read -r -a samplelist <<< "$samples"
fi	
IFS=',' read -r -a arr <<< "${samplelist[$(( ${SLURM_ARRAY_TASK_ID:-2} - 1 ))]}"        

subdir=$1; shift # under $pubDir
binLen=${1:-500}; shift

## Process one sample
prefix=${arr[12]}
indir=$pubDir/$subdir/02_alignment/bowtie2/target
outdir=$pubDir/$subdir/04_reporting/qc/bin_reads
[[ -d $outdir ]] || mkdir -p $outdir
cd $outdir

if [[ -f $indir/${prefix}.target.dedup.bam  ]]; then     
	tgt=${prefix}.target.dedup        
else
	tgt=${prefix}.target.markdup
fi
      			
## Filter and keep the mapped read pairs
## Convert into bed file format
## Only extract the fragment related columns
## We use the mid point of each fragment to infer which 500bp bins does this fragment belong to.
samtools view -bf 0x02 $indir/$tgt.bam | samtools sort -@ 6 -n - | bedtools bamtobed -i - -bedpe | cut -f 1,2,6 | sort -k1,1 -k2,2n -k3,3n | awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  >$tgt.fragmentsCount.bin$binLen.bed


