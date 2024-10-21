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
module load BEDTools/2.29.2-GCC-9.3.0
module load Bowtie2/2.3.4.2-foss-2018b
module load SAMtools/1.11-GCC-10.2.0
module load Kent_tools/20201201-linux.x86_64
module load fhR/4.1.2-foss-2021b; export R_LIBS=/fh/fast/_SR/Genomics/user/fwu/R/x86_64-pc-linux-gnu-library/4.1 


## Source params and check whether all required variables are defined
json=$( realpath $1 ); shift
nfDir=$( cat $json | jq '.indir' | sed "s/\"//g; s/\'//g" )
pubDir=$( cat $json | jq '.outdir' | sed "s/\"//g; s/\'//g" )
chromSize=$( cat $json | jq '.chromSize' | sed "s/\"//g; s/\'//g" )

samples=$( realpath $1 ); shift
if [[ -f $samples ]]; then
	IFS=$'\n\r' read -d '' -r -a samplelist < $samples
else
	IFS=' ' read -r -a samplelist <<< "$samples"
fi	
IFS=',' read -r -a arr <<< "${samplelist[$(( ${SLURM_ARRAY_TASK_ID:-2} - 1 ))]}"        

subdir=$1; shift # under $pubDir

## Process one sample
prefix=${arr[12]}
indir=$pubDir/$subdir/02_alignment/bowtie2/target
outdir=$pubDir/$subdir/05_genome_coverage
[[ -d $outdir ]] || mkdir -p $outdir
cd $outdir

tgt=${prefix}.target.markdup # markdup.bam
if [[ -f $indir/$tgt.bam  ]]; then     
	## unnormalized bdg
	bedtools genomecov -bg -pc -ibam $indir/$tgt.bam | sort -k1,1 -k2,2n - | egrep -v chrEBV > $tgt.bedgraph
	bedGraphToBigWig $tgt.bedgraph $chromSize $tgt.bw
		
	## normalized by CPM
	bedtools genomecov -scale $( echo -e "1000000/$(samtools view -cf 64 -F 260 $indir/$tgt.bam)" | bc -l) -bg -pc -ibam $indir/$tgt.bam | sort -k1,1 -k2,2n - | egrep -v chrEBV >$tgt.CPM.bedgraph
	bedGraphToBigWig $tgt.CPM.bedgraph $chromSize $tgt.CPM.bw
	R --vanilla "--args $tgt.CPM.bedgraph" <<code
	options(stringsAsFactor=F)
	options(scipen=999)
	library(GenomicRanges)
	library(dplyr)
	args <- commandArgs(T)
	read.delim(args[1], header = F, col.names = c('chrom', 'start', 'end', 'coverage')) %>% 
	makeGRangesFromDataFrame(ignore.strand = T, keep.extra.columns = T, seqnames.field = 'chrom', start.field = 'start', end.field = 'end', starts.in.df.are.0based = T) %>% 
	saveRDS(gsub('.bedgraph$', '.rds', args[1]))
code
	
	rm $tgt.CPM.bedgraph

fi

tgt=${prefix}.target.dedup  # dedup.bam
if [[ -f $indir/$tgt.bam  ]]; then     	     
	## unnormalized bdg
	bedtools genomecov -bg -pc -ibam $indir/$tgt.bam | sort -k1,1 -k2,2n - | egrep -v chrEBV > $tgt.bedgraph
	bedGraphToBigWig $tgt.bedgraph $chromSize $tgt.bw
		
	## normalized by CPM
	bedtools genomecov -scale $( echo -e "1000000/$(samtools view -cf 64 -F 260 $indir/$tgt.bam)" | bc -l) -bg -pc -ibam $indir/$tgt.bam | sort -k1,1 -k2,2n - | egrep -v chrEBV >$tgt.CPM.bedgraph
	bedGraphToBigWig $tgt.CPM.bedgraph $chromSize $tgt.CPM.bw
	R --vanilla "--args $tgt.CPM.bedgraph" <<code
	options(stringsAsFactor=F)
	options(scipen=999)
	library(GenomicRanges)
	library(dplyr)
	args <- commandArgs(T)
	read.delim(args[1], header = F, col.names = c('chrom', 'start', 'end', 'coverage')) %>% 
	makeGRangesFromDataFrame(ignore.strand = T, keep.extra.columns = T, seqnames.field = 'chrom', start.field = 'start', end.field = 'end', starts.in.df.are.0based = T) %>% 
	saveRDS(gsub('.bedgraph$', '.rds', args[1]))
code
	
	rm $tgt.CPM.bedgraph
	
fi

