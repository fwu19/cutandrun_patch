#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --mem=24G
#SBATCH -n 4
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
        echo -e "USAGE: bash `basename $0` <params> <samples.csv> <path/to/output/peak/dir>"
        exit 1;
fi

## Load required modules
module purge all
module load fhR/4.1.2-foss-2021b
export R_LIBS=/fh/fast/_SR/Genomics/user/fwu/R/x86_64-pc-linux-gnu-library/4.1 # 
module load BEDTools/2.29.2-GCC-9.3.0
sh=/fh/fast/_SR/Genomics/user/fwu/src/SEACR-1.3/SEACR_1.3.sh

## Source params and check whether all required variables are defined
json=$( realpath $1 ); shift
nfDir=$( cat $json | jq '.indir' | sed "s/\"//g; s/\'//g" )
pubDir=$( cat $json | jq '.outdir' | sed "s/\"//g; s/\'//g" )

samples=$1; shift # a file to store sample IDs, one per line
if [[ -f $samples ]]; then
	IFS=$'\n\r' read -d '' -r -a samplelist < $samples
else
	IFS=' ' read -r -a samplelist <<< "$samples"
fi	
IFS=',' read -r -a arr <<< "${samplelist[$(( ${SLURM_ARRAY_TASK_ID:-2} - 1 ))]}"        

subdir=$1; shift # under $pubDir

chbdg=$(realpath $pubDir)/$subdir/05_genome_coverage/${arr[12]}.target.markdup.bedgraph
inbdg=$(realpath $pubDir)/mergedIgG/05_genome_coverage/${arr[13]}.target.dedup.bedgraph
type=${arr[6]}
prefix=${arr[12]}

if [[ $type != "NA" && $type != "" && -f $chbdg ]]; then # IgG with type of ""
    outdir=$pubDir/$subdir/03_peak_calling/SEACR/$prefix
    [[ -d $outdir ]] || mkdir -p $outdir
    cd $outdir
    log=$prefix.seacr.log

    ## SEACR, unormalized & IgG & norm ############################
	bash $sh $chbdg 0.01 non stringent $prefix.top01 2>&1 | tee -a $log
	
    if [[ -f $chbdg && -f $inbdg ]]; then
	    bash $sh $chbdg $inbdg norm stringent $prefix.IgG.norm 2>&1 | tee -a $log
	    bedtools intersect -a $prefix.IgG.norm.stringent.bed -b $prefix.top01.stringent.bed -u >$prefix.IgG.norm.stringent.filtered.bed
    fi

fi

