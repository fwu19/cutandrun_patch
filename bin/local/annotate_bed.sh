#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --mem=12G
#SBATCH -n 1


###################################################################################
#### HOMER analysis
###################################################################################
## Load required modules
module purge all
module load fhR/4.1.2-foss-2021b; export R_LIBS=/fh/fast/_SR/Genomics/user/fwu/R/x86_64-pc-linux-gnu-library/4.1 
export PATH=/home/fwu/src/homer/bin/:$PATH

bed=$(realpath $1); shift # a file to store sample IDs, one per line
genome=$1; shift
gtf4homer=$(realpath $1); shift

if [[ -f $bed ]]; then
    dir=$( echo $bed | sed "s/.bed$//" )
    prefix=$(basename $dir)
    [[ -d $dir ]] || mkdir -p $dir
    cd $dir
    cat $bed | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$4,".","+"}' >$prefix.homer.bed
    annotatePeaks.pl $prefix.homer.bed $genome -gtf $gtf4homer >$prefix.annotation.txt 2>$prefix.annotation.log
    rm $prefix.homer.bed
    
R --vanilla "--args $bed $prefix" <<code
options(stringsAsFactor=F); 
library(dplyr)

args <- commandArgs(T); bed = args[1]; prefix = args[2]; 

ann <- read.delim(paste(prefix, 'annotation.txt', sep = '.'))[c(1,8:13,16,19)]
colnames(ann)[1] <- 'PeakID'

if (grepl('stringent.bed$', bed)){
    peak <- read.delim(bed, header = F, col.names = c('chrom', 'start', 'end', 'Total.Signal', 'Max.Signal', 'Max.Signal.Region')) 
}else if(grepl('\\.bed$', bed)){
    peak <- read.delim(bed, header = F, 
    col.names = c('chrom', 'start', 'end', 'PeakID', 'Sample.Groups', 'Peak.Count')) 
}else if (grepl('narrowPeak$', bed)){
    peak <- read.delim(bed, header = F, 
    col.names = c('chrom', 'start', 'end', 'name', 'score', 'strand', 'fold.change', 'p.value', 'q.value', 'summit')) 
}else if (grepl('broadPeak$', bed)){
    peak <- read.delim(bed, header = F, 
    col.names = c('chrom', 'start', 'end', 'name', 'score', 'strand', 'fold.change', 'p.value', 'q.value')) 
}else{
    peak <- read.delim(bed, header =F)
    colnames(peak)[1:3] <- c('chrom', 'start', 'end')
}

if(length(args) > 2){
	colv <- eval(parse(text=args[3]))
	colnames(peak)[1:length(colv)] <- colv
}
peak <- left_join(peak, ann, by = c('PeakID' = 'PeakID'))

write.table(peak[order(peak[,1], peak[,2]),], paste(prefix, 'withAnnotation.txt', sep = '.'), sep = '\t', quote = F, row.names = F) 
code

fi

