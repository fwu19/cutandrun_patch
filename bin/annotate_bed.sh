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

bed=$1; shift # a file to store sample IDs, one per line
genome=$1; shift
gtf4homer=$(realpath $1); shift

if [[ -f $bed ]]; then
    dir=$( dirname $bed )
    bed=$( basename $bed )
    prefix=$(echo $bed | sed "s/.bed$//")
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


peak <- read.delim(bed, header = F) %>% 
    mutate(PeakID = paste0(V1, ":", V2+1, "-", V3)) %>% 
    left_join(ann, by = 'PeakID')

write.table(peak[order(peak[,1], peak[,2]),], paste(prefix, 'withAnnotation.txt', sep = '.'), sep = '\t', quote = F, row.names = F) 
code

fi

