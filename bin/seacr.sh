#!/usr/bin/env bash

## Load required modules
# module purge all
# module load fhR/4.1.2-foss-2021b
# export R_LIBS=/fh/fast/_SR/Genomics/user/fwu/R/x86_64-pc-linux-gnu-library/4.1 # 
# module load BEDTools/2.29.2-GCC-9.3.0

prefix=$1; shift     
chBdg=$1; shift
inBdg=$1; shift

log=$prefix.log

if [[ -f $chBdg ]]; then
    SEACR_1.3.sh $chBdg 0.01 non stringent $prefix.top01 2>&1 >> $log
	
fi

if [[ -f $chBdg && -f $inBdg ]]; then
    SEACR_1.3.sh $chBdg $inBdg norm stringent $prefix.IgG.norm 2>&1 >> $log
	bedtools intersect -a $prefix.IgG.norm.stringent.bed -b $prefix.top01.stringent.bed -u >$prefix.IgG.norm.stringent.filtered.bed
fi

