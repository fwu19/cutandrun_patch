#!/usr/bin/env bash

## Load required modules
# module load BEDTools/2.29.2-GCC-9.3.0
# module load Bowtie2/2.3.4.2-foss-2018b
# module load SAMtools/1.11-GCC-10.2.0
# module load fhR/4.1.2-foss-2021b; export R_LIBS=/fh/fast/_SR/Genomics/user/fwu/R/x86_64-pc-linux-gnu-library/4.1 

bam=$1; shift
filebase=$(basename $bam | sed "s/\.bam$//")
chromSize=$1; shift

## unnormalized bdg
bedtools genomecov -bg -pc -ibam $bam | sort -k1,1 -k2,2n - | egrep -v chrEBV > $filebase.bedgraph
bedGraphToBigWig $filebase.bedgraph $chromSize $filebase.bw
		
## normalized by CPM
bedtools genomecov -scale $( echo -e "1000000/$(samtools view -cf 64 -F 260 $bam)" | bc -l) -bg -pc -ibam $bam | sort -k1,1 -k2,2n - | egrep -v chrEBV >$filebase.CPM.bedgraph	
bedGraphToBigWig $filebase.CPM.bedgraph $chromSize $filebase.CPM.bw

