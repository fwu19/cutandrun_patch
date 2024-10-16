#!/usr/bin/env bash

## Load required modules
# module load SAMtools/1.11-GCC-10.2.0
# module load BEDTools/2.29.2-GCC-9.3.0

prefix=$1; shift
bam=$1; shift; [[ -f $bam ]] || exit 1

for peak in "$@"; do
	echo -e $(basename $peak ),$( bedtools intersect -a $bam -b $peak -u | samtools view -F 256 -f 64 -c - )
done >$prefix.readsInPeak.csv # count only R1 reads
