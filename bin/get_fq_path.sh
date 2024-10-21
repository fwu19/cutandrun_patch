#!/usr/bin/env bash

echo "seq.id,fastq_1,fastq_2" >fq.csv
for fqdir in $@; do
	paste -d "," <(find $fqdir/Unaligned/ | egrep "_S[0-9]+_(L[0-9]+_)?R1_.*fastq.gz" | sort | sed "s/.*\///g; s/_S[0-9]\+.*_R1_.*//g" ) <(find $fqdir/Unaligned/ | egrep "_S[0-9]+_(L[0-9]+_)?R1_.*fastq.gz" | sort ) <(find $fqdir/Unaligned/ | egrep "_S[0-9]+_(L[0-9]+_)?R2_.*fastq.gz" |sort )
done >>fq.csv

