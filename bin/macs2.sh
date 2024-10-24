#!/usr/bin/env bash

## Load required modules
# module load MACS2/2.2.6-foss-2019b-Python-3.7.4
# module load BEDTools/2.29.2-GCC-9.3.0

prefix=$1; shift     
chBam=$1; shift
inBam=$1; shift
type=$1; shift
genomeSize=$1; shift

###### call narrow peaks using keep-dup all ######
if [[ $type =~ 'narrow' || $type =~ 'both' ]]; then
	log=$prefix.narrow_peaks.log

	## call peaks with no control
	echo -e "############\n####Call narrow peaks without control\n" | tee $log
	macs2 callpeak --name $prefix.noControl.narrow  --treatment $chBam --outdir ./  --format BAMPE --gsize $genomeSize --keep-dup all -q 0.05 2>> $log
	
	if [[ -f $inBam ]]; then
		## call peaks with respect to IgG ######
		echo -e "############\n####Call narrow peaks with control\n" | tee -a $log
		macs2 callpeak --name $prefix.toIgG.narrow  --treatment $chBam --control $inBam --outdir ./  --format BAMPE --gsize $genomeSize --keep-dup all -q 0.05 2>> $log

		echo -e "############\n####Intersect peaks with and without control\n" | tee -a $log
		bedtools intersect -a $prefix.toIgG.narrow_peaks.narrowPeak -b $prefix.noControl.narrow_peaks.narrowPeak -u | awk -F "\t" '{if ($7 > 2 && $9 > 2) {print} }' >$prefix.toIgG.narrow_peaks.filtered.narrowPeak 2>> $log

	fi

fi


###### call broad peaks using keep-dup all ######
if [[ $type =~ 'broad' || $type =~ 'both' ]]; then
	log=$prefix.broad_peaks.log

	echo -e "############\n####Call broad peaks without control\n" | tee $log
	macs2 callpeak --name $prefix.noControl.broad  --treatment $chBam --outdir ./  --format BAMPE --gsize $genomeSize --keep-dup all --broad --broad-cutoff 0.1 -q 0.05 --nomodel --nolambda 2>> $log

	if [[ -f $inBam ]]; then

		echo -e "############\n####Call broad peaks with control\n" | tee -a $log
		macs2 callpeak --name $prefix.toIgG.broad  --treatment $chBam --control $inBam --outdir ./  --format BAMPE --gsize $genomeSize --keep-dup all --broad --broad-cutoff 0.1 -q 0.05 2>> $log

		echo -e "############\n####Intersect peaks with and without control\n" >>$log
		bedtools intersect -a $prefix.noControl.broad_peaks.broadPeak -b $prefix.toIgG.broad_peaks.broadPeak -u | awk -F "\t" '{if ($7 > 2 && $9 > 2) {print} }'  >$prefix.noControl.broad_peaks.filtered.broadPeak 2>> $log
	
	fi

fi