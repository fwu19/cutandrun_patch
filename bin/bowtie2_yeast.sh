#!/usr/bin/env bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fwu@fredhutch.org
#SBATCH --output=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.out
#SBATCH --error=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.err
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH -n 6

set -uo pipefail

json=params.json
localConfig=nextflow.config

samples=$( cat $json | jq '.input' | sed "s/\"//g; s/\'//g" )
samples=${1:-$samples}; shift

pubDir=$( cat $json | jq '.outdir' | sed "s/\"//g; s/\'//g" )
pubDir=${1:-$pubDir}; shift
[[ -d $pubDir ]] || mkdir -p $pubDir
pubDir=$(realpath $pubDir)

bt2RefYeast=$(realpath $( cat $json | jq '.bt2RefYeast' | sed "s/\"//g; s/\'//g" ))

if [[ -f $samples ]]; then
	IFS=$'\n\r' read -d '' -r -a samplelist < $samples
else
	IFS=' ' read -r -a samplelist <<< "$samples"
fi	


## Load required modules
module purge all
module load Bowtie2/2.5.4-GCC-13.2.0
module load SAMtools/1.19.2-GCC-13.2.0


## cd working directory
workDir=$pubDir/02_alignment/bowtie2/yeast
[[ -d $workDir ]] || mkdir -p $workDir
cd $workDir
[[ -d log ]] || mkdir -p log

## Process one sample
IFS=',' read -r -a arr <<< "${samplelist[$(( ${SLURM_ARRAY_TASK_ID:-2} - 1 ))]}"        
id=${arr[0]}_R${arr[1]}
fq1=$pubDir/01_prealign/trimgalore/${id}_1.trimmed.fastq.gz
fq2=$pubDir/01_prealign/trimgalore/${id}_2.trimmed.fastq.gz

bowtie2 --local --very-sensitive --no-overlap --no-dovetail --no-unal --no-mixed --no-discordant -q -I 10 -X 700 --phred33 --threads 6 -x $bt2RefYeast -1 $fq1 -2 $fq2 2> log/$id.yeast.bowtie2.log | samtools view -@ 6 -bh - > $id.yeast.bam

						
