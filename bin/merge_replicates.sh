#!/usr/bin/env bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fwu@fredhutch.org
#SBATCH --output=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.out
#SBATCH --error=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.err
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH -n 6

## Load required modules
module purge all
module load BEDTools/2.29.2-GCC-9.3.0
module load SAMtools/1.11-GCC-10.2.0


## set up output directory and file base
outprefix=$1; shift

## merge bam files
samtools merge -@ 6 $outprefix.bam $@
samtools index $outprefix.bam


## unnormalized bdg
bedtools genomecov -bg -pc -ibam $outprefix.bam | sort -k1,1 -k2,2n - | egrep -v chrEBV > $outprefix.bedgraph
		
## normalized by CPM (failed?)
#bedtools genomecov -scale $( echo -e "1000000/$(samtools view -cf 64 -F 256 $outprefix.bam)" | bc -l) -bg -pc -ibam $outprefix.bam | sort -k1,1 -k2,2n - | egrep -v chrEBV >$outprefix.CPM.bedgraph
#~/src/UCSCtools/bedGraphToBigWig $outprefix.CPM.bedgraph $chromSize $outprefix.CPM.bw
#rm $outprefix.CPM.bedgraph

						
