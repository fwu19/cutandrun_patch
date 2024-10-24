#!/usr/bin/env bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fwu@fredhutch.org
#SBATCH --output=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.out
#SBATCH --error=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.err
#SBATCH --nodes=1
#SBATCH --mem=24G
#SBATCH -n 6
#SBATCH --job-name=featCounts

samples="$1"; shift # a file to store sample IDs, one per line
IFS=$'\n\r' read -d '' -r -a samplelist < $samples
IFS=',' read -r -a col <<< "${samplelist[$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))]}"        

sample=${col[0]}
bam=$pubDir/$subdir/02_alignment/bowtie2/target/${sample}.target.markdup.bam
conp="${col[1]}"
prefix="$(echo $conp | sed 's/\.bed$//')"
outdir="$indir/featureCounts/$prefix"

cd "$outdir"
module purge all
module load Subread/2.0.0-GCC-8.3.0  
featureCounts -T 6 -p -C -M --minOverlap 1 -a "$prefix.saf" -F SAF -s 0 -o $sample.fragmentCounts.txt $bam 2>$sample.log