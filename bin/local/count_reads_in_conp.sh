#!/usr/bin/env bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fwu@fredhutch.org
#SBATCH --output=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.out
#SBATCH --error=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.err
#SBATCH --nodes=1
#SBATCH --mem=24G
#SBATCH -n 6

# check if script is started via SLURM or bash
# if with SLURM: there variable '$SLURM_JOB_ID' will exist
# `if [ -n $SLURM_JOB_ID ]` checks if $SLURM_JOB_ID is not an empty string
if [[ -n $SLURM_JOB_ID ]];  then
    # check the original location through scontrol and $SLURM_JOB_ID
    SCRIPT_PATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
else
    # otherwise: started with bash. Get the real location.
    SCRIPT_PATH=$(realpath $0)
fi
shDir=$(dirname $SCRIPT_PATH)

## Source params and check whether all required variables are defined
json=$( realpath $1 ); shift
pubDir=$( cat $json | jq '.outdir' | sed "s/\"//g; s/\'//g" )
subdir=$1; shift # under $pubDir
samples=$(realpath $1); shift

## Source params and check whether all required variables are defined
indir=$pubDir/$subdir/07_consensus_peaks/

for conp in $(ls $indir/*.bed); do
  outdir=$indir/featureCounts/$(basename $conp | sed "s/.bed$//")
  [[ -d $outdir ]] || mkdir -p $outdir

    cat $conp | awk 'BEGIN {FS=OFS="\t"} { $2=$2+1; print $1":"$2"-"$3, $1, $2, $3, ".", $4, $5 }' >$outdir/$(basename $conp | sed "s/.bed$/.saf/" )

done

slurm="--parsable --output=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.out --error=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.err --nodes=1 --mem=24G -n 6"
tail -n +2 $samples | while IFS='' read -r line || [[ -n "$line" ]]; do
  IFS=',' read -r -a col <<< "$line"
  sample=${col[0]}
  bam=$pubDir/$subdir/02_alignment/bowtie2/target/${sample}.target.markdup.bam
  conp=${col[1]}
  prefix=$(echo $conp | sed "s/\.bed$//")
  outdir=$indir/featureCounts/$prefix

  echo -e $(sbatch $slurm --wrap="cd $outdir; module purge all; module load Subread/2.0.0-GCC-8.3.0;  featureCounts -T 6 -p -C -M --minOverlap 1 -a $prefix.saf -F SAF -s 0 -o $sample.fragmentCounts.txt $bam 2>$sample.log"),$line,count_reads_in_conp.sh
done


