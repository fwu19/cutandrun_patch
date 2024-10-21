#!/usr/bin/env bash

shdir=$(dirname $0)

slurm="--mail-type=END,FAIL --mail-user=fwu@fredhutch.org --output=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A_%a.%x.out --error=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A_%a.%x.err --parsable"

json=$( realpath $1 ); shift
nfDir=$( cat $json | jq '.indir' | sed "s/\"//g; s/\'//g" )
pubDir=$( cat $json | jq '.outdir' | sed "s/\"//g; s/\'//g" )

samples=$1; shift 
nstart=2 # given the header line
nend=$(egrep -c "^" $samples )

log=cutandrun.log
echo -e "###### $( date ) ######" | tee -a $log
echo -e "Run Cut&Run workflow on $samples" | tee -a $log

## get short fragments for TFs
tf=$(sbatch $slurm --array=${nstart}-${nend} $shdir/filter_short_fragments.sh $json $samples shortFragments)
echo -e "filter_short_fragments: $tf" | tee -a $log

macs2=$(sbatch $slurm --dependency=afterok:$tf --array=${nstart}-${nend} $shdir/macs2.sh $json $samples shortFragments )
echo -e "MACS2_IgG_markdup: $macs2" | tee -a $log

gcov=$(sbatch $slurm --dependency=afterok:$tf --array=${nstart}-${nend} $shdir/genome_coverage.sh $json $samples shortFragments )
echo -e "Genome Coverage: $gcov" | tee -a $log

seacr=$(sbatch $slurm --dependency=afterok:$gcov --array=${nstart}-${nend} $shdir/seacr.sh $json $samples shortFragments )
echo -e "SEACR_IgG_dedup: $seacr" | tee -a $log

rip=$(sbatch $slurm --dependency=afterok:$macs2:$seacr --array=${nstart}-${nend} $shdir/reads_in_peak.sh $json $samples shortFragments )
echo -e "Reads in Peak: $rip" | tee -a $log

binr=$(sbatch $slurm --dependency=afterok:$tf --array=${nstart}-${nend} $shdir/bin_reads.sh $json $samples shortFragments )
echo -e "Bin Reads: $binr" | tee -a $log

readms=$(sbatch $slurm --dependency=afterok:$tf $shdir/read_metrics_shortFragments.sh $json $samples shortFragments )
echo -e "Read metrics: $readms" | tee -a $log
# missing a plot to show short fragment counts?

peakms=$(sbatch $slurm --dependency=afterok:$macs2:$seacr $shdir/collect_peak_metrics.sh $json shortFragments )
echo -e "Peak metrics: $peakms" | tee -a $log

peakplot=$(sbatch $slurm --dependency=afterok:$peakms $shdir/plot_peak_metrics.sh $json shortFragments )
echo -e "Plot peak metrics: $peakplot" | tee -a $log

ann=$(sbatch $slurm --dependency=afterok:$peakms $shdir/annotate_conp.sh $json shortFragments ) # need to create a temp dir
echo -e "Annotate consensus peaks: $ann" | tee -a $log

fc=$(sbatch $slurm --dependency=afterok:$peakms $shdir/count_reads_in_conp.sh params.json shortFragments samples.DP.csv) # need to create sample sheet separately
echo -e "Count reads in consensus peaks: $fc" | tee -a $log

dp=$(sbatch $slurm --dependency=afterok:$fc scripts/shortFragments/sbatch_R.sh scripts/shortFragments/call_differential_peaks.r samples.targetOnly.csv samples.DP.csv comparisons.csv Analysis/shortFragments)
echo -e "Call differential peaks: $dp" | tee -a $log

