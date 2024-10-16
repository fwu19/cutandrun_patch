process READS_IN_PEAK {
    time = '1d'
    cpus = 1
    memory = '12G'

    tag "READS_IN_PEAK on ${sample_id}"

    publishDir "${outdir}/04_reporting/qc/reads_in_peak/", mode: 'copy'

    input:
    tuple val( sample_id), path( macs2_peaks ), path( seacr_peaks )
    path( indir )
    path( outdir )

    output:
    tuple val( sample_id ), path( "*" )

    script:
    """
    reads_in_peak.sh ${sample_id} $indir/02_alignment/bowtie2/target/${sample_id}.target.markdup.bam ${macs2_peaks} ${seacr_peaks}
    """
}
