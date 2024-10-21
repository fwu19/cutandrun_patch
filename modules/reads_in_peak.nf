process READS_IN_PEAK {
    label "process_single"
    tag "READS_IN_PEAK on ${sample_id}"

    input:
    tuple val( id), path( peak )

    output:
    tuple val( id ), path( "*.reads_in_peak.csv" )

    script:
    """
    reads_in_peak.sh ${params.indir}/02_alignment/bowtie2/${params.genome}/markdup/${id}.target.markdup.sorted.bam ${peak}
    """
}
