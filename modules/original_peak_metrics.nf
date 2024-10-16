process ORIGINAL_PEAK_METRICS {
    time = '1d'
    cpus = 1
    memory = '12G'

    tag "Collect QC of original peaks"

    publishDir "${outdir}/04_reporting/qc/", mode: 'copy'

    input:
    path( read_metrics )
    path( peak_list )
    path( rip_list )
    path( outdir )

    output:
    tuple path( "*" )

    script:
    """
    original_peak_metrics.r $read_metrics *.narrowPeak *.broadPeak *.bed *.readsInPeak.csv
    
    """
}
