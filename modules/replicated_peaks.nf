process REPLICATED_PEAKS {
    time = '1d'
    cpus = 1
    memory = '12G'

    tag "Generate replicated peaks and collect metrics"

    publishDir "${outdir}/04_reporting/qc/", pattern: '*.{csv,rds}', mode: 'copy'
    publishDir "${outdir}/06_replicated_peaks/", pattern: '*.bed', mode: 'copy'

    input:
    path( samplesheet )
    path( original_peaks )
    val( min_reps )
    path( outdir )

    output:
    tuple path( "*" )

    script:
    """
    replicated_peaks.r $samplesheet original_peaks.rds original_peak_metrics.csv $min_reps
    
    """
}
