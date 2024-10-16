process CONSENSUS_PEAKS {
    time = '1d'
    cpus = 1
    memory = '12G'

    tag "Generate consensus peaks and collect metrics "

    publishDir "${outdir}/04_reporting/qc/", pattern: '*.{csv,rds}', mode: 'copy'
    publishDir "${outdir}/07_consensus_peaks/", pattern: '*.bed', mode: 'copy'

    input:
    path( replicated_peaks )
    path( outdir )

    output:
    tuple path( "*" )

    script:
    """
    consensus_peaks.r replicated_peaks.rds replicated_peak_metrics.csv
    
    """
}
