process PLOT_METRICS {
    time = '1d'
    cpus = 1
    memory = '12G'

    tag "Make plots of read and peak metrics "

    publishDir "${outdir}/04_reporting/qc/", pattern: '*.{csv,rds}', mode: 'copy'

    input:
    path( read_metrics )
    path( original_peaks)
    path( replicated_peaks )
    path( consensus_peaks )
    path( outdir )

    output:
    tuple path( "*" )

    script:
    """
    plot_metrics.r
    
    """
}
