process PLOT_METRICS {
    label "process_single"

    tag "Make plots of read and peak metrics "

    input:
    path( read_metrics )
    path( original_peaks)
    path( replicated_peaks )
    path( consensus_peaks )

    output:
    tuple path( "*" )

    script:
    """
    plot_metrics.r
    
    """
}
