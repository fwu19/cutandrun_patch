process CONSENSUS_PEAKS {
    label "process_single"

    tag "Generate consensus peaks and collect metrics "

    input:
    path( replicated_peaks )

    output:
    tuple path( "*" )

    script:
    """
    consensus_peaks.r # replicated_peaks.rds replicated_peak_metrics.csv
    
    """
}
