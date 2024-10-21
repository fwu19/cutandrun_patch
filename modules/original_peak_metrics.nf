process ORIGINAL_PEAK_METRICS {
    label "process_single"

    tag "Collect QC of original peaks"

    input:
    path( read_metrics )
    path( peak_list, stageAs: "peaks/*" )
    path( rip_list, stageAs: "rip/*" )

    output:
    tuple path( "*.{rds,csv}" )

    script:
    """
    original_peak_metrics.r # $read_metrics *.narrowPeak *.broadPeak *.bed *.reads_in_peak.csv
    
    """
}
