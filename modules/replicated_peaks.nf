process REPLICATED_PEAKS {
    label "process_single"

    tag "Generate replicated peaks and collect metrics"

    input:
    path( original_peaks )
    val( min_reps )
    

    output:
    tuple path( "*" )

    script:
    """
    replicated_peaks.r ${params.input} $min_reps # original_peaks.rds original_peak_metrics.csv 
    
    """
}
