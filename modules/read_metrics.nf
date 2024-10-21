process READ_METRICS {
    label "process_single"
    tag "Collect reads QC"

    input:
    path (frag_lens)

    output:
    tuple path( "*" )

    script:
    """
    read_metrics.r ${params.input} ${params.indir} ${params.genome} ${params.spikein_genome} ${params.spikein_genome_2}
    
    """
}
