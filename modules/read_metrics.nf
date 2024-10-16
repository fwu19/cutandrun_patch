process READ_METRICS {
    time = '1d'
    cpus = 1
    memory = '12G'

    tag "Collect reads QC"

    publishDir "${outdir}/04_reporting/qc/", mode: 'copy'

    input:
    path (sample_sheet)
    path (indir)
    path (outdir)
    path (frag_lens)

    output:
    tuple path( "*" )

    script:
    """
    read_metrics.r $sample_sheet $indir/04_reporting/meta_table_ctrl.csv 
    
    """
}
