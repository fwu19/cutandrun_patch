process SEACR {
    time = '1d'
    cpus = 1
    memory = '12G'

    tag "SEACR on $target_id"

    publishDir "${outdir}/03_peak_calling/SEACR", mode: 'copy'

    input:
    tuple val (target_id), val (control_id), path (target_bdg), path (control_bdg)
    path (outdir)

    output:
    tuple val (target_id), path( "*.bed" ), path( "*.log" )

    script:
    """
    seacr.sh ${target_id} ${target_bdg} ${control_bdg} // target_id is output prefix
    """
}
