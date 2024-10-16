process FRAGMENT_LENGTH {
    time = '1d'
    cpus = 1
    memory = '20G'

    tag "Compute fragment lengths on ${meta.sample_id}"

    publishDir "${outdir}/04_reporting/qc/fragment_length/", mode: 'copy'

    input:
    val (meta)
    path (indir)
    path (outdir)

    output:
    tuple path("*")

    script:
    """
    fragment_length.sh ${meta.sample_id} $indir/02_alignment/bowtie2/target/${meta.sample_id}.target.markdup.bam
    """
}
