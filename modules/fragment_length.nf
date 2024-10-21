process FRAGMENT_LENGTH {
    label "process_single"

    tag "Compute fragment lengths on ${meta.sample_id}"

    input:
    val (meta)

    output:
    tuple path("*")

    script:
    """
    fragment_length.sh ${meta.id} ${params.indir}/02_alignment/bowtie2/${params.genome}/markdup/${meta.id}.target.markdup.sorted.bam
    """
}
