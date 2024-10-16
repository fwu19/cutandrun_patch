process MACS2 {
    time = '1d'
    cpus = 1
    memory = '12G'
    module = ['MACS2/2.2.6-foss-2019b-Python-3.7.4', 'BEDTools/2.30.0-GCC-10.2.0']

    tag "MACS2 on ${meta.sample_id}"

    publishDir "${outdir}/03_peak_calling/MACS2", mode: 'copy'

    input:
    val (meta)
    path (indir)
    path (outdir)
    val (genome_size)

    output:
    tuple val (meta.sample_id), path ("*.{narrowPeak,broadPeak}"), path( "*.log" )

    script:
    """
    macs2.sh ${meta.sample_id} $indir/02_alignment/bowtie2/target/${meta.sample_id}.target.markdup.bam $indir/02_alignment/bowtie2/target/${meta.control_id}.target.markdup.bam ${meta.peak_type} $genome_size
    
    """
}
