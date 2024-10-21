process MERGE_REPLICATES {
    labe "process_single"
//    module = ['Bowtie2/2.4.2-GCC-10.2.0', 'BEDTools/2.30.0-GCC-10.2.0', 'SAMtools/1.11-GCC-10.2.0']


    tag "Compute genome coverage on ${meta.sample_id}"

    publishDir "${outdir}/05_genome_coverage/", mode: 'copy'

    input:
    val (meta)
    path (indir)
    val (suffix) // target.markdup or target.dedup
    path (outdir)
    path (chromSize)

    output:
    tuple val (meta), path ( ${outbase}.{bam,bai} ), path ( ${outbase}.{markdup,dedup}.bedgraph ), path ( ${outbase}.CPM.{bedgraph,bw}" ) 
    
    script:
    """
    merge_replicates.sh $indir/02_alignment/bowtie2/target/${meta.sample_id}.${suffix}.bam $chromSize
    
    """
}
