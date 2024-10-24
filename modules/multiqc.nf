process MULTIQC {
    module = ['MultiQC/1.9-foss-2019b-Python-3.7.4']

    label "process_high"

    tag "MultiQC on all samples"

    output:
    path ('multiqc_data/')
    path ('multiqc_report.html')

    script:
    """
    multiqc -o . -f ${params.indir}/01_prealign/ ${params.indir}/02_alignment/bowtie2/*/log/ ${params.indir}/02_alignment/bowtie2/*/markdup/ ${params.indir}/03_peak_calling/*/

    """
}