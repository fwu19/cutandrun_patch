// #!/usr/bin/env nextflow 

//include { CHECK_INPUT  } from './modules/check_input.nf'
include { MULTIQC  } from './modules/multiqc.nf'
include { FRAGMENT_LENGTH  } from './modules/fragment_length.nf'
include { READS_IN_PEAK  } from './modules/reads_in_peak.nf'
include { READ_METRICS  } from './modules/read_metrics.nf'
include { ORIGINAL_PEAK_METRICS  } from './modules/original_peak_metrics.nf'
include { REPLICATED_PEAKS } from './modules/replicated_peaks.nf'
include { CONSENSUS_PEAKS } from './modules/consensus_peaks.nf'
include { PLOT_METRICS } from './modules/plot_metrics.nf'
include { OUTPUT_PARAMS  } from './modules/output_params.nf'


Channel
        .fromPath (params.input, checkIfExists: true)
        .splitCsv (header: true)
        .set { ch_samples }

ch_samples
        .filter { it.is_control.toBoolean() == true }
        .set { ch_control }
        
ch_samples
        .filter { it.is_control.toBoolean() == false }
        .set { ch_target }


Channel
        .fromPath (params.peaks, checkIfExists: true)
        .splitCsv (header: false)
        .set { ch_peaks }
        // [sample_id, peak_file]

workflow {
         //CHECK_INPUT(params.input)

        /* Run MultiQC */
        MULTIQC()

        /* Compute fragment length */
        FRAGMENT_LENGTH(ch_samples)
        .flatten()
        .collect()
        .set{ frag_lens }

         /* Compute reads in peak */
        READS_IN_PEAK( ch_peaks )
        .map{ it -> it[1] }
        .collect()
        .set{ rip_list }


        /* Collect reads metrics */
        read_metrics = READ_METRICS( frag_lens )


        /* Collect metrics for original peaks */
        ch_peaks
                .map{ it -> it[1]  }
                .collect()
                .set{ peak_list }

        ORIGINAL_PEAK_METRICS( read_metrics, peak_list, rip_list)
                .set{ original_peaks }


        /* Generate replicated peaks and collect metrics  */
        replicated_peaks = REPLICATED_PEAKS( original_peaks, params.minReplicates )


        /* Generate consensus peaks and collect metrics  */
        consensus_peaks = CONSENSUS_PEAKS( replicated_peaks )


        /* Make plots for report */
        plots = PLOT_METRICS( read_metrics, original_peaks, replicated_peaks, consensus_peaks )


        /* Report params used for the workflow */
        OUTPUT_PARAMS( params.outdir )
}
