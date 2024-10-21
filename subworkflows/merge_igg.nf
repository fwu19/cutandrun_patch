// #!/usr/bin/env nextflow 

//include { CHECK_INPUT  } from './modules/check_input.nf'
include { GENOME_COVERAGE as GENOMECOV_MARKDUP  } from './modules/genome_coverage.nf'
include { GENOME_COVERAGE as GENOMECOV_DEDUP  } from './modules/genome_coverage.nf'
include { MACS2  } from './modules/macs2.nf'
include { SEACR  } from './modules/seacr.nf'
include { FRAGMENT_LENGTH  } from './modules/fragment_length.nf'
include { READS_IN_PEAK  } from './modules/reads_in_peak.nf'
include { READ_METRICS  } from './modules/read_metrics.nf'
include { ORIGINAL_PEAK_METRICS  } from './modules/original_peak_metrics.nf'
include { REPLICATED_PEAKS } from './modules/replicated_peaks.nf'
include { CONSENSUS_PEAKS } from './modules/consensus_peaks.nf'
include { PLOT_METRICS } from './modules/plot_metrics.nf'
include { OUTPUT_PARAMS  } from './modules/output_params.nf'
include { TEST  } from './modules/test.nf'

Channel
        .fromPath (params.input, checkIfExists: true)
        .splitCsv (header: true)
        .set { meta }

meta
        .filter { it.is_control.toInteger() == 1 }
        .set { control_meta }
        
meta
        .filter { it.is_control.toInteger() == 0 }
        .set { target_meta }


workflow {
    //CHECK_INPUT(params.input)
    
    /* Compute genome coverage profile for target samples */
    genomecov_markdup = GENOMECOV_MARKDUP(meta, params.indir, "target.markdup", params.outdir, params.chromSize)
    genomecov_markdup
                    .filter { it[0].is_control.toInteger() == 0 }
                    .map { it -> [ it[0].control_id, it[0].sample_id, it[1] ]}
                    .set { target_markdup_bdg }

    /* Compute genome coverage profile for IgG samples */
    genomecov_dedup = GENOMECOV_DEDUP(control_meta, params.indir, "target.dedup", params.outdir, params.chromSize)
    genomecov_dedup
                    .map { it -> [ it[0].sample_id, it[1] ]}
                    .set { control_dedup_bdg }
    
    
    /* Merge IgG replicates */
    control_meta
                .map { it -> [ it.sample_group, it.group ] }
                .groupTuple ()
                .view()


    /* Pair genome coverage profile between target samples and IgG samples */
    control_dedup_bdg
            .cross(target_markdup_bdg)
            .map { it -> [ it[1][1], it[1][0], it[1][2], it[0][1] ] } // sample_id, control_id, target_bdg, control_bdg
            .set { paired_bdg }
    
    //paired_bdg.view { "PAIRED: $it" }
    
    
    /* Call SEACR peaks */
    seacr = SEACR( paired_bdg, params.outdir) 

    /* Call MACS2 peaks */
    macs2 = MACS2(target_meta, params.indir, params.outdir, params.genomeSize)

    /* Compute fragment length */
    FRAGMENT_LENGTH(meta, params.indir, params.outdir)
        .flatten()
        .collect()
        .set{ frag_lens }
    
    /* Compute reads in peak */
    macs2
        .cross( seacr )
        .map{ it -> [ it[0][0], it[0][1], it[1][1] ] }
        .set{ peaks_ch }
        
    //peaks_ch.view()
    
    READS_IN_PEAK( peaks_ch, params.indir, params.outdir )
        .map{ it -> it[1] }
        .flatten()
        .collect()
        .set{ rip_list }
    
    
    /* Collect reads metrics */
    read_metrics = READ_METRICS( params.input, params.indir, params.outdir, frag_lens )
    
    
    /* Collect metrics for original peaks */
    peaks_ch
            .map{ it -> [ it[1], it[2] ] }
            .flatten()
            .collect()
            .set{ peak_list } // paths to all peak files

    ORIGINAL_PEAK_METRICS( read_metrics, peak_list, rip_list, params.outdir)
            .collect()
            .set{ original_peaks } // orginal_peaks.rds, original_peak_metrics.csv, original_peak_widths.rds


    /* Generate replicated peaks and collect metrics */
    replicated_peaks = REPLICATED_PEAKS( original_peaks, params.minReplicates, params.outdir ) // replicated_peaks.rds, replicated_peaks_metrics.csv


    /* Generate consensus peaks and collect metrics  */
    consensus_peaks = CONSENSUS_PEAKS( replicated_peaks, params.outdir )


    /* Make plots for report */
    plots = PLOT_METRICS( read_metrics, original_peaks, replicated_peaks, consensus_peaks, params.outdir )
    
    
    /* Report params used for the workflow */
    OUTPUT_PARAMS( params.outdir )
}
