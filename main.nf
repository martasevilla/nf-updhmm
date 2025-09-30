#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PREPROCESS_VCF } from './subworkflows/local/preprocess_vcf/main'
include { CALCULATE_EVENTS } from './modules/local/updhmm_analysis/main'
include { MERGEANDFILTERUPD } from './modules/local/mergeAndFilterUPD'

params.upd_events = "NDT1544_updEvents.txt"

workflow {
    
    input_file = params.input ?: 'sample_sheet.csv'
    PREPROCESS_VCF(input_file)
    
    final_vcfs = PREPROCESS_VCF.out.vcfs
    final_vcfs.view { "Final processed VCF: $it" }
    
    //def meta = [id: 'NDT1544']

    //upd_events_ch = Channel.of([meta, file(params.upd_events)])
    //MERGEANDFILTERUPD(upd_events_ch)
    //MERGEANDFILTERUPD.out.filtered_events.view { "Filtered UPD events: $it" }
    
    // Apply UPD analysis to the processed VCFs
    CALCULATE_EVENTS(final_vcfs)
    
    // View the UPD results
    CALCULATE_EVENTS.out.upd_results.view { "UPD results: $it" }
    CALCULATE_EVENTS.out.upd_events.view { "UPD events: $it" }
}