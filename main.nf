#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PREPROCESS_VCF } from './subworkflows/local/preprocess_vcf/main'
include { UPDHMM_ANALYSIS } from './modules/local/updhmm_analysis/main'

workflow {
    
    input_file = params.input ?: 'sample_sheet.csv'
    PREPROCESS_VCF(input_file)
    
    final_vcfs = PREPROCESS_VCF.out.vcfs
    final_vcfs.view { "Final processed VCF: $it" }
    
    // Apply UPD analysis to the processed VCFs
    UPDHMM_ANALYSIS(final_vcfs)
    
    // View the UPD results
    UPDHMM_ANALYSIS.out.upd_results.view { "UPD results: $it" }
    UPDHMM_ANALYSIS.out.upd_events.view { "UPD events: $it" }
}
