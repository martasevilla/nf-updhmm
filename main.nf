#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PREPROCESS_VCF } from './subworkflows/local/preprocess_vcf/main'
include { VCF_CHECK } from './modules/local/vcf_check/main'
include { CALCULATE_EVENTS } from './modules/local/calculate_events/main' 
include { COLLAPSE_EVENTS } from './modules/local/collapse_events/main'
//include { CALCULATE_EVENTS } from './modules/local/updhmm_analysis/main'
//include { RECURRENT_REGIONS  } from './modules/local/recurrent_regions/main'

workflow {
    
    // Step 0: Define input samplesheet (default to 'sample_sheet.csv' if not provided)
    input_file = params.input ?: 'sample_sheet.csv'
    
    // Step 1: Preprocess VCFs (validation, annotation removal, merging, filtering)
    PREPROCESS_VCF(input_file)

    // Step 2: VCF Check (validate and prepare for UPD analysis)
    VCF_CHECK(final_vcfs)
    
    // Step 3: Calculate Events (compute UPD events)
    CALCULATE_EVENTS(VCF_CHECK.out.processed_vcf)
    CALCULATE_EVENTS.out.upd_events_txt.view { "UPD events (TXT): $it" }
    
    // Step 4: Collapse Events (merge adjacent/overlapping events)
    COLLAPSE_EVENTS(CALCULATE_EVENTS.out.upd_events_rds)
    COLLAPSE_EVENTS.out.upd_collapsed_txt.view { "Collapsed events (TXT): $it" }

    
    // Step 5: POST-PROCESSING - Recurrent regions analysis
    /**
    all_collapsed_rds = CALCULATE_EVENTS.out.upd_collapsed_rds
        .map { meta, file -> file } 
        .collect()  // Wait for ALL samples to complete
    
    all_collapsed_rds
        .filter { it.size() > 1 }  
        .set { filtered_rds_files }
    
    RECURRENT_REGIONS(filtered_rds_files)
    **/
}

