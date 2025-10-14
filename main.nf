#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PREPROCESS_VCF } from './subworkflows/local/preprocess_vcf/main'
include { CALCULATE_EVENTS } from './modules/local/updhmm_analysis/main'
//include { RECURRENT_REGIONS  } from './modules/local/recurrent_regions/main'

workflow {
    
    // Step 0: Define input samplesheet (default to 'sample_sheet.csv' if not provided)
    input_file = params.input ?: 'sample_sheet.csv'
    
    // Step 1: Preprocess VCFs (validation, annotation removal, merging, filtering)
    PREPROCESS_VCF(input_file)
    
    // Step 2: Apply UPD analysis to the processed VCFs
    CALCULATE_EVENTS(PREPROCESS_VCF.out.vcfs)
    
    // Step 3: POST-PROCESSING - Mark recurrent regions using RDS files
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

