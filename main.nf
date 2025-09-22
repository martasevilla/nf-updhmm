#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PREPROCESS_VCF } from './subworkflows/local/preprocess_vcf/main'

workflow {
    
    input_file = params.input ?: 'sample_sheet.csv'
    PREPROCESS_VCF(input_file)
    
    final_vcfs = PREPROCESS_VCF.out.vcfs
    final_vcfs.view { "Final processed VCF: $it" }
}