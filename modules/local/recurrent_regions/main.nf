#!/usr/bin/env nextflow

process RECURRENT_REGIONS {
    tag "recurrent_analysis"
    label 'process_high'
    
    container "/home/u0030001/nf-updhmm_zenodo/updhmm-new_1.3.2.sif"

    input:
    path(collapsed_rds_files)

    output:
    path("*.recurrent_regions.txt"), emit: recurrent_regions
    path("*.recurrent_regions.rds"),  emit: recurrent_rds

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "cohort"
    def min_support = params.min_support ?: 2
    def error_threshold = params.error_threshold ? params.error_threshold : "NULL"
    def verbose = task.ext.verbose ? "--verbose" : ""
    
    def file_list = collapsed_rds_files instanceof List ? collapsed_rds_files.join(',') : collapsed_rds_files
    
    """
    recurrent_regions.r \\
        --input_files "${file_list}" \\
        --output_prefix ${prefix} \\
        --min_support ${min_support} \\
        ${error_threshold != "NULL" ? "--error_threshold ${error_threshold}" : ""} \\
        ${verbose} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "cohort"
    """
    touch ${prefix}.recurrent_regions.txt
    touch ${prefix}.recurrent_regions.rds
    
    """
}
