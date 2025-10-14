#!/usr/bin/env nextflow

process CALCULATE_EVENTS {
    tag "$meta.id"
    label 'process_high'
    
    container "/home/u0030001/nf-updhmm_zenodo/updhmm-new_1.3.2.sif"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.upd_results.txt"),   emit: upd_results
    tuple val(meta), path("*.upd_collapsed.txt"), emit: upd_collapsed
    tuple val(meta), path("*.upd_events.rds"),    emit: upd_events 
    tuple val(meta), path("*.upd_collapsed.rds"), emit: upd_collapsed_rds

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genome_build = meta.genome_build ?: params.genome_build ?: "hg38"
    def verbose = task.ext.verbose ? "--verbose" : ""
    
    """
    updhmm_analysis.r \\
        --input ${vcf} \\
        --output_prefix ${prefix} \\
        --genome_build ${genome_build} \\
        ${verbose} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.upd_results.txt
    touch ${prefix}.upd_collapsed.txt
    touch ${prefix}.upd_events.rds  
    touch ${prefix}.upd_collapsed.rds
    
    """
}
