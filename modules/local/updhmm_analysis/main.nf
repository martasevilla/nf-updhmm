#!/usr/bin/env nextflow

process CALCULATE_EVENTS {
    tag "$meta.id"
    label 'process_medium'
    
    conda "bioconda::bioconductor-updhmm=1.4.0 conda-forge::r-optparse=1.7.3"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.upd_results.txt"), emit: upd_results
    tuple val(meta), path("*.upd_events.rds"),  emit: upd_events
    path "versions.yml", emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/R version \\([0-9.]*\\).*/\\1/')
        bioconductor-updhmm: \$(Rscript -e "cat(as.character(packageVersion('UPDhmm')))")
        bioconductor-variantannotation: \$(Rscript -e "cat(as.character(packageVersion('VariantAnnotation')))")
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.upd_results.txt
    touch ${prefix}.upd_events.rds
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.0
        bioconductor-updhmm: 1.4.0
        bioconductor-variantannotation: 1.48.0
        r-optparse: 1.7.3
    END_VERSIONS
    """
}