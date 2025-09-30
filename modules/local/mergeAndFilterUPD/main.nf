#!/usr/bin/env nextflow

process MERGEANDFILTERUPD {
    tag "$meta.id"
    label 'process_medium'
    
    conda "bioconda::bioconductor-genomicranges=1.54.1 conda-forge::r-dplyr=1.1.4 conda-forge::r-optparse=1.7.3"

    input:
    tuple val(meta), path(upd_events_file)

    output:
    tuple val(meta), path("*_filtered_upd_events.txt"), emit: filtered_events
    path "versions.yml", emit: versions
   
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: (meta.id ?: "sample")
    def min_errors = task.ext.min_errors ?: 2
    def min_size = task.ext.min_size ?: 500000
    
    """
    merge_and_filter_upd.r \\
        --input ${upd_events_file} \\
        --output ${prefix}_filtered_upd_events.txt \\
        --min_errors ${min_errors} \\
        --min_size ${min_size} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/R version \\([0-9.]*\\).*/\\1/')
        bioconductor-genomicranges: \$(Rscript -e "cat(as.character(packageVersion('GenomicRanges')))")
        r-dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: (meta.id ?: "sample")
    """
    touch ${prefix}_filtered_upd_events.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.0
        bioconductor-genomicranges: 1.54.1
        r-dplyr: 1.1.4
        r-optparse: 1.7.3
    END_VERSIONS
    """
}