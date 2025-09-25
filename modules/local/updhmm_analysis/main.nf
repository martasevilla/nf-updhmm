process UPDHMM_ANALYSIS {
    tag "$meta.id"
    label 'process_medium'
    
    conda "bioconda::bioconductor-updhmm=1.4.0"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.upd_results.txt"), emit: upd_results
    tuple val(meta), path("*.upd_events.rds"),  emit: upd_events
    path "*.log",                               optional: true, emit: log
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genome_build = meta.genome_build ?: params.genome_build ?: "hg38"
    
    """
    #!/usr/bin/env Rscript

    # Load required libraries
    
    library(UPDhmm)
    library(VariantAnnotation)
    

    # Log file setup
    log_file <- "${prefix}.log"
    cat("Starting UPD analysis for sample: ${meta.id}\\n", file = log_file)
    
    tryCatch({
        # Read VCF file
        cat("Reading VCF file: ${vcf}\\n", file = log_file, append = TRUE)
        vcf <- readVcf("${vcf}", "${genome_build}")
        
        # Get sample names
        sample_names <- colnames(geno(vcf)\$GT)
        cat("Found samples:", paste(sample_names, collapse = ", "), "\\n", 
            file = log_file, append = TRUE)
        
        # Validate trio structure
        if(length(sample_names) != 3) {
            stop("VCF must contain exactly 3 samples for trio analysis. Found: ", 
                 length(sample_names), " samples")
        }
        
        # Define trio roles (assuming order: proband, mother, father)
        # This could be made configurable via meta or params
        proband_id <- sample_names[1]
        mother_id <- sample_names[2] 
        father_id <- sample_names[3]
        
        cat("Trio configuration:\\n", file = log_file, append = TRUE)
        cat("  Proband:", proband_id, "\\n", file = log_file, append = TRUE)
        cat("  Mother:", mother_id, "\\n", file = log_file, append = TRUE)
        cat("  Father:", father_id, "\\n", file = log_file, append = TRUE)
        
        # Process VCF for UPD analysis
        cat("Processing VCF for UPD analysis...\\n", file = log_file, append = TRUE)
        processedVcf <- vcfCheck(vcf, 
                                proband = proband_id, 
                                mother = mother_id, 
                                father = father_id)
        
        # Calculate UPD events
        cat("Calculating UPD events...\\n", file = log_file, append = TRUE)
        updEvents <- calculateEvents(processedVcf)
        
        # Write results
        cat("Writing results...\\n", file = log_file, append = TRUE)
        write.table(updEvents, 
                    file = "${prefix}.upd_results.txt", 
                    sep = "\\t", 
                    row.names = FALSE, 
                    quote = FALSE)
        
        # Save R object
        saveRDS(updEvents, file = "${prefix}.upd_events.rds")
        
        cat("UPD analysis completed successfully!\\n", file = log_file, append = TRUE)
        cat("Results written to: ${prefix}.upd_results.txt\\n", file = log_file, append = TRUE)
        cat("R object saved to: ${prefix}.upd_events.rds\\n", file = log_file, append = TRUE)
        
    }, error = function(e) {
        cat("ERROR:", e\$message, "\\n", file = log_file, append = TRUE)
        stop("UPD analysis failed: ", e\$message)
    })
    
    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-updhmm: \$(Rscript -e "library(UPDhmm); cat(as.character(packageVersion('UPDhmm')))" 2>/dev/null)
        bioconductor-variantannotation: \$(Rscript -e "library(VariantAnnotation); cat(as.character(packageVersion('VariantAnnotation')))" 2>/dev/null)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create mock output files for testing
    echo -e "chromosome\\tstart\\tend\\tgroup\\tn_snps\\tlog_likelihood\\tp_value\\tn_mendelian_error" > ${prefix}.upd_results.txt
    echo -e "chr1\\t1000000\\t2000000\\tmatUPD\\t150\\t-45.2\\t0.001\\t5" >> ${prefix}.upd_results.txt
    
    # Create empty RDS file
    touch ${prefix}.upd_events.rds
    
    # Create log file
    echo "Stub run completed for ${meta.id}" > ${prefix}.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.0
        bioconductor-updhmm: 1.4.0
        bioconductor-variantannotation: 1.52.0
    END_VERSIONS
    """
}
