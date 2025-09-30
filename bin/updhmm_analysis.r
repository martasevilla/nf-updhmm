#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(UPDhmm)
  library(VariantAnnotation)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"), 
              type = "character", 
              default = NULL,
              help = "Archivo VCF de entrada", 
              metavar = "character"),
  make_option(c("-o", "--output_prefix"), 
              type = "character", 
              default = "sample",
              help = "Prefijo para archivos de salida [default= %default]", 
              metavar = "character"),
  make_option(c("-g", "--genome_build"), 
              type = "character", 
              default = "hg38",
              help = "Versión del genoma (hg19, hg38) [default= %default]", 
              metavar = "character"),
  make_option(c("-v", "--verbose"), 
              action = "store_true", 
              default = TRUE,
              help = "Imprimir mensajes detallados [default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Debe especificar un archivo VCF de entrada con --input", call. = FALSE)
}

if (!file.exists(opt$input)) {
  stop(paste("El archivo VCF no existe:", opt$input), call. = FALSE)
}


tryCatch({

  vcf <- readVcf(opt$input, opt$genome_build)
  sample_names <- colnames(geno(vcf)$GT)
  
  processedVcf <- vcfCheck(vcf, 
                           proband = sample_names[1], 
                           mother = sample_names[2], 
                           father = sample_names[3])
  
  updEvents <- calculateEvents(processedVcf)
  
  results_file <- paste0(opt$output_prefix, ".upd_results.txt")
  rds_file <- paste0(opt$output_prefix, ".upd_events.rds")
  
  write.table(updEvents, file = results_file, sep = "\t", row.names = FALSE, quote = FALSE)
  saveRDS(updEvents, file = rds_file)
  
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  quit(status = 1)
})