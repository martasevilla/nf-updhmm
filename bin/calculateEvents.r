#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(UPDhmm)
  library(BiocParallel)
  library(parallel)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"), 
              type = "character", 
              default = NULL,
              help = "Input processed VCF RDS file from vcfCheck", 
              metavar = "character"),
  make_option(c("-o", "--output_prefix"), 
              type = "character", 
              default = "sample",
              help = "Prefix for output files [default= %default]", 
              metavar = "character"),
  make_option(c("-c", "--cpus"), 
              type = "integer", 
              default = 1,
              help = "Number of CPUs for parallel processing [default= %default]", 
              metavar = "integer"),
  make_option(c("-v", "--verbose"), 
              action = "store_true", 
              default = TRUE,
              help = "Print detailed messages [default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("You must specify an input processed VCF RDS file with --input", call. = FALSE)
}

if (!file.exists(opt$input)) {
  stop(paste("The input RDS file does not exist:", opt$input), call. = FALSE)
}

tryCatch({

  processedVcf <- readRDS(opt$input)

  # Setup parallel processing based on available CPUs
  if (opt$cpus > 1) {
    # Detect available cores and use the minimum of requested vs available
    available_cores <- detectCores()
    workers <- min(opt$cpus, available_cores, na.rm = TRUE)
    bp_param <- MulticoreParam(workers = workers)
  } else {
    bp_param <- SerialParam()
  }

  start_time <- Sys.time()
  updEvents <- calculateEvents(processedVcf, BPPARAM = bp_param)
  end_time <- Sys.time()

  if (opt$verbose) {
    cat("=== UPD ANALYSIS COMPLETED ===\n")
    cat("Processing time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")}
  
  results_file <- paste0(opt$output_prefix, ".upd_events.txt")
  rds_file <- paste0(opt$output_prefix, ".upd_events.rds")

  write.table(updEvents, file = results_file, sep = "\t", row.names = FALSE, quote = FALSE)
  saveRDS(updEvents, file = rds_file)
  
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  quit(status = 1)
})