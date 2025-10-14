#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(UPDhmm)
  #library(data.table)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input_files"), 
              type = "character", 
              default = NULL,
              help = "Comma-separated list of .upd_collapsed.rds files", 
              metavar = "character"),
  make_option(c("-o", "--output_prefix"), 
              type = "character", 
              default = "cohort",
              help = "Prefix for output files [default= %default]", 
              metavar = "character"),
  make_option(c("-s", "--min_support"), 
              type = "integer", 
              default = 2,
              help = "Minimum number of samples to consider a region recurrent [default= %default]", 
              metavar = "integer"),
  make_option(c("-e", "--error_threshold"), 
              type = "numeric", 
              default = NULL,
              help = "Mendelian error threshold to exclude regions [default= NULL]", 
              metavar = "numeric"),
  make_option(c("-v", "--verbose"), 
              action = "store_true", 
              default = TRUE,
              help = "Print detailed messages [default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_files)) {
  print_help(opt_parser)
  stop("You must specify input RDS files with --input_files", call. = FALSE)
}

tryCatch({

  input_files <- trimws(unlist(strsplit(opt$input_files, ",")))
  missing_files <- input_files[!file.exists(input_files)]
  if (length(missing_files) > 0) {
    stop(paste("The following files do not exist:", paste(missing_files, collapse = ", ")), call. = FALSE)
  }

  all_collapsed_data <- lapply(input_files, readRDS)
  #combined_collapsed <- rbindlist(all_collapsed_data, use.names = TRUE, fill = TRUE)
  combined_collapsed <- do.call(rbind, all_collapsed_data)

  recurrent <- markRecurrentRegions(
    subset_df = combined_collapsed,
    error_threshold = opt$error_threshold,
    min_support = opt$min_support
  )

  results_file <- paste0(opt$output_prefix, ".recurrent_regions.txt")
  rds_file <- paste0(opt$output_prefix, ".recurrent_regions.rds")

  #fwrite(recurrent, file = results_file, sep = "\t", quote = FALSE)
  write.table(recurrent, file = results_file, sep = "\t", row.names = FALSE, quote = FALSE)
  saveRDS(recurrent, file = rds_file)

}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  quit(status = 1)
})