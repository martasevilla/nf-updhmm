#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(dplyr)
    library(optparse)
})

option_list <- list(
    make_option(c("-i", "--input"), 
                type = "character", 
                default = NULL,
                help = "Archivo de entrada con eventos UPD (.txt o .rds)", 
                metavar = "character"),
    make_option(c("-o", "--output"), 
                type = "character", 
                default = "filtered_upd_events.txt",
                help = "Archivo de salida [default= %default]", 
                metavar = "character"),
    make_option(c("--min_errors"), 
                type = "integer", 
                default = 2,
                help = "Número mínimo de errores mendelianos [default= %default]", 
                metavar = "integer"),
    make_option(c("--min_size"), 
                type = "integer", 
                default = 500000,
                help = "Tamaño mínimo del evento en bp [default= %default]", 
                metavar = "integer"),
    make_option(c("-v", "--verbose"), 
                action = "store_true", 
                default = TRUE,
                help = "Imprimir mensajes detallados [default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Debe especificar un archivo de entrada con --input", call. = FALSE)
}

if (!file.exists(opt$input)) {
    stop(paste("El archivo de entrada no existe:", opt$input), call. = FALSE)
}


read_upd_file <- function(file_path) {
    if (grepl("\\.rds$", file_path, ignore.case = TRUE)) {
        return(readRDS(file_path))
    } else if (grepl("\\.txt$", file_path, ignore.case = TRUE)) {
        return(read.delim(file_path, stringsAsFactors = FALSE))
    } else {
        stop("Formato de archivo no soportado. Use .rds o .txt")
    }
}


mergeOverlappingEvents <- function(df) {
    
    gr <- GRanges(
        seqnames = df$seqnames,
        ranges = IRanges(start = df$start, end = df$end),
        group = df$group,
        idx = 1:nrow(df)  
    )
    
    result_list <- list()
    
    for (g in unique(gr$group)) {
        
        gr_group <- gr[gr$group == g]
        
        ol <- reduce(gr_group, with.revmap = TRUE)
        
        for (i in seq_along(ol)) {
            revmap_idx <- ol$revmap[[i]]

            fused_range <- data.frame(
                seqnames = as.character(seqnames(ol[i])),
                start = start(ol[i]),
                end = end(ol[i]),
                group = g
            )
            
            combined_info <- df[gr_group$idx[revmap_idx], ]
        
            num_cols <- c("n_snps", "log_likelihood", "p_value", "n_mendelian_error")
            num_combined <- sapply(num_cols, function(col) {
                if (col == "n_mendelian_error") {
                    sum(as.numeric(combined_info[[col]]), na.rm = TRUE)
                } else {
                    paste(combined_info[[col]], collapse = ";")
                }
            })
            
            fused_row <- cbind(fused_range, as.data.frame(t(num_combined), stringsAsFactors = FALSE))
            
            result_list[[length(result_list) + 1]] <- fused_row
        }
    }
    

    result_df <- do.call(rbind, result_list)
    for (col in num_cols) {
        result_df[[col]] <- as.character(result_df[[col]])
    }
    
    return(result_df)
}

filterLowConfidenceEvents <- function(df, min_errors, min_size) {

    filtered_df <- df %>%
        filter(
            as.numeric(n_mendelian_error) >= min_errors,
            (end - start + 1) >= min_size
        )
    
    return(filtered_df)
}

main <- function() {
    
    df <- read_upd_file(opt$input)

    merged_df <- mergeOverlappingEvents(df)

    filtered_df <- filterLowConfidenceEvents(merged_df, opt$min_errors, opt$min_size)

    write.table(filtered_df, opt$output, sep = "\t", row.names = FALSE, quote = FALSE)
}


tryCatch({
    main()
}, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    quit(status = 1)
})