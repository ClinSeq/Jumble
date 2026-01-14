#!/usr/bin/env Rscript
# Jumble Reference - Command-line wrapper for build_reference()
# Usage: Rscript jumble-reference.R -i <counts_dir> -a <annotation_source> -o <output_file> -c <cores>

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Jumble))

# Define command-line options
option_list <- list(
    make_option(c("-i", "--input"),
        type = "character", default = NULL,
        help = "Directory containing .counts.RDS files (required)", metavar = "DIR"
    ),
    make_option(c("-a", "--annotation"),
        type = "character", default = "biomart",
        help = "Annotation source: 'biomart' or path to annotation folder [default: %default]", metavar = "SOURCE"
    ),
    make_option(c("-g", "--genome"),
        type = "character", default = NULL,
        help = "Genome version ('hg19' or 'hg38'). If omitted, detected from counts [default: auto]", metavar = "GENOME"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = NULL,
        help = "Output file path [default: auto-generated based on BED file]", metavar = "FILE"
    ),
    make_option(c("-c", "--cores"),
        type = "integer", default = 1,
        help = "Number of cores for parallel processing [default: %default]", metavar = "N"
    )
)

opt_parser <- OptionParser(
    usage = "Usage: %prog [options]",
    option_list = option_list,
    description = "\nBuild a reference panel from count files.\n\nExample:\n  Rscript jumble-reference.R -i reference_counts/ -a biomart -o reference.RDS -c 4"
)

opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$input)) {
    stop("Input directory (--input) is required")
}

if (!dir.exists(opt$input)) {
    stop("Input directory not found: ", opt$input)
}

# Find count files
count_files <- list.files(opt$input, pattern = "\\.counts\\.RDS$", full.names = TRUE)

if (length(count_files) == 0) {
    stop("No .counts.RDS files found in: ", opt$input)
}

cat("Found", length(count_files), "count files\n")

if (length(count_files) < 4) {
    warning("Recommended minimum is 4 reference samples. Found: ", length(count_files))
}

# Determine output file
if (is.null(opt$output)) {
    # Auto-generate based on first count file's BED file
    first_count <- readRDS(count_files[1])
    bed_name <- if (!is.null(first_count$target_bed_file)) {
        basename(first_count$target_bed_file)
    } else {
        "wgs"
    }
    opt$output <- paste0(tools::file_path_sans_ext(bed_name), ".reference.RDS")
}

cat("Output file:", opt$output, "\n")
cat("Annotation source:", opt$annotation, "\n")
cat("Genome:", opt$genome, "\n")
cat("Cores:", opt$cores, "\n\n")

# Build reference
cat("Building reference...\n")

reference <- build_reference(
    count_files = count_files,
    annotation_source = opt$annotation,
    genome = opt$genome,
    output_file = opt$output,
    cores = opt$cores
)

cat("\n✓ Reference built successfully\n")
cat("  - Bins:", nrow(reference$target_template), "\n")
cat("  - Samples:", length(reference$samples), "\n")
cat("  - File:", opt$output, "\n")

cat("\nDone.\n")
