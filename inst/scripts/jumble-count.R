#!/usr/bin/env Rscript
# Jumble Count - Command-line wrapper for generate_counts()
# Usage: Rscript jumble-count.R -t <target_BED_file> -b <input_bam_file>
#    or: Rscript jumble-count.R -t <target_BED_file> -c <cores>

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Jumble))

# Define command-line options
option_list <- list(
    make_option(c("-t", "--target"),
        type = "character", default = NULL,
        help = "Target BED file (required for targeted sequencing, omit for WGS)", metavar = "FILE"
    ),
    make_option(c("-b", "--bam"),
        type = "character", default = NULL,
        help = "Input BAM file (sorted and indexed)", metavar = "FILE"
    ),
    make_option(c("-c", "--cores"),
        type = "integer", default = 1,
        help = "Number of cores for parallel processing [default: %default]", metavar = "N"
    ),
    make_option(c("-w", "--wgs-bin-size"),
        type = "integer", default = 10000,
        help = "Bin size for WGS mode [default: %default]", metavar = "N"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = NULL,
        help = "Output directory [default: current directory]", metavar = "DIR"
    )
)

opt_parser <- OptionParser(
    usage = "Usage: %prog [options]",
    option_list = option_list,
    description = "\nGenerate read count files for Jumble analysis.\n\nExamples:\n  Rscript jumble-count.R -t targets.bed -b sample.bam\n  Rscript jumble-count.R -t targets.bed -c 4  # Process all BAMs in current dir\n  Rscript jumble-count.R -b sample.bam -w 10000  # WGS mode"
)

opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$target) && is.null(opt$`wgs-bin-size`)) {
    stop("Either --target (for targeted sequencing) or --wgs-bin-size (for WGS) must be specified")
}

# Set output directory
output_dir <- if (!is.null(opt$output)) opt$output else getwd()
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Process single BAM or all BAMs in directory
if (!is.null(opt$bam)) {
    # Single BAM mode
    if (!file.exists(opt$bam)) {
        stop("BAM file not found: ", opt$bam)
    }

    cat("Processing BAM file:", opt$bam, "\n")

    counts <- generate_counts(
        bam_file = opt$bam,
        target_bed = opt$target,
        wgs_bin_size = if (is.null(opt$target)) opt$`wgs-bin-size` else NULL
    )

    # Save counts
    output_file <- file.path(output_dir, paste0(basename(opt$bam), ".counts.RDS"))
    saveRDS(counts, output_file)
    cat("Saved counts to:", output_file, "\n")
} else {
    # Multi-BAM mode (process all BAMs in current directory)
    bam_files <- list.files(pattern = "\\.bam$", full.names = TRUE)

    if (length(bam_files) == 0) {
        stop("No BAM files found in current directory")
    }

    cat("Found", length(bam_files), "BAM files\n")
    cat("Processing with", opt$cores, "cores\n\n")

    if (opt$cores > 1) {
        suppressPackageStartupMessages(library(parallel))

        mclapply(bam_files, function(bam_file) {
            cat("Processing:", bam_file, "\n")

            counts <- generate_counts(
                bam_file = bam_file,
                target_bed = opt$target,
                wgs_bin_size = if (is.null(opt$target)) opt$`wgs-bin-size` else NULL
            )

            output_file <- file.path(output_dir, paste0(basename(bam_file), ".counts.RDS"))
            saveRDS(counts, output_file)
            cat("Saved:", output_file, "\n")

            invisible(NULL)
        }, mc.cores = opt$cores)
    } else {
        # Sequential processing
        for (bam_file in bam_files) {
            cat("Processing:", bam_file, "\n")

            counts <- generate_counts(
                bam_file = bam_file,
                target_bed = opt$target,
                wgs_bin_size = if (is.null(opt$target)) opt$`wgs-bin-size` else NULL
            )

            output_file <- file.path(output_dir, paste0(basename(bam_file), ".counts.RDS"))
            saveRDS(counts, output_file)
            cat("Saved:", output_file, "\n")
        }
    }

    cat("\nProcessed", length(bam_files), "BAM files\n")
}

cat("\nDone.\n")
