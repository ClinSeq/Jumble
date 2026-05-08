#!/usr/bin/env Rscript
# Jumble Run - Command-line wrapper for run_jumble()
# Usage: Rscript jumble-run.R -r <reference_file> -b <input_bam_file> -o <output_dir>

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Jumble))

# Define command-line options
option_list <- list(
    make_option(c("-r", "--reference"),
        type = "character", default = NULL,
        help = "Reference file (required)", metavar = "FILE"
    ),
    make_option(c("-b", "--bam"),
        type = "character", default = NULL,
        help = "Input BAM file or .counts.RDS file (required)", metavar = "FILE"
    ),
    make_option(c("-v", "--vcf"),
        type = "character", default = NULL,
        help = "VCF file for SNP analysis (optional)", metavar = "FILE"
    ),
    make_option(c("-s", "--somatic"),
        type = "character", default = NULL,
        help = "Somatic VCF file for MSI analysis (optional)", metavar = "FILE"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = NULL,
        help = "Output directory [default: current directory]", metavar = "DIR"
    ),
    make_option(c("-m", "--model"),
        type = "character", default = NULL,
        help = "HRD model file for GIS Custom HRD scoring (optional)", metavar = "FILE"
    ),
    make_option(c("-a", "--alpha"),
        type = "double", default = 0.001,
        help = "Segmentation alpha parameter [default: %default]", metavar = "N"
    ),
    make_option(c("-c", "--correction"),
        type = "character", default = "optim",
        help = "Correction method to use: 'optim' (L1+TV, default) or 'rlm' (Robust LM)", metavar = "METHOD"
    ),
    make_option("--exclude-long-fragments",
        action = "store_true", default = FALSE,
        help = "Exclude long fragments (>300bp) from depth signal, using only fragments ≤300bp. Recommended for clipoverlap BAMs."
    )
)

opt_parser <- OptionParser(
    usage = "Usage: %prog [options]",
    option_list = option_list,
    description = "\nRun Jumble copy number analysis.\n\nExamples:\n  Rscript jumble-run.R -r reference.RDS -b sample.bam -c optim -o output/\n  Rscript jumble-run.R -r reference.RDS -b sample.counts.RDS -v sample.vcf.gz"
)

opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$reference)) {
    stop("Reference file (--reference) is required")
}

if (is.null(opt$bam)) {
    stop("Input BAM or counts file (--bam) is required")
}

if (!file.exists(opt$reference)) {
    stop("Reference file not found: ", opt$reference)
}

if (!file.exists(opt$bam)) {
    stop("Input file not found: ", opt$bam)
}

if (!is.null(opt$vcf) && !file.exists(opt$vcf)) {
    stop("VCF file not found: ", opt$vcf)
}

if (!is.null(opt$somatic) && !file.exists(opt$somatic)) {
    stop("Somatic VCF file not found: ", opt$somatic)
}

if (!is.null(opt$model) && !file.exists(opt$model)) {
    warning("HRD Model file not found: ", opt$model)
}

# Set output directory
output_dir <- if (!is.null(opt$output)) opt$output else getwd()
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

cat("Jumble Copy Number Analysis\n")
cat("===========================\n\n")
cat("Reference:", opt$reference, "\n")
cat("Input:", opt$bam, "\n")
if (!is.null(opt$vcf)) {
    cat("VCF:", opt$vcf, "\n")
}
if (!is.null(opt$somatic)) {
    cat("Somatic VCF:", opt$somatic, "\n")
}
if (!is.null(opt$model)) {
    cat("Model:", opt$model, "\n")
}
cat("Output:", output_dir, "\n")
cat("Alpha:", opt$alpha, "\n")
cat("Correction:", opt$correction, "\n")
if (opt$exclude_long_fragments) {
    cat("Fragment length filter: excluding >300bp\n")
}
cat("\n")

# Run Jumble
result <- run_jumble(
    bam_file = opt$bam,
    reference_file = opt$reference,
    output_dir = output_dir,
    snp_vcf = opt$vcf,
    somatic_vcf = opt$somatic,
    alpha = opt$alpha,
    correction = opt$correction,
    hrd_model_file = opt$model,
    exclude_long_fragments = opt$exclude_long_fragments
)

cat("\n✓ Analysis complete\n")
cat("  - Segments:", nrow(result$segments), "\n")
cat("  - X chromosome segments:", sum(result$segments$chromosome == "X", na.rm = TRUE), "\n")

# List output files
output_files <- list.files(output_dir, pattern = basename(tools::file_path_sans_ext(opt$bam)), full.names = FALSE)
cat("\nOutput files:\n")
for (f in output_files) {
    cat("  -", f, "\n")
}

cat("\nDone.\n")
