#!/usr/bin/env Rscript
# Jumble Frankenplot - Command-line wrapper for frankenplot()
# Usage: Rscript jumble-frankenplot.R -j <tumor.jumble.csv> -c <tumor.cns> -o <output.html>

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Jumble))

# When this wrapper is run directly from an uninstalled source checkout, the
# installed Jumble package may lag behind the source tree. Fall back to loading
# the source checkout if the installed package does not yet export frankenplot().
if (!exists("frankenplot", mode = "function")) {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg) > 0) {
        script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
        package_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)
        description_file <- file.path(package_root, "DESCRIPTION")
        if (file.exists(description_file) && requireNamespace("devtools", quietly = TRUE)) {
            suppressPackageStartupMessages(devtools::load_all(package_root, quiet = TRUE))
        }
    }
}

if (!exists("frankenplot", mode = "function")) {
    stop("frankenplot() is not available. Install the current Jumble package or run the wrapper from a source checkout with devtools available.")
}

# Define command-line options
option_list <- list(
    make_option(c("-j", "--jumble-csv"),
        type = "character", default = NULL,
        help = "Tumor .jumble.csv file (required)", metavar = "FILE"
    ),
    make_option(c("-c", "--cns"),
        type = "character", default = NULL,
        help = "Tumor .cns segment file (required)", metavar = "FILE"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = NULL,
        help = "Output HTML report file (required)", metavar = "FILE"
    ),
    make_option("--normal-jumble-csv",
        type = "character", default = NULL,
        help = "Normal .jumble.csv file (optional)", metavar = "FILE"
    ),
    make_option("--normal-cns",
        type = "character", default = NULL,
        help = "Normal .cns segment file (optional)", metavar = "FILE"
    ),
    make_option("--tumor-snp-vcf",
        type = "character", default = NULL,
        help = "Tumor SNP VCF for allele-ratio plotting (optional)", metavar = "FILE"
    ),
    make_option("--normal-snp-vcf",
        type = "character", default = NULL,
        help = "Normal SNP VCF for allele-ratio plotting (optional)", metavar = "FILE"
    ),
    make_option(c("-s", "--somatic-vcf"),
        type = "character", default = NULL,
        help = "Somatic mutation VCF for variant overlay (optional)", metavar = "FILE"
    ),
    make_option("--germline-vcf",
        type = "character", default = NULL,
        help = "Germline mutation VCF for variant overlay (optional)", metavar = "FILE"
    ),
    make_option("--tumor-sample-name",
        type = "character", default = NULL,
        help = "Tumor sample name hint for VCF sample-column selection (optional)", metavar = "NAME"
    ),
    make_option("--normal-sample-name",
        type = "character", default = NULL,
        help = "Normal sample name hint for VCF sample-column selection (optional)", metavar = "NAME"
    ),
    make_option(c("-g", "--genome"),
        type = "character", default = NULL,
        help = "Genome build: 'hg19' or 'hg38' [default: auto-detect]", metavar = "GENOME"
    ),
    make_option("--hrdtable",
        type = "character", default = NULL,
        help = "Precomputed GIS/HRD table file readable by data.table::fread() (optional)", metavar = "FILE"
    ),
    make_option("--qc-file",
        type = "character", default = NULL,
        help = "QC metrics CSV file (optional; auto-detected from tumor .jumble.csv when omitted)", metavar = "FILE"
    ),
    make_option("--dpyd-json",
        type = "character", default = NULL,
        help = "DPYD JSON result file (optional)", metavar = "FILE"
    ),
    make_option("--dpyd-csv",
        type = "character", default = NULL,
        help = "DPYD CSV evidence file (optional)", metavar = "FILE"
    )
)

opt_parser <- OptionParser(
    usage = "Usage: %prog [options]",
    option_list = option_list,
    description = "\nGenerate a Frankenplot HTML genome report from Jumble output files.\n\nExamples:\n  Rscript jumble-frankenplot.R -j sample.jumble.csv -c sample.cns -o sample.frankenplot.html\n  Rscript jumble-frankenplot.R -j tumor.jumble.csv -c tumor.cns -o report.html --tumor-snp-vcf sample.vcf.gz --somatic-vcf sample.somatic.vep.vcf.gz\n  Rscript jumble-frankenplot.R -j tumor.jumble.csv -c tumor.cns -o report.html --normal-jumble-csv normal.jumble.csv --normal-cns normal.cns --tumor-snp-vcf paired.vcf.gz --normal-snp-vcf paired.vcf.gz"
)

opt <- parse_args(opt_parser)

# Validate required inputs
if (is.null(opt$`jumble-csv`)) {
    stop("Tumor .jumble.csv file (--jumble-csv) is required")
}

if (is.null(opt$cns)) {
    stop("Tumor .cns segment file (--cns) is required")
}

if (is.null(opt$output)) {
    stop("Output HTML report file (--output) is required")
}

# Validate required file existence
if (!file.exists(opt$`jumble-csv`)) {
    stop("Tumor .jumble.csv file not found: ", opt$`jumble-csv`)
}

if (!file.exists(opt$cns)) {
    stop("Tumor .cns segment file not found: ", opt$cns)
}

# Validate optional file existence
optional_files <- list(
    `--normal-jumble-csv` = opt$`normal-jumble-csv`,
    `--normal-cns` = opt$`normal-cns`,
    `--tumor-snp-vcf` = opt$`tumor-snp-vcf`,
    `--normal-snp-vcf` = opt$`normal-snp-vcf`,
    `--somatic-vcf` = opt$`somatic-vcf`,
    `--germline-vcf` = opt$`germline-vcf`,
    `--hrdtable` = opt$hrdtable,
    `--qc-file` = opt$`qc-file`,
    `--dpyd-json` = opt$`dpyd-json`,
    `--dpyd-csv` = opt$`dpyd-csv`
)

for (nm in names(optional_files)) {
    path <- optional_files[[nm]]
    if (!is.null(path) && !file.exists(path)) {
        stop(nm, " file not found: ", path)
    }
}

# Require paired normal copy-number inputs when normal data are supplied.
has_normal_csv <- !is.null(opt$`normal-jumble-csv`)
has_normal_cns <- !is.null(opt$`normal-cns`)
if (xor(has_normal_csv, has_normal_cns)) {
    stop("Normal copy-number input requires both --normal-jumble-csv and --normal-cns")
}

# Validate genome if supplied.
if (!is.null(opt$genome) && !(opt$genome %in% c("hg19", "hg38"))) {
    stop("Genome must be 'hg19' or 'hg38' when supplied")
}

# Create output directory when needed.
output_dir <- dirname(opt$output)
if (!identical(output_dir, ".") && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Read externally supplied HRD/GIS table if supplied.
hrdtable <- NULL
if (!is.null(opt$hrdtable)) {
    hrdtable <- data.table::fread(opt$hrdtable)
}

cat("Jumble Frankenplot Report\n")
cat("=========================\n\n")
cat("Tumor CSV:", opt$`jumble-csv`, "\n")
cat("Tumor CNS:", opt$cns, "\n")
if (has_normal_csv) {
    cat("Normal CSV:", opt$`normal-jumble-csv`, "\n")
    cat("Normal CNS:", opt$`normal-cns`, "\n")
}
if (!is.null(opt$`tumor-snp-vcf`)) {
    cat("Tumor SNP VCF:", opt$`tumor-snp-vcf`, "\n")
}
if (!is.null(opt$`normal-snp-vcf`)) {
    cat("Normal SNP VCF:", opt$`normal-snp-vcf`, "\n")
}
if (!is.null(opt$`somatic-vcf`)) {
    cat("Somatic VCF:", opt$`somatic-vcf`, "\n")
}
if (!is.null(opt$`germline-vcf`)) {
    cat("Germline VCF:", opt$`germline-vcf`, "\n")
}
if (!is.null(opt$hrdtable)) {
    cat("HRD/GIS table:", opt$hrdtable, "\n")
}
if (!is.null(opt$`qc-file`)) {
    cat("QC file:", opt$`qc-file`, "\n")
}
if (!is.null(opt$`dpyd-json`)) {
    cat("DPYD JSON:", opt$`dpyd-json`, "\n")
}
if (!is.null(opt$`dpyd-csv`)) {
    cat("DPYD CSV:", opt$`dpyd-csv`, "\n")
}
cat("Genome:", if (is.null(opt$genome)) "auto" else opt$genome, "\n")
cat("Output:", opt$output, "\n\n")

# Generate report.
frankenplot(
    tumor_jumble_csv = opt$`jumble-csv`,
    tumor_cns = opt$cns,
    output_file = opt$output,
    normal_jumble_csv = opt$`normal-jumble-csv`,
    normal_cns = opt$`normal-cns`,
    tumor_snp_vcf = opt$`tumor-snp-vcf`,
    normal_snp_vcf = opt$`normal-snp-vcf`,
    somatic_vcf = opt$`somatic-vcf`,
    germline_vcf = opt$`germline-vcf`,
    tumor_sample_name = opt$`tumor-sample-name`,
    normal_sample_name = opt$`normal-sample-name`,
    genome = opt$genome,
    hrdtable = hrdtable,
    qc_file = opt$`qc-file`,
    dpyd_json = opt$`dpyd-json`,
    dpyd_csv = opt$`dpyd-csv`
)

cat("\n✓ Frankenplot report generated\n")
cat("  - File:", opt$output, "\n")
cat("\nDone.\n")
