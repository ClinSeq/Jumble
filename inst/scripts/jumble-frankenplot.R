#!/usr/bin/env Rscript
# Jumble Frankenplot - Command-line wrapper for frankenplot()
# Usage: Rscript jumble-frankenplot.R -j <tumor.jumble.csv> -c <tumor.cns> -o <output.html>

suppressPackageStartupMessages(library(optparse))

# Prefer the installed package, but allow source-checkout validation before
# installation. This keeps wrapper tests active during local development.
if (requireNamespace("Jumble", quietly = TRUE)) {
    suppressPackageStartupMessages(library(Jumble))
}

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
    make_option("--output-png",
        type = "character", default = NULL,
        help = "Static tumor overview PNG file (optional; defaults to output HTML path with .png extension)", metavar = "FILE"
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
    make_option("--compute-gis",
        action = "store_true", default = FALSE,
        help = "Compute the GIS table from the tumor .jumble.csv and --tumor-snp-vcf (requires --tumor-snp-vcf and --genome; multi-sample VCFs also require --tumor-sample-name). Cannot be combined with --hrdtable"
    ),
    make_option("--gis-output",
        type = "character", default = NULL,
        help = "Output CSV for the computed GIS table (only used with --compute-gis; defaults to output HTML path with .gis.csv extension)", metavar = "FILE"
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
    description = "\nGenerate a Frankenplot HTML genome report from Jumble output files.\n\nExamples:\n  Rscript jumble-frankenplot.R -j sample.jumble.csv -c sample.cns -o sample.frankenplot.html\n  Rscript jumble-frankenplot.R -j tumor.jumble.csv -c tumor.cns -o report.html --tumor-snp-vcf sample.vcf.gz --somatic-vcf sample.somatic.vep.vcf.gz\n  Rscript jumble-frankenplot.R -j tumor.jumble.csv -c tumor.cns -o report.html --normal-jumble-csv normal.jumble.csv --normal-cns normal.cns --tumor-snp-vcf paired.vcf.gz --normal-snp-vcf paired.vcf.gz\n  Rscript jumble-frankenplot.R -j tumor.jumble.csv -c tumor.cns -o report.html --tumor-snp-vcf sample.vcf.gz --genome hg38 --compute-gis --gis-output report.gis.csv"
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

# GIS is only computed when explicitly requested with --compute-gis; without
# the flag the wrapper behaves exactly as before (no extra validation or output).
compute_gis <- isTRUE(opt$`compute-gis`)
if (!compute_gis && !is.null(opt$`gis-output`)) {
    stop("--gis-output requires --compute-gis")
}
gis_output <- NULL
if (compute_gis) {
    if (!is.null(opt$hrdtable)) {
        stop("--compute-gis cannot be combined with --hrdtable")
    }
    if (is.null(opt$`tumor-snp-vcf`)) {
        stop("--tumor-snp-vcf is required with --compute-gis")
    }
    if (is.null(opt$genome)) {
        stop("--genome is required with --compute-gis")
    }

    # The tumor sample must be identified unambiguously: a supplied name must
    # exactly match one VCF column, and multi-sample VCFs require a name.
    vcf_samples <- VariantAnnotation::samples(
        VariantAnnotation::scanVcfHeader(opt$`tumor-snp-vcf`)
    )
    tumor_name <- opt$`tumor-sample-name`
    if (length(vcf_samples) > 1 && is.null(tumor_name)) {
        stop("--tumor-sample-name is required with --compute-gis for a multi-sample VCF. ",
             "Available samples: ", paste(vcf_samples, collapse = ", "))
    }
    if (!is.null(tumor_name) && sum(vcf_samples == tumor_name) != 1) {
        stop("--tumor-sample-name '", tumor_name, "' must exactly match one sample in ",
             opt$`tumor-snp-vcf`, ". Available samples: ", paste(vcf_samples, collapse = ", "))
    }

    gis_output <- opt$`gis-output`
    if (is.null(gis_output)) {
        gis_output <- sub("\\.html?$", ".gis.csv", opt$output, ignore.case = TRUE)
        if (identical(gis_output, opt$output)) {
            gis_output <- paste0(opt$output, ".gis.csv")
        }
    }
}

# Create output directories when needed.
output_dir <- dirname(opt$output)
if (!identical(output_dir, ".") && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}
if (!is.null(opt$`output-png`)) {
    output_png_dir <- dirname(opt$`output-png`)
    if (!identical(output_png_dir, ".") && !dir.exists(output_png_dir)) {
        dir.create(output_png_dir, recursive = TRUE)
    }
}
if (!is.null(gis_output)) {
    gis_output_dir <- dirname(gis_output)
    if (!identical(gis_output_dir, ".") && !dir.exists(gis_output_dir)) {
        dir.create(gis_output_dir, recursive = TRUE)
    }
}

# Read externally supplied HRD/GIS table, or compute it when --compute-gis is set.
hrdtable <- NULL
gis_written <- FALSE
if (!is.null(opt$hrdtable)) {
    hrdtable <- data.table::fread(opt$hrdtable)
} else if (compute_gis) {
    # Mirror run_jumble(): raw .jumble.csv targets + process_snps() SNP table.
    # A separate copy is used so frankenplot() reads the CSV unmodified.
    gis_targets <- data.table::fread(opt$`jumble-csv`)
    gis_targets[, chromosome := sub("^chr", "", as.character(chromosome))]
    snp_table <- Jumble:::process_snps(opt$`tumor-snp-vcf`, gis_targets,
                                       sample_name = opt$`tumor-sample-name`)
    if (is.null(snp_table) || nrow(snp_table) == 0) {
        message("No usable SNPs in ", opt$`tumor-snp-vcf`, "; skipping GIS computation.")
    } else {
        gis_res <- Jumble:::compute_gis_and_maf(
            gis_targets, snp_table,
            reference = list(genome = opt$genome),
            hrd_model = NULL
        )
        hrdtable <- gis_res$gis_table
        if (!is.null(hrdtable)) {
            data.table::fwrite(hrdtable, gis_output, sep = ",")
            gis_written <- TRUE
        }
    }
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
} else if (compute_gis) {
    cat(paste0("HRD/GIS table: computed from ", opt$`tumor-snp-vcf`,
        " (sample ", if (is.null(opt$`tumor-sample-name`)) "<single-sample VCF>" else opt$`tumor-sample-name`,
        ", genome ", opt$genome, ")\n"))
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
cat("Output:", opt$output, "\n")
cat("Static PNG:", if (is.null(opt$`output-png`)) "auto" else opt$`output-png`, "\n\n")

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
    dpyd_csv = opt$`dpyd-csv`,
    output_png = opt$`output-png`
)

cat("\n✓ Frankenplot report generated\n")
cat("  - File:", opt$output, "\n")
cat("  - Static PNG:", if (is.null(opt$`output-png`)) sub("\\.html?$", ".png", opt$output, ignore.case = TRUE) else opt$`output-png`, "\n")
if (gis_written) {
    cat("  - GIS CSV:", gis_output, "\n")
}
cat("\nDone.\n")
