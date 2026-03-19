# Jumble

**Copy Number Analysis of Short Read Sequencing Data**

Jumble is an R package for copy number analysis, offering functions for counting reads, building reference panels, and calling copy number alterations (CNA). It supports both Targeted Sequencing (gene panels) and Whole Genome Sequencing (WGS).

*   **[NEWS](NEWS.md)**: Changelog and version history.


# Prototype

This package is a re-implementation of legacy Jumble R-scripts. It is currently in a prototype state and may not be stable. It may also be incomplete. Use at your own risk.



## Features

*   **Copy Number Calling**: Segmented copy number plotting and calling.
*   **Genomic Instability Score (GIS)**: Calculation and visualization (`.gis.png`).
*   **Pipeline Ready**: Includes command-line wrapper scripts for easy integration into Nextflow/Snakemake pipelines.


## Installation

```r
# Install from GitHub
devtools::install_github("ClinSeq/Jumble")

# Or install from a local clone
devtools::install("/path/to/Jumble")
```

## Usage Workflow

The standard Jumble workflow consists of three steps: generating binned counts, building a reference panel from normal samples, and running the analysis on tumor data.

### 1. Generate Count Files
Generate binned read counts from your BAM files. This produces `.counts.RDS` files.

**Targeted Sequencing (Gene Panel):**
```r
library(Jumble)
generate_counts(
  bam_file = "sample.bam",
  target_bed_file = "targets.bed",
  output_dir = "counts_output"
)
```

**Whole Genome Sequencing (WGS):**
```r
generate_counts(
  bam_file = "sample.bam",
  bin_size = 10000,
  output_dir = "counts_output"
)
```

### 2. Build Reference Panel
Create a reference panel using a set of "normal" (CNV-negative) count files.

```r
count_files <- list.files("normal_counts/", pattern = ".counts.RDS$", full.names = TRUE)

reference <- build_reference(
  count_files = count_files,
  annotation_source = "biomart", # Fetches gene data from Ensembl
  genome = "hg19"                 # Supports hg19, hg38
)

saveRDS(reference, "my_reference.RDS")
```

### 3. Run Jumble Analysis
Normalize and segment your tumor sample against the reference panel.

```r
results <- run_jumble(
  bam_file = "tumor_sample.bam", # Can also be a pre-computed .counts.RDS
  reference_file = "my_reference.RDS",
  snp_vcf = "germline_snps.vcf.gz", # Optional: for GIS/LOH calculation
  output_dir = "jumble_results"
)

# Access segments and QC metrics
print(head(results$segments))
print(results$qc)
```

### Command Line Interface
Jumble also includes wrapper scripts for pipeline integration:
*   `jumble-count.R`: Generate counts from BAM.
*   `jumble-reference.R`: Build a reference panel.
*   `jumble-run.R`: Run the full analysis.

Find script locations with: `system.file("scripts", package = "Jumble")`

## Output Files

`run_jumble()` generates the following files in `output_dir`:

| File | Format | Description |
|------|--------|-------------|
| `<sample>.jumble.csv` | CSV | Bin-level targets with log2 ratios and counts |
| `<sample>.jumble_snps.csv` | CSV | SNP allele ratios (requires SNP VCF) |
| `<sample>.jumble_gis.csv` | CSV | GIS arm-level scores (requires SNP VCF) |
| `<sample>.qc.csv` | CSV | QC metrics (18 columns, see [QC_METRICS.md](docs/QC_METRICS.md)) |
| `<sample>.cnr` | TSV | CNVkit-compatible bin-level file |
| `<sample>.cns` | TSV | CNVkit-compatible segment file |
| `<sample>.seg` | TSV | IGV/GISTIC segment file |
| `<sample>.jumble.png` | PNG | Copy number plot |
| `<sample>.msi.png` | PNG | MSI VAF plot (requires somatic VCF) |

## Quick Verification (Internal Test Data)
To verify your installation using the small dataset included in the package:

```r
testdata_dir <- system.file("testdata", package = "Jumble")
ref_file <- file.path(testdata_dir, "gene_panel/reference.RDS")
sample_file <- file.path(testdata_dir, "gene_panel/samples/test_sample_1.counts.RDS")

results <- run_jumble(sample_file, ref_file, output_dir = "test_run")
```

