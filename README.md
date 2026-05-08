# Jumble

**Copy Number Analysis of Short Read Sequencing Data**

Jumble is an R package for copy number analysis, offering functions for counting reads, building reference panels, and calling copy number alterations (CNA). It supports both Targeted Sequencing (gene panels) and Whole Genome Sequencing (WGS).

*   **[NEWS](NEWS.md)**: Changelog and version history.
*   **[Methods](docs/METHODS.md)**: Detailed algorithmic methods documentation.

> **Citation**: Mayrhofer, M. et al. "Sensitive detection of copy number alterations in samples with low circulating tumor DNA fraction." *MedRxiv* (2024). [doi:10.1101/2024.05.04.24306860](https://www.medrxiv.org/content/10.1101/2024.05.04.24306860v1)


## Features

*   **Copy Number Calling**: Segmented copy number plotting and calling.
*   **Genomic Instability Score (GIS)**: Calculation and visualization (`.gis.png`).
*   **MSI Classification**: Microsatellite instability calling from somatic VCFs.
*   **TMB Estimation**: Tumor mutational burden with rare germline SNP rejection.
*   **Custom HRD Models**: Support for user-supplied HRD classification models.
*   **Pipeline Ready**: Includes command-line wrapper scripts for easy integration into Nextflow/Snakemake pipelines.


## Installation

You can install the development version of Jumble from GitHub (once hosted) or locally:

```r
# Install from local source
devtools::install(".")
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
  target_bed = "targets.bed"
)
```

**Whole Genome Sequencing (WGS):**
```r
generate_counts(
  bam_file = "sample.bam",
  wgs_bin_size = 10000
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
print(head(results$targets))
```

### Command Line Interface
Jumble also includes wrapper scripts for pipeline integration:
*   `jumble-count.R`: Generate counts from BAM.
*   `jumble-reference.R`: Build a reference panel.
*   `jumble-run.R`: Run the full analysis.

Find script locations with: `system.file("scripts", package = "Jumble")`

## Quick Verification (Internal Test Data)
To verify your installation using the small dataset included in the package:

```r
testdata_dir <- system.file("testdata", package = "Jumble")
ref_file <- file.path(testdata_dir, "gene_panel/reference.RDS")
sample_file <- file.path(testdata_dir, "gene_panel/samples/test_sample_1.counts.RDS")

results <- run_jumble(sample_file, ref_file, output_dir = "test_run")
```

