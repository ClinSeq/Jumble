# Jumble

**Copy Number Analysis of Short Read Sequencing Data**

Jumble is an R package for copy number analysis, offering functions for counting reads, building reference panels, and calling copy number alterations (CNA). It supports both Targeted Sequencing (gene panels) and Whole Genome Sequencing (WGS).

*   **[NEWS](NEWS.md)**: Changelog and version history.


# PROTOTYPE

This package is a re-implementation of the legacy *Jumble* R-scripts. It is currently in a prototype state and may not be stable. It may also be incomplete. Use at your own risk.



## Features

*   **Copy Number Calling**: Segmented copy number plotting and calling.
*   **Genomic Instability Score (GIS)**: Calculation and visualization (`.gis.png`).
*   **Pipeline Ready**: Includes command-line wrapper scripts for easy integration into Nextflow/Snakemake pipelines.


## Installation

You can install the development version of Jumble from GitHub (once hosted) or locally:

```r
# Install from local source
devtools::install(".")
```

## Quick Start (R Interactive)

Jumble comes with small dataset examples to get you started.

```r
library(Jumble)

# 1. Locate test data
testdata_dir <- system.file("testdata", package = "Jumble")
ref_file <- file.path(testdata_dir, "gene_panel/reference.RDS")
sample_file <- file.path(testdata_dir, "gene_panel/samples/test_sample_1.counts.RDS")
vcf_file <- file.path(testdata_dir, "gene_panel/samples/test_sample_1.vcf.gz")

# 2. Run Jumble
# Generates plots and tables in the output directory
results <- run_jumble(
  bam_file = sample_file,      # Can be BAM or pre-computed counts
  reference_file = ref_file,   # Reference panel
  snp_vcf = vcf_file,          # Optional: VCF for allele-specific visualization
  output_dir = "jumble_output"
)

# 3. Explore Results
print(head(results$segments))
```

