# Jumble

**Copy Number Analysis of Short Read Sequencing Data**

Jumble is an R package for copy number analysis, offering functions for counting reads, building reference panels, and calling copy number alterations (CNA). It supports both Targeted Sequencing (gene panels) and Whole Genome Sequencing (WGS).

*   **[NEWS](NEWS.md)**: Changelog and version history.
*   **[Methods](docs/METHODS.md)**: Detailed algorithmic methods documentation.

> **Citation**: Mayrhofer, M. et al. "Sensitive detection of copy number alterations in samples with low circulating tumor DNA fraction." *MedRxiv* (2024). [doi:10.1101/2024.05.04.24306860](https://www.medrxiv.org/content/10.1101/2024.05.04.24306860v1)


## Features

*   **Copy Number Calling**: Segmented copy number plotting and calling.
*   **Genomic Instability Score (GIS)**: Calculation and visualization (`.gis.png`).
*   **Contamination Estimation**: Automatic DNA contamination detection from germline SNP allele ratios using a Random Forest model.
*   **MSI Classification**: Microsatellite instability calling from somatic VCFs.
*   **TMB Estimation**: Tumor mutational burden with rare germline SNP rejection.
*   **Custom HRD Models**: Support for user-supplied HRD classification models.
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

### 4. Generate Frankenplot Report (Optional)
Create an interactive HTML genome report from Jumble output files. GIS/HRD curves
are displayed when a precomputed Jumble GIS table is supplied with `hrdtable`.

```r
gis_table <- data.table::fread("jumble_results/tumor_sample.jumble_gis.csv")

frankenplot(
  tumor_jumble_csv = "jumble_results/tumor_sample.jumble.csv",
  tumor_cns = "jumble_results/tumor_sample.cns",
  output_file = "tumor_sample_frankenplot.html",
  tumor_snp_vcf = "germline_snps.vcf.gz",       # Optional: SNP allele ratios
  somatic_vcf = "somatic_mutations.vcf.gz",      # Optional: somatic overlay
  germline_vcf = "germline_mutations.vcf.gz",    # Optional: germline overlay
  hrdtable = gis_table                           # Optional: GIS/HRD curves
)
```

### Command Line Interface
Jumble also includes wrapper scripts for pipeline integration:
*   `jumble-count.R`: Generate counts from BAM.
*   `jumble-reference.R`: Build a reference panel.
*   `jumble-run.R`: Run the full analysis.
*   `jumble-frankenplot.R`: Generate a Frankenplot HTML report from Jumble output files.

Find script locations with: `system.file("scripts", package = "Jumble")`

Example Frankenplot wrapper call:

```sh
Rscript $(Rscript -e 'cat(system.file("scripts", "jumble-frankenplot.R", package = "Jumble"))') \
  --jumble-csv jumble_results/tumor_sample.jumble.csv \
  --cns jumble_results/tumor_sample.cns \
  --output tumor_sample_frankenplot.html \
  --tumor-snp-vcf germline_snps.vcf.gz \
  --somatic-vcf somatic_mutations.vcf.gz \
  --hrdtable jumble_results/tumor_sample.jumble_gis.csv
```

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
