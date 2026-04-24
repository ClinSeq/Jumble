# Command-Line Usage

The Jumble package provides command-line wrapper scripts that match the original workflow, allowing use in pipelines and production environments.

## Installation

After installing the package, the scripts are available in the package installation:

```r
# Find script location
system.file("scripts", package = "Jumble")
```

Or add to your PATH:

```bash
# Add to ~/.bashrc or ~/.zshrc
export PATH="$PATH:$(Rscript -e 'cat(system.file("scripts", package="Jumble"))')"
```

## Scripts

### 1. jumble-count.R

Generate read count files from BAM files.

**Single BAM:**
```bash
Rscript jumble-count.R -t targets.bed -b sample.bam
```

**All BAMs in directory (parallel):**
```bash
Rscript jumble-count.R -t targets.bed -c 4
```

**WGS mode:**
```bash
Rscript jumble-count.R -b sample.bam -w 10000
```

**Options:**
- `-t, --target FILE`: Target BED file (for targeted sequencing)
- `-b, --bam FILE`: Input BAM file (sorted and indexed)
- `-c, --cores N`: Number of cores for parallel processing [default: 1]
- `-w, --wgs-bin-size N`: Bin size for WGS mode [default: 10000]
- `-o, --output DIR`: Output directory [default: current directory]

**Output:** `<input_bam>.counts.RDS`

---

### 2. jumble-reference.R

Build a reference panel from count files.

**Usage:**
```bash
Rscript jumble-reference.R -i reference_counts/ -a biomart -o reference.RDS -c 4
```

**Options:**
- `-i, --input DIR`: Directory containing .counts.RDS files (required)
- `-a, --annotation SOURCE`: Annotation source: 'biomart' or path to folder [default: biomart]
- `-g, --genome GENOME`: Genome version: 'hg19' or 'hg38' [default: hg19]
- `-o, --output FILE`: Output file path [default: auto-generated]
- `-c, --cores N`: Number of cores [default: 1]

**Output:** `<target_bed>.reference.RDS`

**Notes:**
- Requires 4+ reference samples (10-100 recommended)
- Use "normal-like" samples without CNV alterations
- All samples must use the same BED file

---

### 3. jumble-run.R

Run Jumble copy number analysis.

**Basic usage:**
```bash
Rscript jumble-run.R -r reference.RDS -b sample.bam -o output/
```

**With SNP VCF:**
```bash
Rscript jumble-run.R -r reference.RDS -b sample.counts.RDS -v sample.vcf.gz
```

**With SNP VCF and Somatic VCF:**
```bash
Rscript jumble-run.R -r reference.RDS -b sample.bam -v germline.vcf.gz -s somatic.vcf.gz -o output/
```

**Options:**
- `-r, --reference FILE`: Reference file (required)
- `-b, --bam FILE`: Input BAM or .counts.RDS file (required)
- `-v, --vcf FILE`: Germline SNP VCF file for GIS/LOH analysis (optional)
- `-s, --somatic FILE`: Somatic VCF file for mutation overlay on plots (optional)
- `-o, --output DIR`: Output directory [default: current directory]
- `-a, --alpha N`: Segmentation alpha parameter [default: 0.001]

**Output files:**
- `<sample>.jumble.csv` - Normalized bin-level copy number data
- `<sample>.jumble_snps.csv` - SNP/LOH data (if SNP VCF provided)
- `<sample>.jumble_gis.csv` - GIS/HRD scores (if SNP VCF provided)
- `<sample>.qc.csv` - QC metrics (21 columns, always generated)
- `<sample>.png` - Copy number plot

---

## Complete Workflow Example

### Targeted Sequencing

```bash
# 1. Generate counts for all samples
Rscript jumble-count.R -t targets.bed -c 8

# 2. Build reference from normal samples
mkdir reference_counts
mv normal*.counts.RDS reference_counts/
Rscript jumble-reference.R -i reference_counts/ -o panel.reference.RDS -c 4

# 3. Run analysis on tumor samples
for sample in tumor*.counts.RDS; do
    vcf="${sample%.counts.RDS}.vcf.gz"
    Rscript jumble-run.R -r panel.reference.RDS -b "$sample" -v "$vcf" -o results/
done
```

### Whole Genome Sequencing

```bash
# 1. Generate counts (WGS mode)
Rscript jumble-count.R -w 10000 -c 8

# 2. Build reference
mkdir wgs_reference_counts
mv normal*.counts.RDS wgs_reference_counts/
Rscript jumble-reference.R -i wgs_reference_counts/ -o wgs.reference.RDS -c 4

# 3. Run analysis
Rscript jumble-run.R -r wgs.reference.RDS -b tumor.counts.RDS -o results/
```

---

## Comparison: Scripts vs R Package

| Feature | Command-Line Scripts | R Package Functions |
|---------|---------------------|---------------------|
| **Use Case** | Pipelines, production | Interactive analysis, development |
| **Installation** | Included in package | `library(Jumble)` |
| **Syntax** | `Rscript jumble-run.R -r ref.RDS -b sample.bam` | `run_jumble(bam_file, reference_file)` |
| **Options** | Command-line flags | Function arguments |
| **Parallelization** | `-c` flag | `cores` argument |
| **Output** | Files to disk | R objects + files |
| **Integration** | Shell scripts, workflows | R scripts, notebooks |

---

## Pipeline Integration

### Nextflow Example

```groovy
process jumble_count {
    input:
    path bam
    path bed
    
    output:
    path "${bam}.counts.RDS"
    
    script:
    """
    Rscript jumble-count.R -t ${bed} -b ${bam}
    """
}

process jumble_run {
    input:
    path counts
    path reference
    path vcf
    
    output:
    path "*.png"
    path "*.RDS"
    
    script:
    """
    Rscript jumble-run.R -r ${reference} -b ${counts} -v ${vcf}
    """
}
```

### Snakemake Example

```python
rule jumble_count:
    input:
        bam = "{sample}.bam",
        bed = "targets.bed"
    output:
        "{sample}.bam.counts.RDS"
    shell:
        "Rscript jumble-count.R -t {input.bed} -b {input.bam}"

rule jumble_run:
    input:
        counts = "{sample}.counts.RDS",
        reference = "reference.RDS",
        vcf = "{sample}.vcf.gz"
    output:
        "{sample}.jumble.RDS",
        "{sample}.png"
    shell:
        "Rscript jumble-run.R -r {input.reference} -b {input.counts} -v {input.vcf}"
```

---

## Notes

- All scripts require the Jumble package to be installed
- BAM files must be sorted and indexed (.bai)
- VCF files should be bgzipped and indexed (.tbi) for best performance
- Use `--help` flag for detailed options: `Rscript jumble-run.R --help`
