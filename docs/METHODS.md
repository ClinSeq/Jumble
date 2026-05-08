# Jumble: Algorithmic Methods

> **Version**: Jumble 0.5.3
> **Author**: Markus Mayrhofer
> **Preprint**: Mayrhofer, M. et al. "Sensitive detection of copy number alterations in samples with low circulating tumor DNA fraction." *MedRxiv* (2024). [doi:10.1101/2024.05.04.24306860](https://www.medrxiv.org/content/10.1101/2024.05.04.24306860v1)

This document describes the algorithmic methods implemented in the Jumble R package for copy number analysis of short-read sequencing data. It covers the normalization pipeline, genomic instability scoring, microsatellite instability classification, and somatic variant integration.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Read Counting and Bin Construction](#2-read-counting-and-bin-construction)
3. [Reference Panel Construction](#3-reference-panel-construction)
4. [Normalization Pipeline](#4-normalization-pipeline)
5. [Segmentation](#5-segmentation)
6. [Genomic Instability Score (GIS)](#6-genomic-instability-score-gis)
7. [MSI Classification](#7-msi-classification)
8. [Somatic Variant Integration](#8-somatic-variant-integration)
9. [Quality Control Metrics](#9-quality-control-metrics)

---

## 1. Overview

Jumble is a copy number analysis method designed for both targeted sequencing (gene panels) and whole-genome sequencing (WGS) data. The core pipeline consists of:

1. **Counting** — Binned read counting from BAM files with fragment length stratification
2. **Reference construction** — Building a reference panel from normal samples with PCA decomposition
3. **Normalization** — L1+TV-penalized PCA correction and GC bias correction
4. **Segmentation** — Circular Binary Segmentation (CBS) via DNAcopy
5. **GIS scoring** — Purity-aware genomic instability feature extraction and scoring
6. **Variant integration** — MSI classification, TMB estimation, and somatic annotation

The method is implemented in R and uses Bioconductor infrastructure (GenomicRanges, bamsignals, DNAcopy, VariantAnnotation).

---

## 2. Read Counting and Bin Construction

### 2.1 Bin Definition

For **targeted sequencing**, bins are constructed from the input BED file:

1. Target regions are extended by 50 bp on each side and merged (reduced) to eliminate overlaps.
2. Merged regions are snapped to 200 bp bin boundaries (center-aligned) and split into fixed 200 bp bins.
3. **Background (antitarget) bins** are constructed by inverting the target regions: target regions are extended by 1 kb on each side, and the genomic gaps between them are split into 1 Mb bins (minimum 300 kb).

For **WGS**, the genome is tiled into fixed-size bins (default 10 kb).

### 2.2 Fragment Length Stratification

Jumble counts reads in three fragment length categories per bin:

| Count | Fragment Length | Purpose |
|-------|----------------|---------|
| `count` | All fragments | Primary signal (default) |
| `count_short` | ≤150 bp | Short-fragment signal; enriched for ctDNA in liquid biopsy |
| `count_medium` | ≤300 bp | Medium-fragment signal; excludes long fragments |

All counting uses **midpoint** mode (fragment midpoint determines bin assignment), MAPQ ≥ 20, and filters QC-fail, optical duplicate, supplementary, secondary, and unmapped reads (SAM flag 2816).

The `exclude_long_fragments` parameter switches the primary signal from `count` (all fragments) to `count_medium` (≤300 bp). This is designed for BAMs processed with `clipoverlap` (overlap clipping), where the TLEN field may be inflated by the clipping process, causing the midpoint calculation to misplace long fragments.

### 2.3 Tiled Region Detection

Bins are classified as "tiled" if the distance to the nearest neighboring bin is less than 250 bp. This identifies densely targeted regions (e.g., exon-by-exon tiling) which receive special handling during normalization.

---

## 3. Reference Panel Construction

### 3.1 Reference Composition

A reference panel is built from a set of CNV-negative ("normal") samples. The reference object stores:

- **Raw counts** for all reference samples (all three fragment length categories)
- **Target template** with bin coordinates, GC content, and gene annotations
- **Gene annotations** from Ensembl (protein-coding genes and exons)
- **Cancer gene list** from a bundled CSV with oncogene/TSG/ambiguous classification
- **Cytoband coordinates** for genomic region labeling

### 3.2 GC Content Calculation

Per-bin GC content is computed from the reference genome sequence (BSgenome) using the `letterFrequency` function from Biostrings. This is stored in the reference and used for GC bias correction during normalization.

### 3.3 Leave-Me-Out

When a query sample is found to be present in the reference panel (detected by exact count vector matching), it is automatically excluded before PCA computation to prevent self-normalization bias.

---

## 4. Normalization Pipeline

The normalization pipeline transforms raw read counts into corrected log2 ratios suitable for segmentation and copy number calling. The pipeline operates in the following stages:

### 4.1 Log Ratio Calculation

Raw counts are converted to log2 values with a floor of 1 (to avoid log(0)):

```
rawLR = log2(max(count, 1))
```

This is computed separately for the primary count and the short-fragment count.

### 4.2 Backbone-Weighted Median Correction

Each sample's log2 values are centered by subtracting the **weighted median** of the backbone bins. The backbone weighting scheme prevents densely targeted genes from dominating the centering:

| Condition | Weight |
|-----------|--------|
| Non-autosomal bin (X, Y) | 0.0 |
| Bin with no gene annotation | 1.0 |
| Bin in gene with ≤ 25 bins | 1.0 |
| Bin in gene with > 25 bins | 25 / n_bins |

This ensures that a gene covered by 200 bins contributes the same total weight as a gene covered by 25 bins, preventing large genes from biasing the median.

### 4.3 Sex Chromosome Correction

Sex chromosome correction is applied independently for X and Y:

1. **Gender detection**: Inferred from the median log2 of X chromosome bins (male if 2^median < 0.75) and Y chromosome bins (male if 2^median > 0.25), independently.
2. **X correction**: For males, non-pseudoautosomal X bins (2.70 Mb – 154.93 Mb) are shifted up by +1 log2 unit (compensating for hemizygosity).
3. **Y correction**: For males, Y bins below 28.79 Mb are shifted up by +1 log2 unit. For females, Y bins are set to NA.

### 4.4 Reference Median Subtraction

The bin-wise median across all reference samples is subtracted from each bin's log2 value. This removes systematic bin-level biases (e.g., mappability, GC content) that are shared across samples.

### 4.5 PCA Decomposition

Principal Component Analysis is performed on the reference panel's log2 matrix (bins × samples) to identify latent technical factors. PCA is computed **separately** for:

- Target bins (on-target regions)
- Background bins (off-target regions)
- Primary counts and short-fragment counts (4 PCA decompositions total)

The PCA scores (loadings per bin) are stored and used to correct the query sample.

### 4.6 L1 + Total Variation PCA Correction

This is the core normalization innovation in Jumble. Rather than using standard linear regression to remove PCA components (which assumes Gaussian residuals and can over-correct focal CNAs), Jumble uses an **L1 norm + Total Variation (TV) penalty** optimization:

#### Objective Function

For a query sample with log2 ratio vector **y** and PCA score matrix **P**, the correction coefficients **c** are found by minimizing:

```
minimize  Σ|yᵢ - (Pc)ᵢ|  +  λ · Σ|Δ(y - Pc)ᵢ|
    c

where:
  yᵢ         = log2 ratio of bin i
  (Pc)ᵢ      = PCA-predicted correction for bin i
  Δ(·)ᵢ      = first difference operator: (·)ᵢ₊₁ - (·)ᵢ
  λ           = TV penalty ratio (default: 1.0)
```

The first term (L1 norm) is robust to outliers — focal copy number alterations appear as outliers relative to the PCA model and are preserved rather than corrected away. The second term (Total Variation) penalizes roughness in the corrected signal, encouraging smooth corrections that respect the spatial structure of the genome.

#### Progressive Coefficient Expansion

The optimization uses Nelder-Mead (a derivative-free simplex method) and builds coefficients progressively:

1. Start with the first 3 PCs and optimize their coefficients
2. Expand by 10 PCs at a time, using the previous solution as the starting point
3. Continue until all available PCs are included (up to 12 for the RLM fallback)

This progressive strategy avoids local minima that can trap the optimizer when starting with many parameters simultaneously.

#### Comparison with Robust Linear Model (RLM)

Jumble also supports a fallback correction method using MASS::rlm (robust linear regression with Huber weights). The `correction` parameter controls which method is used:

- `"optim"` (default): L1+TV optimization — better preservation of focal events
- `"rlm"`: Robust linear model — faster, adequate for samples with few focal events

### 4.7 GC Bias Correction

After PCA correction, residual GC bias is removed using **loess smoothing**:

1. A loess model is fitted: `log2 ~ gc`, using backbone-weighted bins as training data
2. The loess prediction is subtracted from all bins
3. Training is capped at 10,000 points for computational efficiency
4. The `family = "symmetric"` option provides robustness to outliers

### 4.8 X Chromosome Panel Correction

For targeted sequencing panels, an additional X chromosome correction compensates for systematic differences between on-target and off-target (background) bins on chrX. The correction computes the median difference between background and target log2 values in 1 Mb windows, smooths with a running median (k=11), and applies the median offset to target bins.

### 4.9 Value Clamping

Final log2 values are clamped to [-5, +7] for the primary signal and [-4, +7] for the short-fragment signal to prevent extreme outliers from distorting visualization and downstream analysis.

---

## 5. Segmentation

### 5.1 Circular Binary Segmentation

Jumble uses DNAcopy's Circular Binary Segmentation (CBS) algorithm with `smooth.CNA` preprocessing:

1. **smooth.CNA**: Removes single-bin outliers before segmentation, improving detection of focal deletions that might otherwise be masked by noisy neighbors.
2. **CBS parameters**:
   - α = 0.01 for targeted sequencing (permissive, suitable for sparse data)
   - α = 10⁻⁵ for WGS (stringent, appropriate for dense data)
   - Trim = 0.05 (5% trimmed mean for segment values)

### 5.2 Segment Processing

After CBS, segment means are recalculated as the **median** of constituent bin log2 values (more robust than the CBS trimmed mean for small segments). Segments are annotated with:

- Overlapping cancer genes (with oncogene/TSG classification)
- Exonic overlap status
- Cytoband labels
- Gene counts (rolled up to "N genes" when > 10)

---

## 6. Genomic Instability Score (GIS)

The GIS module estimates Homologous Recombination Deficiency (HRD) from copy number and allele ratio data. It is designed to work with both gene panel and WGS data.

### 6.1 Multi-Scale Binning

The GIS computation operates on a coarser resolution than the primary analysis:

1. **1 Mb binning**: Target-level data is aggregated into 1 Mb bins. Values are averaged by segment ID to use the segmented signal.

2. **5 Mb binning**: A genome-wide 5 Mb bin template is constructed from chromosome arm definitions (centromeric bins excluded). The 1 Mb values are mapped to 5 Mb bins via genomic overlap, taking the median per bin.

This multi-scale approach ensures that the GIS features capture arm-level and sub-arm-level events rather than bin-level noise, and makes the scoring applicable to both sparse panel data and dense WGS data.

### 6.2 Tumor Fraction Sweep

GIS features are computed across a sweep of 100 assumed tumor fractions (0.01 to 1.00 in steps of 0.01). This produces a GIS score profile as a function of purity, allowing the user to:

- Assess GIS at a known tumor fraction
- Evaluate the sensitivity of the GIS call to purity uncertainty

### 6.3 Purity-Aware Thresholds

All feature thresholds are scaled by the assumed tumor fraction:

```
threshold = max(tumor_fraction / 3, 0.05)
threshold_log2 = log2(1 + threshold)
```

This ensures that at low tumor fractions, the thresholds are relaxed to detect diluted signals, while at high fractions, the thresholds are tightened to avoid false positives from noise.

### 6.4 Feature Definitions

Nine genomic instability features are computed per tumor fraction. All features except `transitions` are counted as the **number of distinct chromosome arms** exhibiting the feature (max 46 arms for autosomes + X):

#### Copy Number Features

| Feature | Definition | Aggregation |
|---------|-----------|-------------|
| **transitions** | Number of 5 Mb bins where `\|Δlog2\| > threshold_log2` between adjacent bins | Sum of bins (both sides of each transition) |
| **long_cnv** | Arms where the 25 Mb running median deviates from the arm median by > threshold_log2 | Count of unique arms |
| **local_cnv** | Arms with any bin showing local gain OR local loss | Count of unique arms |
| **local_gain** | Arms with any bin where `log2 - arm_median > threshold_log2` | Count of unique arms |
| **local_loss** | Arms with any bin where `log2 - arm_median < -threshold_log2` | Count of unique arms |
| **focal_gain** | Arms with any bin where `log2 - long_median > threshold_log2` AND `log2 - arm_median > threshold_log2` | Count of unique arms |
| **focal_loss** | Arms with any bin where `log2 - long_median < -threshold_log2` AND `log2 - arm_median < -threshold_log2` | Count of unique arms |

The `long_median` is a 25 Mb running median (k=5 bins × 5 Mb) computed per chromosome, capturing large-scale trends. The `arm_median` is the median log2 across all bins on a chromosome arm.

**Focal** events are defined as deviations from BOTH the long-range trend and the arm-level baseline — this dual requirement distinguishes true focal events from arm-level shifts.

#### Allele-Based Features

| Feature | Definition | Aggregation |
|---------|-----------|-------------|
| **loh** | Arms with > 3 bins showing LOH AND > 3 bins without LOH (partial LOH) | Count of unique arms |
| **tai** | Arms with allelic imbalance at a telomeric bin | Count of unique arms |

### 6.5 Purity-Aware LOH Calling

LOH is detected using a purity-aware threshold that accounts for the expected allele ratio given the copy number state and tumor fraction:

```
Step 1: Estimate copy ratio from log2 and tumor fraction (tfr):
  copyratio = 1 + (2^log2 - 1) / tfr

Step 2: Compute expected DNA ratio (fraction of tumor DNA at this locus):
  dnaratio = copyratio × tfr / (copyratio × tfr + 1 × (1 - tfr))

Step 3: Compute LOH threshold:
  thr = dnaratio × 0.8
  thr_maf = 0.5 + thr × 0.5
  thr_maf = max(thr_maf, 0.55)

Step 4: A bin is LOH if MAF > thr_maf
```

The formula models the expected allele frequency shift under single-copy loss: at a given tumor fraction and copy number state, the expected MAF for a heterozygous SNP in a region of LOH is determined by the ratio of tumor to normal DNA. The 0.8 scaling factor provides a margin for noise, and the 0.55 floor prevents calling LOH in regions with minimal allelic imbalance.

LOH calls are smoothed with a running median (k=5 bins = 25 Mb) per arm to reduce noise.

### 6.6 Allelic Imbalance and TAI

Telomeric Allelic Imbalance (TAI) is computed as:

1. MAF values are cleaned (NA → 0.5) and smoothed with a running median (k=5) per arm
2. A bin has allelic imbalance if the smoothed MAF > 0.6
3. TAI counts the number of arms where a **telomeric** bin (first bin of a p-arm or last bin of a q-arm) shows allelic imbalance

### 6.7 Linear Prediction Model

The GIS score is predicted from three features using a hard-coded linear model:

```
GIS = 2.632 + (-0.852 × focal_gain) + (1.799 × local_cnv) + (2.727 × loh)
```

The coefficients were derived from a training set of samples with known HRD status (Myriad myChoice GIS scores). The model uses only three of the nine features:

- **loh** (coefficient 2.727): Strongest predictor — LOH is a hallmark of HRD
- **local_cnv** (coefficient 1.799): Captures arm-level copy number changes characteristic of HRD
- **focal_gain** (coefficient -0.852): Negative coefficient — focal gains are more characteristic of non-HRD tumors (e.g., oncogene amplification)

### 6.8 Custom HRD Model Interface

Users can supply their own HRD model via the `hrd_model_file` parameter. The model receives a data frame with all or some of the 9 GIS features and can be:

- A **classification model** (e.g., randomForest): Returns the positive-class probability × 100
- A **regression model** (e.g., glm): Returns the raw numeric prediction
- A **plain R function**: Receives the feature data frame and returns a numeric score

The interface auto-detects the model type by attempting `predict(model, type = "prob")` first, falling back to regression prediction. For classification models, the positive class is identified by column name priority: GIS > HRD > Pos > 1.

---

## 7. MSI Classification

### 7.1 Overview

Jumble classifies somatic indels as microsatellite instability (MSI)-like based on whether they occur within repeat tracts. The MSI score for a sample is the count of MSI-like indels.

### 7.2 Flanking Sequence Extraction

For each indel variant, left and right flanking sequences (25 bp each) are extracted from the reference genome. Processing is done per-chromosome for memory efficiency. Chromosome name styles are automatically harmonized between the VCF and the BSgenome reference.

### 7.3 Repeat Tract Classification Algorithm

For each indel, the algorithm tests repeat periods k ∈ {1, 2, 3} (shortest first — first match wins):

1. **Motif extraction**: The first k characters of the variant sequence (anchor base stripped) form the candidate motif.
2. **Dominance check**: The variant sequence must be > 80% composed of repeats of the motif.
3. **Bidirectional extension**: Count consecutive motif repeats extending backward into the left flank and forward into the right flank.
4. **Tract length**: `tract_len = left_repeats + variant_repeats + right_repeats` (in motif units).
5. **Threshold check**: If `tract_len ≥ threshold[k]`, classify as MSI period k.

### 7.4 Thresholds

| Period | Motif Example | Threshold (repeat units) | Threshold (bases) | Validation |
|--------|--------------|--------------------------|-------------------|------------|
| 1 (mono) | A, T, G, C | ≥ 6 | ≥ 6 bp | AUC = 0.930 (Youden-optimal) |
| 2 (di) | CA, AT, etc. | ≥ 4 | ≥ 8 bp | Provisional |
| 3 (tri) | CAG, etc. | ≥ 3 | ≥ 9 bp | Provisional |

Mononucleotide repeat indels (period 1) are indicative of MLH1/MSH2/PMS2 loss (classic MSI). Dinucleotide repeat indels (period 2) are indicative of MSH3 loss.

---

## 8. Somatic Variant Integration

### 8.1 Somatic VCF Processing

Somatic variants are extracted from VCF files with support for multiple caller formats:

- **Allele frequency**: Extracted from VAF, AF, or AD fields (in priority order)
- **VEP annotation**: Parsed from CSQ INFO field when available, with fallback to INFO/GENE, INFO/EFFECT, INFO/CLINVAR
- **Population frequency filtering**: Variants with MAX_AF > 0.01 (1% population frequency) are removed as likely germline

### 8.2 Tumor Mutational Burden (TMB) Estimation

TMB is computed as mutations per megabase of callable target region:

1. **Callable region**: Target bins with count > max(50, 0.2 × median_target_count) are considered callable.
2. **Variant filtering**: Only variants with AF ≥ 0.05 in callable bins are counted.
3. **Rare germline SNP rejection**: Variants flagged as `is_rare_snp` (see §8.3) are excluded from the TMB count.


### 8.3 Rare Germline SNP Rejection via LOH Integration

A key challenge in TMB estimation from panel data is distinguishing true somatic variants from rare germline SNPs (population frequency < 1%) that are not in common databases. Jumble addresses this by integrating LOH information:

1. **Local SNP background computation**: For each genomic segment, the median MAF of germline heterozygous SNPs is computed and smoothed with a running median (k=9 within segment).

2. **Cross-segment extrapolation**: Segments without germline SNP data inherit the background MAF from the nearest adjacent segment, weighted by log2 ratio similarity:
   ```
   For segment S without SNP data:
     - Find nearest left segment L and right segment R with SNP data
     - Compute logR distance: dist_left = |logR(S) - logR(L)|, dist_right = |logR(S) - logR(R)|
     - Assign background from the segment with smaller logR distance
   ```
   This logR-distance weighting ensures that segments at similar copy number states share allelic background, which is biologically appropriate (segments at the same copy number are likely to have similar LOH status).

3. **Flagging**: A somatic variant is flagged as `is_rare_snp = TRUE` if:
   ```
   |somatic_MAF - background_MAF| ≤ 0.05
   ```
   This identifies variants whose allele frequency is consistent with the local germline heterozygous background — i.e., they look like germline SNPs rather than somatic mutations.

---

## 9. Quality Control Metrics

Jumble computes 21 QC metrics per sample. The key algorithmic metrics are:

### 9.1 GC Bias

```
gc_bias = log2(mean_count[GC ∈ 0.5-0.6] / mean_count[GC ∈ 0.3-0.4])
```

Values near 0 indicate no GC bias. Positive values indicate high-GC enrichment.

### 9.2 Noise (Linear MAPD)

```
noise = 2^(median(|Δlog2|)) - 1
```

where Δlog2 is the first difference of adjacent bin log2 values. This is a linearized version of the Median Absolute Pairwise Difference (MAPD), expressed as a fraction of the signal.

### 9.3 Waviness

Large-scale systematic bias is quantified as:

1. Log2 values are smoothed with a running median (k=11)
2. Standard deviation is computed in 1 Mb windows (requiring ≥ 5 bins per window)
3. The weighted median of window SDs (weighted by bin count) is reported

### 9.4 Sex Inference

Sex is inferred from the ratio of chrX heterozygous SNP density to autosomal heterozygous SNP density, excluding pseudoautosomal regions (PAR1 and PAR2, using the union of hg19 and hg38 boundaries for robustness):

```
ratio = (het_SNPs_on_nonPAR_X / n_X_targets) / (het_SNPs_on_autosomes / n_auto_targets)

If ratio < 0.01: male
If ratio > 0.05: female
Otherwise: ambiguous (NA)
```

### 9.5 Contamination

The contamination estimation algorithm is currently retired (returns NA) pending further development.

---

## References

- Mayrhofer, M. et al. "Sensitive detection of copy number alterations in samples with low circulating tumor DNA fraction." *MedRxiv* (2024). [doi:10.1101/2024.05.04.24306860](https://www.medrxiv.org/content/10.1101/2024.05.04.24306860v1)
- Olshen, A.B. et al. "Circular binary segmentation for the analysis of array-based DNA copy number data." *Biostatistics* 5(4):557-572 (2004).
- Venkatraman, E.S. & Olshen, A.B. "A faster circular binary segmentation algorithm for the analysis of array CGH data." *Bioinformatics* 23(6):657-663 (2007).
