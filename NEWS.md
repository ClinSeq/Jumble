# Jumble 0.5.3

*   **Segmentation Engine**: Switched from PSCBS to DNAcopy CBS with `smooth.CNA`
    preprocessing. The smoothing step removes single-bin outliers before
    segmentation, improving detection of focal deletions. Added `trim` parameter
    (default 0.05) for robust segment mean estimation.
*   **Fragment Length Filtering**: New `exclude_long_fragments` parameter in
    `run_jumble()` and CLI. When TRUE, uses `count_medium` (fragments ≤300 bp)
    instead of all fragments as the main depth signal — designed for clipoverlap
    BAMs where TLEN inflation can misplace long-fragment midpoints.
*   **Methods Documentation**: Added comprehensive algorithmic methods document
    (`docs/METHODS.md`) covering the full normalization pipeline, GIS/HRD feature
    computation, MSI classification, and TMB estimation. Enhanced roxygen
    documentation for key functions (`correct_by_optim`, `gis_model`,
    `compute_gis_table`, `comp_gis_for_fraction`, `compute_gis_and_maf`).
*   **README**: Added Methods Overview section, citation block, and link to
    METHODS.md.
*   **DESCRIPTION**: Added URL field with GitHub repository and preprint DOI.

# Jumble 0.5.2

*   **Leave-Me-Out Guard**: Fixed crash when query sample matches all reference
    samples (single-sample reference). Exclusion is now skipped with a warning
    rather than leaving an empty reference.
*   **GC Correction**: `correct_by_gc` now trains loess on backbone-weighted bins
    only, and caps training set at 10,000 points.
*   **QC Output**: `compute_qc_metrics` returns 21 columns — added `TMB_snv`,
    `TMB_indel`, `TMB_score` (mutations per Mb; `NA` without somatic VCF).
*   **Plot Title**: Fixed scalar extraction from QC metric columns.

# Jumble 0.5.1

*   **Leave-Me-Out Restored**: Test samples present in the reference are now
    correctly excluded before PCA normalization, preventing self-normalization bias.
*   **Legacy Count Sanitization**: Restored the full two-step `sanitize_legacy_counts`
    function (exact-duplicate removal + phantom grid fingerprinting via target overlap).
    Replaces a broken interim implementation.
*   **Somatic VCF Fixes**:
    *   Restored VEP fallback parser for VCFs lacking CSQ annotation (parses INFO/GENE,
        INFO/EFFECT, INFO/CLINVAR).
    *   Restored MAX_AF population frequency filter to remove common germline
        variants from somatic calls.
    *   Restored chromosome name harmonization (`chr` prefix handling) in
        `map_variants_to_bins`.
*   **Local SNP Background**: Restored segment-based `local_snp_bg` smoothing
    (within-segment `runmed`, cross-segment extrapolation) for somatic filtering.
*   **Normalize Fix**: Restored `clean_chrom_names` call in `compute_reference_pca`.
*   **Strand Filter**: Fixed `gaps()` in `create_background_bins` to filter for
    `strand == "*"`, preventing phantom stranded background bins in new count files.

# Jumble 0.5.0

*   **Custom HRD Model**: Support for user-supplied HRD models (randomForest, glm,
    or plain function) via `hrd_model_file` parameter and `-m` CLI flag. Adds a
    `custom_HRD` column to GIS output.
*   **Segment Annotation**: Enhanced segment annotation with all protein-coding genes
    (not just cancer genes) and cytoband labeling. Gene lists are rolled up to
    "N genes" when a segment spans more than 10 genes.
*   **Gene Table Integration**: Full Ensembl gene table cached in reference for
    comprehensive segment-level gene annotation.
*   **Ideogram Improvements**: Improved ideogram visualization.

# Jumble 0.4.1

*   **MSI Classification**: Added microsatellite instability calling from somatic VCFs.
    *   New `classify_msi()` engine detects indels in mono-, di-, and trinucleotide repeat tracts.
    *   MSI VAF plot shows repeat category counts across VAF thresholds (log scale).
*   **QC Table Standardization**: Fixed 18-column output, always present (NA when input unavailable).
    *   Added `het_snps`, `hom_snps`, `sex`, `somatic_vcf` columns.
    *   Renamed `total_snvs`/`total_indels` → `somatic_snvs`/`somatic_indels`.
    *   Removed internal diagnostic columns (GC bin counts, target/background counts).
*   **Sex Inference**: Inferred from chrX heterozygosity vs autosomal rate, excluding pseudoautosomal regions.
*   **Somatic Plot Shapes**: Indels shown as filled triangles (▲ insertion, ▼ deletion), SNVs as circles.
*   **Bug fix**: Fixed `generate_counts` handling of pre-computed counts files.

# Jumble 0.4.0

*   **Preview Release**: First preview release of the Jumble R package for copy number analysis of short read sequencing data.
*   **Features**:
    *   Gene panel and WGS support.
    *   Somatic mutation integration and visualization.
    *   Genomic Instability Score (GIS) calculation for HRD estimation.
    *   Normalization and segmentation pipeline.
