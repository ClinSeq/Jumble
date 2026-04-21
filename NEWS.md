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
