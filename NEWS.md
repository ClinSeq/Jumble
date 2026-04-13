# Jumble 0.5.0

*   **Y-Chromosome Artifact Sanitization**: Eliminated "double background bin" copy number artifacts caused by a bug in `generate_counts`. Existing `.counts.RDS` files are handled automatically and should not need to be regenerated.
*   **TMB (Tumor Mutational Burden)**: Added TMB estimation to QC output when a somatic VCF is provided.
    *   New QC columns: `TMB_snv`, `TMB_indel`, `TMB_score` (Poisson CI: e.g. `8.3 (6.1-11.0)` per Mb).
    *   TMB is only computed over well-covered target bins (count ≥ 20% of median and ≥ 50 reads).


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
