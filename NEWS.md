# Jumble 0.2.0

*   **Improvements**:
    *   Transitioned normalization to L1 loss and Total Variation (TV) penalty.
    *   Refined backbone estimation with downweighting.
    *   Added robust QC metrics ("Noise" and "Waviness") and integrated them into output plots.

# Jumble 0.1.0

*   **Initial Release**: First version of the Jumble package for copy number analysis of short read sequencing data.
*   **Features**:
    *   Gene panel and WGS support.
    *   Somatic mutation integration and visualization.
    *   Genomic Instability Score (GIS) calculation for HRD estimation.
    *   Normalization and segmentation pipeline.
*   **Note**: The contamination estimation feature (`estimate_contamination`) is currently disabled (returns `NA`) pending further refinement.
