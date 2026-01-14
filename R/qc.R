#' Compute Quality Control Metrics
#'
#' Calculates QC metrics for a Jumble analysis run.
#'
#' @param targets Targets data.table with bin information and counts.
#' @param bam_file Path to input BAM or counts file.
#' @param reference_file Path to reference file.
#' @param snp_vcf Path to SNP VCF file (optional).
#' @param sample_name Sample name.
#' @param contamination Estimated contamination fraction (optional).
#' @return A data.table with one row containing QC metrics.
#' @importFrom data.table data.table
#' @export
compute_qc_metrics <- function(targets, bam_file, reference_file,
                               snp_vcf = NULL, sample_name = NULL,
                               contamination = NA_real_) {
  # 1. Initialize ------------------------------------------------------------
  qc <- data.table(
    sample = if (!is.null(sample_name)) sample_name else NA_character_,
    bam_file = if (!is.null(bam_file)) basename(bam_file) else NA_character_,
    reference_file = if (!is.null(reference_file)) {
      basename(reference_file)
    } else {
      NA_character_
    },
    vcf_file = if (!is.null(snp_vcf)) basename(snp_vcf) else NA_character_,
    contamination = contamination
  )

  # 2. Median Counts ---------------------------------------------------------
  target_bins <- targets[is_target == TRUE]
  if (nrow(target_bins) > 0) {
    qc$median_target_count <- median(target_bins$count, na.rm = TRUE)
  } else {
    qc$median_target_count <- NA_real_
  }

  # 3. GC Bias ---------------------------------------------------------------
  # This metric assesses whether read coverage is biased by GC content.
  # We compare the mean coverage of high-GC bins (50-60%) to low-GC bins
  # (30-40%).
  # Formula: log2(mean_high_gc / mean_low_gc)
  # A value near 0 indicates no bias. Positive values indicate high-GC bias,
  # negative values indicate low-GC bias.

  # Define GC ranges:
  # High-GC: 0.5 <= gc < 0.6
  # Low-GC:  0.3 <= gc < 0.4

  # Filter for target bins with valid GC content
  gc_bins <- targets[is_target == TRUE & !is.na(gc)]

  if (nrow(gc_bins) > 0) {
    # Define GC ranges
    low_gc_bins <- gc_bins[gc >= 0.3 & gc < 0.4]
    high_gc_bins <- gc_bins[gc >= 0.5 & gc < 0.6]

    # Calculate mean counts
    mean_low_gc <- mean(low_gc_bins$count, na.rm = TRUE)
    mean_high_gc <- mean(high_gc_bins$count, na.rm = TRUE)

    # Compute GC bias (log2 ratio)
    if (is.finite(mean_low_gc) && is.finite(mean_high_gc) &&
      mean_low_gc > 0 && mean_high_gc > 0) {
      qc$gc_bias <- log2(mean_high_gc / mean_low_gc)
    } else {
      qc$gc_bias <- NA_real_
    }

    # Additional metrics for context
    qc$n_low_gc_bins <- nrow(low_gc_bins)
    qc$n_high_gc_bins <- nrow(high_gc_bins)
    qc$mean_low_gc_count <- mean_low_gc
    qc$mean_high_gc_count <- mean_high_gc
  } else {
    qc$gc_bias <- NA_real_
    qc$n_low_gc_bins <- 0L
    qc$n_high_gc_bins <- 0L
    qc$mean_low_gc_count <- NA_real_
    qc$mean_high_gc_count <- NA_real_
  }

  # 4. Global Metrics --------------------------------------------------------
  qc$total_bins <- nrow(targets)
  qc$target_bins <- sum(targets$is_target)
  qc$background_bins <- sum(!targets$is_target)

  # Mean count across all targets
  qc$mean_target_count <- mean(target_bins$count, na.rm = TRUE)

  return(qc)
}

#' Write QC Metrics to CSV
#'
#' Writes QC metrics to a CSV file.
#'
#' @param qc_metrics QC metrics data.table.
#' @param output_file Path to output CSV file.
#' @importFrom data.table fwrite
#' @export
write_qc_metrics <- function(qc_metrics, output_file) {
  fwrite(qc_metrics, output_file, sep = ",", quote = TRUE)
  message("QC metrics saved to: ", output_file)
}
