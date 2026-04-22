#' Compute Quality Control Metrics
#'
#' Calculates QC metrics for a Jumble analysis run. Always returns a fixed set
#' of 18 columns, filling with NA when the required input is unavailable.
#'
#' @param targets Targets data.table with bin information and counts.
#' @param bam_file Path to input BAM or counts file.
#' @param reference_file Path to reference file.
#' @param snp_vcf Path to SNP VCF file (optional).
#' @param somatic_vcf Path to somatic VCF file (optional).
#' @param sample_name Sample name.
#' @param contamination Estimated contamination fraction (optional).
#' @param snp_table SNP data.table from process_snps (optional, for het/hom/sex).
#' @param somatic Somatic variants data.table (optional, for SNV/indel/MSI counts).
#' @return A data.table with one row containing 18 QC metric columns.
#' @importFrom data.table data.table
#' @keywords internal
compute_qc_metrics <- function(targets, bam_file, reference_file,
                               snp_vcf = NULL, somatic_vcf = NULL,
                               sample_name = NULL,
                               contamination = NA_real_,
                               snp_table = NULL,
                               somatic = NULL) {
  # ── 1. Identity & Input Files ──────────────────────────────────────────────
  sample_name <- if (!is.null(sample_name)) sample_name else NA_character_
  qc <- data.table(
    sample         = sample_name,
    bam_file       = if (!is.null(bam_file)) basename(bam_file) else NA_character_,
    reference_file = if (!is.null(reference_file)) basename(reference_file) else NA_character_,
    snp_vcf        = if (!is.null(snp_vcf)) basename(snp_vcf) else NA_character_,
    somatic_vcf    = if (!is.null(somatic_vcf)) basename(somatic_vcf) else NA_character_
  )

  # ── 2. Technical QC (always computable from counts) ────────────────────────
  target_bins <- targets[is_target == TRUE]

  # Median target count
  qc$median_target_count <- if (nrow(target_bins) > 0) {
    median(target_bins$count, na.rm = TRUE)
  } else {
    NA_real_
  }

  # GC Bias
  qc$gc_bias <- compute_gc_bias(targets)

  # Noise (Linear MAPD)
  qc$noise <- compute_noise(target_bins)

  # Waviness
  qc$waviness <- compute_waviness(target_bins)

  # ── 3. Computed Estimates (require VCF inputs) ─────────────────────────────
  # Het/hom SNPs and sex (require snp_table + targets)
  snp_stats <- compute_snp_stats(snp_table, targets)
  qc$het_snps      <- snp_stats$het_snps
  qc$hom_snps      <- snp_stats$hom_snps
  qc$sex           <- snp_stats$sex
  qc$contamination <- contamination

  # Somatic counts (require somatic table)
  qc$somatic_snvs   <- NA_integer_
  qc$somatic_indels  <- NA_integer_
  qc$MSI_mono       <- NA_integer_
  qc$MSI_di         <- NA_integer_
  qc$MSI_tri        <- NA_integer_

  qc$TMB_snv        <- NA_integer_
  qc$TMB_indel      <- NA_integer_
  qc$TMB_score      <- NA_character_

  if (!is.null(somatic) && nrow(somatic) > 0) {
    # 1. Base MSI Filtering (AF >= 0.05)
    high_af_somatic <- somatic[AF >= 0.05]
    
    is_indel <- nchar(high_af_somatic$REF) != nchar(high_af_somatic$ALT)
    qc$somatic_snvs  <- sum(!is_indel)
    qc$somatic_indels <- sum(is_indel)

    if ("MSI" %in% names(high_af_somatic)) {
      msi_vals <- high_af_somatic$MSI[!is.na(high_af_somatic$MSI)]
      qc$MSI_mono <- sum(msi_vals == 1)
      qc$MSI_di   <- sum(msi_vals == 2)
      qc$MSI_tri  <- sum(msi_vals == 3)
    }

    # 2. TMB Footprint Masking
    if ("bin" %in% names(high_af_somatic)) {
      median_count <- median(targets$count[targets$is_target == TRUE], na.rm = TRUE)
      if (is.finite(median_count) && !is.na(median_count)) {
        valid_bin_indices <- which(targets$is_target == TRUE & targets$count > (0.2 * median_count) & targets$count > 50)
        
        if (length(valid_bin_indices) > 0) {
          target_mb <- sum(targets[valid_bin_indices, end - start], na.rm = TRUE) / 1e6
          
          # Reject background rare germline SNPs (LOH tracked)
          if ("is_rare_snp" %in% names(high_af_somatic)) {
            tmb_pool <- high_af_somatic[!is.na(bin) & bin %in% valid_bin_indices & (is_rare_snp == FALSE | is.na(is_rare_snp))]
          } else {
            tmb_pool <- high_af_somatic[!is.na(bin) & bin %in% valid_bin_indices]
          }
          
          tmb_is_indel <- nchar(tmb_pool$REF) != nchar(tmb_pool$ALT)
          
          tmb_indel <- sum(tmb_is_indel)
          tmb_snv <- sum(!tmb_is_indel)
          
          qc$TMB_snv <- tmb_snv
          qc$TMB_indel <- tmb_indel
          
          n_estim <- round(tmb_indel + (tmb_snv * 0.70))
          ptest <- stats::poisson.test(n_estim)
          est <- round(n_estim / target_mb, 1)
          l_ci <- round(ptest$conf.int[1] / target_mb, 1)
          u_ci <- round(ptest$conf.int[2] / target_mb, 1)
          
          qc$TMB_score <- sprintf("%.1f (%.1f-%.1f)", est, l_ci, u_ci)
        }
      }
    }
  }

  return(qc)
}


#' Compute GC Bias Metric
#'
#' @param targets Targets data.table
#' @return Numeric GC bias value (log2 ratio) or NA
#' @keywords internal
compute_gc_bias <- function(targets) {
  gc_bins <- targets[is_target == TRUE & !is.na(gc)]
  if (nrow(gc_bins) == 0) return(NA_real_)

  low_gc_bins  <- gc_bins[gc >= 0.3 & gc < 0.4]
  high_gc_bins <- gc_bins[gc >= 0.5 & gc < 0.6]

  mean_low_gc  <- mean(low_gc_bins$count, na.rm = TRUE)
  mean_high_gc <- mean(high_gc_bins$count, na.rm = TRUE)

  if (is.finite(mean_low_gc) && is.finite(mean_high_gc) &&
    mean_low_gc > 0 && mean_high_gc > 0) {
    return(round(log2(mean_high_gc / mean_low_gc), 2))
  }

  NA_real_
}


#' Compute Noise (Linear MAPD)
#'
#' @param target_bins Target bins data.table
#' @return Numeric noise value or NA
#' @keywords internal
compute_noise <- function(target_bins) {
  if (nrow(target_bins) < 2 || !"log2" %in% names(target_bins)) return(NA_real_)

  valid_log2 <- target_bins$log2[is.finite(target_bins$log2)]
  if (length(valid_log2) < 2) return(NA_real_)

  mapd <- median(abs(diff(valid_log2)))
  round((2^mapd) - 1, 2)
}


#' Compute Waviness (1Mb Window Smoothed SD)
#'
#' @param target_bins Target bins data.table
#' @return Numeric waviness value or NA
#' @keywords internal
compute_waviness <- function(target_bins) {
  if (nrow(target_bins) == 0 || !"log2" %in% names(target_bins)) return(NA_real_)

  w_dt <- target_bins[is.finite(log2), .(chromosome, start, log2)]
  if (nrow(w_dt) == 0) return(NA_real_)

  chrom_levels <- c(as.character(1:22), "X", "Y")
  w_dt[, sort_chr := stringr::str_remove(as.character(chromosome), "^chr")]
  w_dt[, sort_fac := factor(sort_chr, levels = chrom_levels)]
  data.table::setorder(w_dt, sort_fac, start, na.last = TRUE)

  if (nrow(w_dt) >= 11) {
    w_dt[, smoothed_log2 := stats::runmed(log2, k = 11, na.action = "na.omit")]
  } else {
    w_dt[, smoothed_log2 := log2]
  }

  w_dt[, window_seq := floor(start / 1e6)]
  w_dt[, window_id := paste0(sort_chr, "_", window_seq)]

  window_stats <- w_dt[, .(
    N = sum(!is.na(smoothed_log2)),
    sd_smoothed = {
      if (sum(!is.na(smoothed_log2)) >= 5) {
        stats::sd(smoothed_log2, na.rm = TRUE)
      } else {
        NA_real_
      }
    }
  ), by = window_id]

  valid_windows <- window_stats[!is.na(sd_smoothed) & N >= 5]
  if (nrow(valid_windows) == 0) return(NA_real_)

  round(weighted_median(valid_windows$sd_smoothed, valid_windows$N), 2)
}


#' Compute Het/Hom SNP Counts and Infer Sex
#'
#' @param snp_table SNP data.table (NULL if no SNP VCF)
#' @param targets Targets data.table (for counting X target bins)
#' @return List with het_snps, hom_snps, sex
#' @keywords internal
compute_snp_stats <- function(snp_table, targets) {
  result <- list(
    het_snps = NA_integer_,
    hom_snps = NA_integer_,
    sex      = NA_character_
  )

  if (is.null(snp_table) || nrow(snp_table) == 0) return(result)

  # Het: minor allele depth >= 2 AND minor allele ratio >= 1%
  minor_depth <- pmin(snp_table$AD, snp_table$RD)
  minor_ratio <- minor_depth / snp_table$DP
  is_het <- minor_depth >= 2 & minor_ratio >= 0.01
  is_het[is.na(is_het)] <- FALSE

  result$het_snps <- sum(is_het)
  result$hom_snps <- sum(!is_het)

  # Sex inference: compare non-PAR chrX het density to autosomal het density
  # Require >= 100 target bins on chrX to attempt
  chrom_clean <- stringr::str_remove(as.character(targets$chromosome), "^chr")

  # Pseudoautosomal regions (excluded — they behave like autosomes)
  # PAR1 hg19: chrX:60001-2699520,    PAR2 hg19: chrX:154931044-155260560
  # PAR1 hg38: chrX:10001-2781479,    PAR2 hg38: chrX:155701383-156030895
  # Use the union of both genome builds for robustness
  par_start <- c(10001, 154931044)
  par_end   <- c(2781479, 156030895)

  is_x_target <- chrom_clean == "X" & targets$is_target == TRUE
  is_par_target <- is_x_target &
    ((targets$start >= par_start[1] & targets$end <= par_end[1]) |
     (targets$start >= par_start[2] & targets$end <= par_end[2]))
  is_nonpar_x_target <- is_x_target & !is_par_target

  n_x_targets <- sum(is_nonpar_x_target)
  if (n_x_targets < 100) return(result)

  snp_chrom <- stringr::str_remove(as.character(snp_table$chromosome), "^chr")
  autosomes <- as.character(1:22)

  # Exclude PAR SNPs from X het count
  is_par_snp <- snp_chrom == "X" &
    ((snp_table$start >= par_start[1] & snp_table$start <= par_end[1]) |
     (snp_table$start >= par_start[2] & snp_table$start <= par_end[2]))

  het_x    <- sum(is_het & snp_chrom == "X" & !is_par_snp)
  het_auto <- sum(is_het & snp_chrom %in% autosomes)

  n_auto_targets <- sum(chrom_clean %in% autosomes & targets$is_target == TRUE)

  if (n_auto_targets == 0 || het_auto == 0) return(result)

  het_per_target_x    <- het_x / n_x_targets
  het_per_target_auto <- het_auto / n_auto_targets
  ratio <- het_per_target_x / het_per_target_auto

  if (ratio < 0.01) {
    result$sex <- "male"
  } else if (ratio > 0.05) {
    result$sex <- "female"
  }
  # else NA (ambiguous)

  return(result)
}


#' Write QC Metrics to CSV
#'
#' Writes QC metrics to a CSV file.
#'
#' @param qc_metrics QC metrics data.table.
#' @param output_file Path to output CSV file.
#' @importFrom data.table fwrite
#' @keywords internal
write_qc_metrics <- function(qc_metrics, output_file) {
  fwrite(qc_metrics, output_file, sep = ",", quote = TRUE)
  message("QC metrics saved to: ", output_file)
}
