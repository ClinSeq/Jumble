# Package-level environment for caching the contamination RF model
.jumble_contam_env <- new.env(parent = emptyenv())

#' Load Contamination RF Model
#'
#' Loads the pre-trained Random Forest model from inst/extdata, caching it
#' in a package-level environment to avoid re-reading on repeated calls.
#'
#' @return The model list (model, feature_cols, bin_breaks_desc, etc.), or NULL on failure.
#' @keywords internal
load_contamination_model <- function() {
  if (!is.null(.jumble_contam_env$model)) {
    return(.jumble_contam_env$model)
  }

  model_path <- system.file("extdata", "contamination_rf.RDS", package = "Jumble")
  if (model_path == "" || !file.exists(model_path)) {
    return(NULL)
  }

  tryCatch({
    model_data <- readRDS(model_path)
    .jumble_contam_env$model <- model_data
    return(model_data)
  }, error = function(e) {
    warning("Failed to load contamination model: ", e$message)
    return(NULL)
  })
}

#' Estimate DNA Contamination
#'
#' Estimates the level of DNA contamination in a sample using a Random Forest
#' model trained on hom-alt SNP VAF histograms. The algorithm:
#'
#' \enumerate{
#'   \item Selects homozygous-alt SNPs (allele_ratio > 0.75, DP >= 20)
#'   \item Applies 5 Mb balanced-region selection (top 10\% windows by het count)
#'   \item Builds a 10-bin log-spaced VAF histogram (75-100\%)
#'   \item Predicts contamination fraction using a pre-trained Random Forest
#' }
#'
#' Returns \code{NA_real_} when:
#' \itemize{
#'   \item Input is NULL or empty
#'   \item Fewer than 50 hom-alt SNPs with DP >= 20 (e.g., hom SNPs filtered from VCF)
#'   \item The \pkg{randomForest} package is not installed
#'   \item The pre-trained model file is missing
#' }
#'
#' @param snp_table A data.table from \code{process_snps()}, containing at minimum:
#'   \code{chromosome}, \code{start}, \code{allele_ratio}, \code{DP}, \code{AD}, \code{RD}.
#' @return A numeric value in [0, 1] representing the estimated contamination fraction,
#'   or \code{NA_real_} if estimation is not possible.
#' @importFrom data.table data.table
#' @keywords internal
estimate_contamination <- function(snp_table) {
  # Guard: NULL or empty input

  if (is.null(snp_table) || !is.data.frame(snp_table) || nrow(snp_table) == 0) {
    return(NA_real_)
  }

  # Guard: required columns

  required_cols <- c("chromosome", "start", "allele_ratio", "DP", "AD", "RD")
  if (!all(required_cols %in% names(snp_table))) {
    return(NA_real_)
  }

  # Guard: randomForest package
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    warning("Package 'randomForest' is required for contamination estimation. ",
            "Install it with: install.packages('randomForest')")
    return(NA_real_)
  }

  # Guard: load model
  model_data <- load_contamination_model()
  if (is.null(model_data)) {
    warning("Contamination RF model not found in package extdata.")
    return(NA_real_)
  }

  # 1. Select hom-alt SNPs (allele_ratio > 0.75, DP >= 20)
  hom_alt <- snp_table[snp_table$allele_ratio > 0.75 & snp_table$DP >= 20, ]

  if (nrow(hom_alt) < 50) {
    return(NA_real_)
  }

  # 2. 5 Mb balanced-region selection
  chrom_clean <- gsub("^chr", "", as.character(hom_alt$chromosome))
  hom_alt_copy <- data.table::copy(hom_alt)
  hom_alt_copy[, window_id := paste0(chrom_clean, ":", floor(hom_alt$start / 5e6))]

  # Count het-like SNPs per window (from full snp_table, not just hom-alt)
  all_snps_copy <- data.table::copy(snp_table)
  all_chrom <- gsub("^chr", "", as.character(all_snps_copy$chromosome))
  all_snps_copy[, window_id := paste0(all_chrom, ":", floor(all_snps_copy$start / 5e6))]
  all_snps_copy[, maf := 0.5 + abs(allele_ratio - 0.5)]

  balanced_counts <- all_snps_copy[maf < 0.6, .N, by = window_id]
  data.table::setorder(balanced_counts, -N)

  n_windows <- length(unique(hom_alt_copy$window_id))
  n_select <- max(1, floor(n_windows * 0.10))
  selected_windows <- balanced_counts[seq_len(min(n_select, nrow(balanced_counts)))]$window_id

  sel <- hom_alt_copy[window_id %in% selected_windows]
  if (nrow(sel) < 5) sel <- hom_alt_copy  # fallback to all hom-alt
  if (nrow(sel) < 3) return(NA_real_)

  # 3. Build VAF histogram (10 log-spaced bins, 75-100%)
  bin_breaks_asc <- c(75, 83, 88, 92, 95, 97, 98, 99, 99.5, 99.8, 100)
  n_bins <- length(bin_breaks_asc) - 1
  bin_names <- c(
    "vaf_99.8_100", "vaf_99.5_99.8", "vaf_99_99.5",
    "vaf_98_99", "vaf_97_98", "vaf_95_97",
    "vaf_92_95", "vaf_88_92", "vaf_83_88", "vaf_75_83"
  )

  vaf_pct <- sel$allele_ratio * 100
  bin_idx <- findInterval(vaf_pct, bin_breaks_asc, rightmost.closed = TRUE)
  bin_idx_flipped <- length(bin_breaks_asc) - bin_idx
  bin_idx_flipped[bin_idx_flipped < 1] <- 1
  bin_idx_flipped[bin_idx_flipped > n_bins] <- n_bins

  bin_counts <- tabulate(bin_idx_flipped, nbins = n_bins)
  total <- sum(bin_counts)
  bin_fracs <- if (total > 0) bin_counts / total else rep(0, n_bins)

  # 4. Build feature vector
  features <- data.table::data.table(
    mean_DP = mean(sel$DP),
    sd_DP = stats::sd(sel$DP)
  )
  for (i in seq_along(bin_names)) {
    data.table::set(features, j = bin_names[i], value = bin_fracs[i])
  }

  # 5. Predict
  pred <- tryCatch({
    stats::predict(model_data$model, newdata = features)
  }, error = function(e) {
    warning("Contamination prediction failed: ", e$message)
    return(NA_real_)
  })

  if (is.na(pred)) return(NA_real_)

  # Clamp to [0, 1]
  result <- max(0, min(1, as.numeric(pred)))
  return(result)
}
