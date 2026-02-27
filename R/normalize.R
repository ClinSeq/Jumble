#' Reconstruct Targets from Reference Counts
#'
#' Rebuilds target data.table from reference allcounts and template
#'
#' @param reference Reference object with allcounts and target_template
#' @return data.table with reconstructed targets from all reference samples
#' @keywords internal
reconstruct_targets_from_reference <- function(reference) {
  allcounts <- reference$allcounts
  target_template <- reference$target_template
  
  targetlist <- list()
  for (i in seq_along(allcounts)) {
    counts <- allcounts[[i]]
    t <- data.table::copy(target_template)
    t[, sample := paste0("ref_", i)]
    t[, count := counts$count]
    t[, count_short := counts$count_short]
    
    targetlist[[i]] <- t
  }
  
  data.table::rbindlist(targetlist)
}

#' Apply Value Floor
#'
#' Helper to replace values below 1 with 1 (for log safety)
#'
#' @param x Numeric vector
#' @return Vector with floor of 1 and NA replaced
#' @keywords internal
apply_value_floor <- function(x) {
  x[x < 1] <- 1
  x[is.na(x)] <- 1
  x
}

#' Filter Bins by Coverage Thresholds
#'
#' Remove bins with extremely low or high coverage
#'
#' @param targets data.table with count column
#' @return Filtered targets
#' @keywords internal
filter_bins_by_coverage <- function(targets) {
  # Low coverage threshold
  threshold_low <- median(targets[is_target == TRUE]$count) * 0.01
  keep_bins_low <- targets[, median(count), by = bin][V1 > threshold_low]$bin
  targets <- targets[bin %in% keep_bins_low]
  
  # High coverage threshold
  threshold_high <- median(targets[is_target == TRUE]$count) / 0.05
  keep_bins_high <- targets[, median(count), by = bin][V1 < threshold_high]$bin
  targets <- targets[bin %in% keep_bins_high]
  
  targets
}

#' Calculate Log Ratio
#'
#' Convert counts to log2 ratio with safety checks
#'
#' @param targets data.table with count column
#' @param count_col Column name to convert
#' @param output_col Name for output column
#' @return Modified targets with output_col added
#' @keywords internal
calculate_logr <- function(targets, count_col, output_col) {
  targets[[output_col]] <- log2(apply_value_floor(targets[[count_col]]))
  targets
}

weighted_median <- function(x, w) {
  valid <- is.finite(x) & is.finite(w) & w > 0
  if (sum(valid) == 0) return(NA_real_)
  x <- x[valid]
  w <- w[valid]
  o <- order(x)
  x <- x[o]
  w <- w[o]
  p <- cumsum(w) / sum(w)
  x[which(p >= 0.5)[1]]
}

#' Median Correct by Backbone
#'
#' Subtract weighted median of backbone bins from each sample/target group
#'
#' @param targets data.table with backbone_weight and is_target columns
#' @param lr_col Column name with log ratio values
#' @return Modified targets
#' @keywords internal
median_correct_backbone <- function(targets, lr_col) {
  # Ensure we have a proper copy (avoids data.table shallow copy warning)
  targets <- data.table::copy(targets)
  
  targets[, (lr_col) := get(lr_col) - weighted_median(get(lr_col), backbone_weight),
    by = c("sample", "is_target")
  ]
  targets
}

#' Detect Sample Gender
#'
#' Infer sample gender from X chromosome median (for X correction)
#' and Y chromosome median (for Y correction), independently.
#' Follows the original jumble logic where X-based and Y-based gender
#' assignments are used for their respective chromosome corrections.
#'
#' @param targets data.table
#' @param lr_col Column name with log ratio values
#' @return data.table with male_from_x, male_from_y columns
#' @keywords internal
detect_sample_gender <- function(targets, lr_col) {
  # Ensure we have a proper copy (avoids data.table shallow copy warning)
  targets <- data.table::copy(targets)
  
  # X chromosome median (untiled bins only, as in original target bins)
  targets[, xmedian := median(get(lr_col)[chromosome == "X" & is_tiled == FALSE],
    na.rm = TRUE
  ), by = "sample"]
  targets[, male_from_x := 2^xmedian < 0.75]
  
  # Y chromosome median
  targets[, ymedian := median(get(lr_col)[chromosome == "Y"], na.rm = TRUE),
    by = "sample"
  ]
  targets[, male_from_y := 2^ymedian > 0.25]
  
  targets
}

#' Apply Sex Chromosome Corrections
#'
#' Adjust X chromosome using X-based gender detection and
#' Y chromosome using Y-based gender detection, independently.
#' This matches the original jumble logic.
#'
#' @param targets data.table with male_from_x and male_from_y columns
#' @param lr_col Column name with log ratio values
#' @return Modified targets with sex-chromosome-corrected values
#' @keywords internal
correct_sex_chromosomes <- function(targets, lr_col) {
  # Ensure we have a proper copy (avoids data.table shallow copy warning)
  targets <- data.table::copy(targets)
  
  # Hard-coded PAR (Pseudoautosomal Region) boundaries
  targets[, nonPA := chromosome %in% c("X") & end > 2.70e6 & start < 154.93e6]
  
  # X chromosome correction for males (using X-based gender)
  if (any(targets$male_from_x == TRUE, na.rm = TRUE)) {
    targets[chromosome == "X" & male_from_x == TRUE & nonPA, 
            (lr_col) := get(lr_col) + 1]
  }
  
  # Y chromosome correction for males (using Y-based gender)
  if (any(targets$male_from_y == TRUE, na.rm = TRUE)) {
    targets[male_from_y == TRUE & chromosome == "Y" & end < 28.79e6, 
            (lr_col) := get(lr_col) + 1]
  }
  
  # Y values to NA for females (using Y-based gender)
  targets[chromosome == "Y" & male_from_y == FALSE, (lr_col) := NA]
  
  targets
}

#' Clean Temporary Gender Detection Columns
#'
#' @param targets data.table
#' @return Modified targets with temp columns removed
#' @keywords internal
cleanup_gender_columns <- function(targets) {
  # Ensure we have a proper copy (avoids data.table shallow copy warning)
  targets <- data.table::copy(targets)
  cols_to_remove <- intersect(
    c("xmedian", "ymedian", "male_from_x", "male_from_y", "male", "nonPA"),
    names(targets)
  )
  if (length(cols_to_remove) > 0) targets[, (cols_to_remove) := NULL]
  targets
}

#' Impute Missing Values
#'
#' Replace NA in log ratio with random noise near 0
#'
#' @param targets data.table
#' @param lr_col Column name with log ratios
#' @return Modified targets with imputed values
#' @keywords internal
impute_missing_logr <- function(targets, lr_col) {
  # Ensure we have a proper copy (avoids data.table shallow copy warning)
  targets <- data.table::copy(targets)
  targets[is.na(get(lr_col)), 
          (lr_col) := rnorm(.N, mean = 0, sd = 0.01)]
  targets
}

#' Apply Bin Median Correction
#'
#' Subtract bin-wise median (reference correction)
#'
#' @param targets data.table
#' @param lr_col Column name with log ratios
#' @param ref_col Name for output reference median column
#' @return Modified targets with reference-corrected values
#' @keywords internal
apply_bin_median_correction <- function(targets, lr_col, ref_col) {
  # Ensure we have a proper copy (avoids data.table shallow copy warning)
  targets <- data.table::copy(targets)
  targets[, (ref_col) := median(get(lr_col), na.rm = TRUE), by = bin]
  targets[, (lr_col) := get(lr_col) - get(ref_col)]
  targets
}

#' Cast Matrix for PCA
#'
#' Convert long-format data to wide matrix for PCA
#'
#' @param targets data.table with bin, sample, and value columns
#' @param id_cols Identifier columns (typically "bin")
#' @param formula_cols Formula for dcast (e.g., "bin ~ sample")
#' @param value_col Column to spread
#' @return data.table in wide format with 'bin' and sample columns
#' @keywords internal
cast_to_matrix <- function(targets, formula_cols, value_col) {
  data.table::dcast(targets[, c("bin", "sample", value_col), with = FALSE],
                    formula = formula_cols,
                    value.var = value_col)
}

#' Perform PCA on Subset
#'
#' Helper to perform PCA on specific bins (targets or background)
#'
#' @param matrix data.table with bin column and sample columns (wide format)
#' @param bin_list Vector of bins to include in PCA
#' @return data.table with PCA scores (PC1, PC2, ...) and bin column, or NULL
#' @keywords internal
perform_pca_on_bins <- function(matrix, bin_list) {
  if (length(bin_list) == 0) return(NULL)
  
  submatrix <- matrix[bin %in% bin_list, -1, with = FALSE]
  if (ncol(submatrix) < 2) return(NULL)  # Need at least 2 samples
  
  pca_result <- stats::prcomp(submatrix, center = FALSE, scale. = FALSE)
  dt <- as.data.table(pca_result$x)
  dt$bin <- matrix[bin %in% bin_list]$bin
  
  dt
}

#' Prepare Bin Metadata for Correction
#'
#' Extract reference median values for later normalization
#'
#' @param targets data.table with refmedian columns
#' @return data.table with unique bin-level metadata
#' @keywords internal
prepare_bin_metadata <- function(targets) {
  unique(targets[, .(bin, refmedian, refmedian_short)])
}

#' Filter Valid Bins
#'
#' Identify bins that survived QC filtering
#'
#' @param matrix data.table with bin column
#' @return Vector of valid bins
#' @keywords internal
extract_valid_bins <- function(matrix) {
  matrix$bin
}

#' Extract Bin Subset for PCA
#'
#' Get bins of specific type (target or background)
#'
#' @param target_template data.table with bin column
#' @param is_target_bool Logical: TRUE for targets, FALSE for background
#' @param valid_bins Bins that passed QC
#' @return Vector of bins of specified type
#' @keywords internal
get_bins_by_type <- function(target_template, is_target_bool, valid_bins = NULL) {
  bins <- target_template[is_target == is_target_bool]$bin
  if (!is.null(valid_bins)) {
    bins <- bins[bins %in% valid_bins]
  }
  bins
}

#' Apply PCA Correction to Values using Optimisation
#'
#' L1 + Total Variation penalty optimisation using Nelder-Mead
#'
#' @param data data.table with 'lr' and PC columns
#' @param ratio TV penalty ratio (default 1.0 based on benchmarking)
#' @return Corrected log ratio vector
#' @keywords internal
correct_by_optim <- function(data, ratio = 1.0) {
  pcs <- sum(grepl("^PC", colnames(data)))
  if (pcs == 0) return(data$lr)

  # Extract PC matrix
  pc_cols <- paste0("PC", seq_len(pcs))
  pc_mat <- as.matrix(data[, pc_cols, with = FALSE])
  lr_vec <- data$lr

  # Filter out NAs for optimization stability
  keep <- is.finite(lr_vec) & complete.cases(pc_mat)
  pc_mat_opt <- pc_mat[keep, , drop = FALSE]
  lr_vec_opt <- lr_vec[keep]

  # Cost function: L1 norm + Total Variation
  make_objective <- function(pm, lv, r) {
    function(mc) {
      n <- length(mc)
      new <- lv - as.numeric(pm[, seq_len(n), drop = FALSE] %*% mc)
      sum(abs(new)) + r * sum(abs(diff(new)))
    }
  }
  obj_fn <- make_objective(pc_mat_opt, lr_vec_opt, ratio)

  # Progressive coefficient building: start small, expand
  n <- min(3, pcs)
  mc <- rep(0, n)
  mc <- stats::optim(par = mc, fn = obj_fn, method = "Nelder-Mead")$par

  if (pcs > n) {
    while (n < pcs) {
      old_n <- n
      n <- min(n + 10, pcs)
      mc <- c(mc, rep(0, n - old_n))
      mc <- stats::optim(par = mc, fn = obj_fn, method = "Nelder-Mead")$par
    }
  }

  correction <- as.numeric(pc_mat[, seq_len(length(mc)), drop = FALSE] %*% mc)
  data$lr - correction
}

#' Apply PCA Correction to Values
#'
#' Robust linear model correction using PCA components
#'
#' @param data data.table with 'lr' and PC columns
#' @param train_indices Logical or numeric indices for training subset
#' @return Corrected log ratio vector
#' @keywords internal
correct_by_pca <- function(data, train_indices = NULL) {
  if (is.null(train_indices)) {
    train_indices <- rep(TRUE, nrow(data))
  }
  
  # Subsample if too large
  if (is.logical(train_indices) && length(which(train_indices)) > 20000) {
    train_indices <- sample(which(train_indices), 20000)
  } else if (!is.logical(train_indices) && length(train_indices) > 20000) {
    train_indices <- sample(train_indices, 20000)
  }
  
  # Count PCs
  pcs <- sum(grepl("^PC", colnames(data)))
  if (pcs == 0) return(data$lr)
  
  # Build formula
  formula_str <- paste("lr ~ ", paste0("PC", seq_len(min(12, pcs)), collapse = "+"))
  
  # Robust linear regression
  rlm_mod <- tryCatch(
    {
      MASS::rlm(as.formula(formula_str), data = data, subset = train_indices)
    },
    error = function(e) NULL
  )
  
  if (!is.null(rlm_mod)) {
    data[, lr := lr - stats::predict(rlm_mod, data)]
  }
  
  data$lr
}

#' Apply GC Content Correction
#'
#' Loess smoothing correction for GC content bias
#'
#' @param data data.table with 'lr', 'gc' and 'backbone_weight' columns
#' @param train_indices Logical or numeric indices for training
#' @param span Loess span parameter
#' @return Corrected log ratio vector
#' @keywords internal
correct_by_gc <- function(data, train_indices = NULL, span = 0.75) {
  if (is.null(train_indices)) {
    train_indices <- rep(TRUE, nrow(data))
  }
  
  # Identify valid rows
  valid_rows <- is.finite(data$lr) & is.finite(data$gc)
  train_indices_clean <- train_indices & valid_rows
  
  n_points <- sum(train_indices_clean)
  if (n_points <= 5) return(data$lr)  # Not enough points for loess
  
  w <- if ("backbone_weight" %in% names(data)) data$backbone_weight[train_indices_clean] else NULL
  
  # Adjust span for small datasets
  if (n_points < 50) span <- 1.0
  
  # Fit loess on valid training points
  tryCatch(
    {
      if (!is.null(w)) {
        loess_mod <- stats::loess(lr ~ gc, 
          data = data, subset = train_indices_clean, weights = w,
          span = span, family = "symmetric",
          control = stats::loess.control(surface = "interpolate")
        )
      } else {
        loess_mod <- stats::loess(lr ~ gc, 
          data = data, subset = train_indices_clean,
          span = span, family = "symmetric",
          control = stats::loess.control(surface = "interpolate")
        )
      }
      
      # Predict and handle NAs
      preds <- stats::predict(loess_mod, data)
      preds[is.na(preds)] <- 0
      data[, lr := lr - preds]
    },
    error = function(e) warning("Loess failed: ", e$message)
  )
  
  data$lr
}

#' Apply Combined PCA and GC Corrections
#'
#' Orchestrates both PCA and GC corrections
#'
#' @param data data.table with lr, gc, and PCA component columns
#' @param train_indices Logical indices for training subset
#' @param correction String indicating the method: "optim" (L1+TV) or "rlm" (Robust LM)
#' @return Corrected log ratio vector
#' @keywords internal
apply_combined_corrections <- function(data, train_indices = NULL, correction = "optim") {
  if (is.null(train_indices)) {
    train_indices <- rep(TRUE, nrow(data))
  }
  
  # 1. PCA correction
  if (correction == "optim") {
    corrected <- correct_by_optim(data, ratio = 1.0)
  } else {
    corrected <- correct_by_pca(data, train_indices)
  }
  data[, lr := corrected]
  
  # 2. GC correction
  corrected <- correct_by_gc(data, train_indices)
  
  corrected
}

#' Apply Correction for Subset
#'
#' Main function that orchestrates PCA + GC correction for targets or background
#'
#' @param targets data.table with all required columns
#' @param pca_data PCA results data.table
#' @param lr_col Column name with log ratios
#' @param output_col Column name for corrected output
#' @param is_target_bool Logical: TRUE for targets, FALSE for background
#' @param correction String indicating the method: "optim" (L1+TV) or "rlm" (Robust LM)
#' @return Modified targets with output_col updated
#' @keywords internal
apply_normalization_corrections <- function(targets, pca_data, 
                                           lr_col, output_col, 
                                           is_target_bool,
                                           correction = "optim") {
  if (is.null(pca_data)) {
    targets[[output_col]] <- targets[[lr_col]]
    return(targets)
  }
  
  # Filter to subset
  ix <- targets$is_target == is_target_bool
  if (!any(ix)) return(targets)
  
  # Merge with PCA
  data <- merge(targets[ix], pca_data, by = "bin", all.x = TRUE)
  data[, lr := get(lr_col)]
  
  # Apply corrections
  corrected <- apply_combined_corrections(data, data$is_backbone, correction)
  
  # Assign back
  data[, corrected_lr := corrected]
  targets[data, (output_col) := i.corrected_lr, on = "bin"]
  
  targets
}

#' Clamp Log Ratio Values
#'
#' Set min/max bounds on log ratios
#'
#' @param targets data.table
#' @param lr_col Column to clamp
#' @param min_val Minimum value
#' @param max_val Maximum value
#' @return Modified targets
#' @keywords internal
clamp_logr_values <- function(targets, lr_col, min_val, max_val) {
  # Ensure we have a proper copy (avoids data.table shallow copy warning)
  targets <- data.table::copy(targets)
  targets[get(lr_col) < min_val, (lr_col) := min_val]
  targets[get(lr_col) > max_val, (lr_col) := max_val]
  targets
}

#' Sort Targets Genomically
#'
#' Order by chromosome and genomic position
#'
#' @param targets data.table with chromosome and start columns
#' @return Sorted targets
#' @keywords internal
sort_genomically <- function(targets) {
  # Ensure we have a proper copy (avoids data.table shallow copy warning)
  targets <- data.table::copy(targets)
  
  chrom_levels <- c(as.character(1:22), "X", "Y")
  
  # Standardize chromosome names
  targets[, sort_chr := as.character(chromosome)]
  targets[, sort_chr := stringr::str_remove(sort_chr, "^chr")]
  targets[, sort_fac := factor(sort_chr, levels = chrom_levels)]
  
  # Sort by chromosome factor, then start position
  data.table::setorder(targets, sort_fac, start, na.last = TRUE)
  
  # Clean temp columns
  targets[, c("sort_chr", "sort_fac") := NULL]
  
  targets
}

#' Clean Up Temporary Columns
#'
#' Remove intermediate columns used during normalization
#'
#' @param targets data.table
#' @param cols_to_remove Character vector of column names
#' @return Modified targets
#' @keywords internal
cleanup_temp_columns <- function(targets, cols_to_remove) {
  # Ensure we have a proper copy (avoids data.table shallow copy warning)
  targets <- data.table::copy(targets)
  
  cols_exist <- cols_to_remove[cols_to_remove %in% names(targets)]
  if (length(cols_exist) > 0) {
    targets[, (cols_exist) := NULL]
  }
  targets
}

#' Compute Reference PCA (Refactored)
#'
#' Performs PCA on reference dataset to identify latent features for normalization.
#'
#' @param reference The reference object.
#' @return A list containing PCA results for targets and background.
#' @importFrom stats prcomp sd
#' @importFrom data.table dcast as.data.table setorder
#' @export
compute_reference_pca <- function(reference) {
  # 1. Reconstruct and define backbone
  targets <- reconstruct_targets_from_reference(reference)
  targets <- define_backbone(targets)
  
  # 2. Filter bins by coverage
  targets <- filter_bins_by_coverage(targets)
  
  # 3. Calculate log ratios
  targets <- calculate_logr(targets, "count", "rawLR")
  targets <- calculate_logr(targets, "count_short", "rawLR_short")
  
  # 4. Median correct to backbone
  targets <- median_correct_backbone(targets, "rawLR")
  targets <- median_correct_backbone(targets, "rawLR_short")
  
  # 5. Apply sex chromosome corrections
  targets <- detect_sample_gender(targets, "rawLR")
  targets <- correct_sex_chromosomes(targets, "rawLR")
  targets <- detect_sample_gender(targets, "rawLR_short")
  targets <- correct_sex_chromosomes(targets, "rawLR_short")
  targets <- cleanup_gender_columns(targets)
  
  # 6. Impute and median correct
  targets <- impute_missing_logr(targets, "rawLR")
  targets <- impute_missing_logr(targets, "rawLR_short")
  targets <- apply_bin_median_correction(targets, "rawLR", "refmedian")
  targets <- apply_bin_median_correction(targets, "rawLR_short", "refmedian_short")
  
  # 7. Cast to matrices
  mat <- cast_to_matrix(targets, bin ~ sample, "rawLR")
  mat_short <- cast_to_matrix(targets, bin ~ sample, "rawLR_short")
  
  # 8. Get bin lists
  target_template <- reference$target_template
  targetbins <- get_bins_by_type(target_template, TRUE, mat$bin)
  backgroundbins <- get_bins_by_type(target_template, FALSE, mat$bin)
  
  # 9. Perform PCA for latent factors
  set.seed(25)
  tpca <- perform_pca_on_bins(mat, targetbins)
  tpca_short <- perform_pca_on_bins(mat_short, targetbins)
  bgpca <- perform_pca_on_bins(mat, backgroundbins)
  bgpca_short <- perform_pca_on_bins(mat_short, backgroundbins)
  
  # 11. Prepare metadata
  bins_meta <- prepare_bin_metadata(targets)
  valid_bins <- extract_valid_bins(mat)
  
  list(
    tpca = tpca,
    tpca_short = tpca_short,
    bgpca = bgpca,
    bgpca_short = bgpca_short,
    bins_meta = bins_meta,
    valid_bins = valid_bins
  )
}

#' Normalize Sample Data (Refactored)
#'
#' Normalizes the query sample using reference PCA and GC correction.
#'
#' @param targets Query sample data (data.table).
#' @param reference_pca PCA results from compute_reference_pca.
#' @param correction String indicating the method: "optim" (L1+TV) or "rlm" (Robust LM)
#' @return Normalized targets.
#' @importFrom MASS rlm
#' @importFrom stats predict loess loess.control median
#' @export
normalize_sample <- function(targets, reference_pca, correction = "optim") {
  # 1. Filter and prepare targets
  targets <- targets[bin %in% reference_pca$valid_bins]
  targets <- merge(targets, reference_pca$bins_meta, by = "bin", all.x = TRUE)
  
  # 2. Calculate log ratios
  targets <- calculate_logr(targets, "count", "rawLR")
  targets <- calculate_logr(targets, "count_short", "rawLR_short")
  
  # 3. Median correct to backbone
  targets <- median_correct_backbone(targets, "rawLR")
  targets <- median_correct_backbone(targets, "rawLR_short")
  
  # 4. Subtract reference median
  targets[, rawLR := rawLR - refmedian]
  targets[, rawLR_short := rawLR_short - refmedian_short]
  
  # 5. Impute missing values
  targets <- impute_missing_logr(targets, "rawLR")
  targets <- impute_missing_logr(targets, "rawLR_short")
  
  # 6. Apply normalization corrections (targets)
  targets <- apply_normalization_corrections(targets, reference_pca$tpca,
                                             "rawLR", "log2", TRUE, correction)
  targets <- apply_normalization_corrections(targets, reference_pca$tpca_short,
                                             "rawLR_short", "log2_short", TRUE, correction)
  
  # 7. Apply normalization corrections (background)
  targets <- apply_normalization_corrections(targets, reference_pca$bgpca,
                                             "rawLR", "log2", FALSE, correction)
  targets <- apply_normalization_corrections(targets, reference_pca$bgpca_short,
                                             "rawLR_short", "log2_short", FALSE, correction)
  
  # 8. Clamp values
  targets <- clamp_logr_values(targets, "log2", -5, 7)
  targets <- clamp_logr_values(targets, "log2_short", -4, 7)
  
  # 9. Sort genomically
  targets <- sort_genomically(targets)
  
  return(targets)
}
