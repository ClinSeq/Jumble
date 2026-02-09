#' Compute Reference PCA
#'
#' Performs PCA on the reference dataset to identify latent features for
#' normalization.
#'
#' @param reference The reference object.
#' @return A list containing PCA results for targets and background.
#' @importFrom stats prcomp sd
#' @importFrom data.table dcast as.data.table
#' @export
compute_reference_pca <- function(reference) {
  # Extract targets from reference
  # We need to reconstruct the matrix from reference$target_template and
  # reference$allcounts
  # Or maybe reference already has the matrix?
  # In build_reference, we saved 'allcounts'. We need to build the matrix here.

  # Reconstruct targets from allcounts
  # This logic was in jumble-run.R lines 229-247

  # 1. Extract Counts --------------------------------------------------------
  allcounts <- reference$allcounts
  target_template <- reference$target_template

  targetlist <- list()
  for (i in seq_along(allcounts)) {
    counts <- allcounts[[i]]
    # We need to map counts to the template bins
    # Assuming bins match exactly as checked in build_reference

    t <- data.table::copy(target_template)
    t[, sample := paste0("ref_", i)]

    # Add counts
    t[, count := counts$count]
    if (!is.null(counts$count_medium)) {
      # Logic for medium fraglength?
      # In run.R it was conditional. Let's stick to standard count for PCA for
      # now?
      # Or maybe we should support it.
      # For simplicity, let's use 'count' (all fragments) and 'count_short'.
    }
    t[, count_short := counts$count_short]

    targetlist[[i]] <- t
  }
  targets <- data.table::rbindlist(targetlist)

  # Define backbone
  targets <- define_backbone(targets)

  # Filter bins (remove worst) - logic from run.R lines 256-268
  # We should probably do this once and save it in reference?
  # But run.R does it dynamically.

  # 2. Filter Bins -----------------------------------------------------------
  # Low coverage
  threshold_low <- median(targets[is_target == TRUE]$count) * 0.01
  keep_bins_low <- targets[, median(count), by = bin][V1 > threshold_low]$bin
  targets <- targets[bin %in% keep_bins_low]

  # High coverage
  threshold_high <- median(targets[is_target == TRUE]$count) / 0.05
  keep_bins_high <- targets[, median(count), by = bin][V1 < threshold_high]$bin
  targets <- targets[bin %in% keep_bins_high]

  # (Mappability filter removed - not essential for normalization)

  # 3. Initial Corrections ---------------------------------------------------
  # LogR
  min1 <- function(x) {
    x[x < 1] <- 1
    x[is.na(x)] <- 1
    x
  }
  targets[, rawLR := log2(min1(count))]
  targets[, rawLR_short := log2(min1(count_short))]

  # Median correct to backbone
  targets[, rawLR := rawLR - median(rawLR[is_backbone]),
    by = c("sample", "is_target")
  ]
  targets[, rawLR_short := rawLR_short - median(rawLR_short[is_backbone]),
    by = c("sample", "is_target")
  ]

  # X-Y chromosome correct (matching original script lines 299-328)
  # Doctor the X values in samples where their median implies male
  targets[, xmedian := median(rawLR[chromosome == "X" & is_tiled == FALSE],
    na.rm = TRUE
  ), by = c("sample", "is_target")]
  targets[, xmedian_short := median(rawLR_short[chromosome == "X" &
    is_tiled == FALSE], na.rm = TRUE), by = c("sample", "is_target")]
  targets[, male := 2^xmedian < 0.75] # assign gender

  # hard coded PA
  targets[, nonPA := chromosome %in% c("X") & end > 2.70e6 & start < 154.93e6]

  # adjust to the median
  if (length(unique(targets[male == TRUE]$sample)) > 0) { # if at least 1 male
    targets[chromosome == "X" & male == TRUE & nonPA, rawLR := rawLR + 1]
    targets[
      chromosome == "X" & male == TRUE & nonPA,
      rawLR_short := rawLR_short + 1
    ]
  }

  # Doctor Y values where their median implies male
  targets[, ymedian := median(rawLR[chromosome == "Y"], na.rm = TRUE),
    by = "sample"
  ] # median by sample
  targets[, male := 2^ymedian > 0.25] # assign gender based on Y

  if (length(unique(targets[male == TRUE]$sample)) > 0) { # if at least 1 male
    # hard coded PA
    targets[male == TRUE & chromosome == "Y" & end < 28.79e6, rawLR := rawLR + 1]
    targets[
      male == TRUE & chromosome == "Y" & end < 28.79e6,
      rawLR_short := rawLR_short + 1
    ]
  }

  # Y values set to NA where X median implied female
  targets[chromosome == "Y" & male == FALSE, rawLR := NA]
  targets[chromosome == "Y" & male == FALSE, rawLR_short := NA]

  # Clean up temp columns
  targets[, c("xmedian", "xmedian_short", "ymedian", "male", "nonPA") := NULL]

  # 4. Impute & Med Correct --------------------------------------------------
  # Impute missing (BEFORE bin median correct!)
  # If any missing, replace with random value near 0 (matching original script
  # lines 338-340)
  targets[is.na(rawLR), rawLR := rnorm(n = .N, mean = 0, sd = 0.01)]
  targets[is.na(rawLR_short), rawLR_short := rnorm(n = .N, mean = 0, sd = 0.01)]

  # Bin median correct (AFTER X/Y correction!)
  targets[, refmedian := median(rawLR, na.rm = TRUE), by = bin]
  targets[, rawLR := rawLR - refmedian]

  targets[, refmedian_short := median(rawLR_short, na.rm = TRUE), by = bin]
  targets[, rawLR_short := rawLR_short - refmedian_short]

  # Matrix
  mat <- dcast(targets[, .(bin, sample, rawLR)], bin ~ sample,
    value.var = "rawLR"
  )
  mat_short <- dcast(targets[, .(bin, sample, rawLR_short)], bin ~ sample,
    value.var = "rawLR_short"
  )

  # 5. Outlier Removal (PCA 1) -----------------------------------------------
  targetbins <- target_template[is_target == TRUE]$bin
  # Only keep those that survived filtering
  targetbins <- targetbins[targetbins %in% mat$bin]

  # Helper for PCA
  do_pca <- function(m, bins) {
    subm <- m[bin %in% bins, -1, with = FALSE]
    if (ncol(subm) < 2) {
      return(NULL)
    } # Need at least 2 samples
    res <- prcomp(subm, center = FALSE, scale. = FALSE)
    dt <- as.data.table(res$x)
    dt$bin <- m[bin %in% bins]$bin
    return(dt)
  }

  # Initial PCA for outlier detection
  set.seed(25)
  tpca <- do_pca(mat, targetbins)

  remove_bins <- c()

  if (!is.null(tpca)) {
    tpca[, keep := TRUE]
    pcs <- colnames(tpca)[grep("^PC", colnames(tpca))]
    for (pc in pcs[1:min(100, length(pcs))]) {
      fact <- ifelse(pc %in% c("PC1", "PC2", "PC3"), 3, 4)
      sd_val <- sd(tpca[[pc]])
      tpca[get(pc) < -sd_val * fact, keep := FALSE]
      tpca[get(pc) > sd_val * fact, keep := FALSE]
    }
    remove_bins <- mat[bin %in% targetbins]$bin[tpca$keep == FALSE]
  }

  # Background outliers
  backgroundbins <- target_template[is_target == FALSE]$bin
  backgroundbins <- backgroundbins[backgroundbins %in% mat$bin]

  if (length(backgroundbins) > 0) {
    bgpca <- do_pca(mat, backgroundbins)
    if (!is.null(bgpca)) {
      bgpca[, keep := TRUE]
      pcs <- colnames(bgpca)[grep("^PC", colnames(bgpca))]
      for (pc in pcs[1:min(100, length(pcs))]) {
        fact <- ifelse(pc %in% c("PC1", "PC2", "PC3"), 3, 4)
        sd_val <- sd(bgpca[[pc]])
        bgpca[get(pc) < -sd_val * fact, keep := FALSE]
        bgpca[get(pc) > sd_val * fact, keep := FALSE]
      }
      remove_bins <- c(
        remove_bins,
        mat[bin %in% backgroundbins]$bin[bgpca$keep == FALSE]
      )
    }
  }

  # Remove outliers
  if (length(remove_bins) > 0) {
    mat <- mat[!bin %in% remove_bins]
    mat_short <- mat_short[!bin %in% remove_bins]
  }

  valid_bins <- mat$bin

  # 6. Latent Factors (PCA 2) ------------------------------------------------
  set.seed(25)
  tpca <- do_pca(mat, targetbins)
  tpca_short <- do_pca(mat_short, targetbins)

  bgpca <- NULL
  bgpca_short <- NULL
  if (length(backgroundbins) > 0) {
    bgpca <- do_pca(mat, backgroundbins)
    bgpca_short <- do_pca(mat_short, backgroundbins)
  }

  # We need to return the 'bins_for_mediancorrect' as well
  bins_meta <- unique(targets[, .(bin, refmedian, refmedian_short)])

  list(
    tpca = tpca,
    tpca_short = tpca_short,
    bgpca = bgpca,
    bgpca_short = bgpca_short,
    bins_meta = bins_meta,
    valid_bins = valid_bins
  )
}

#' Normalize Sample Data
#'
#' Normalizes the query sample using reference PCA and GC correction.
#'
#' @param targets Query sample data (data.table).
#' @param reference_pca PCA results from compute_reference_pca.
#' @return Normalized targets.
#' @importFrom MASS rlm
#' @importFrom stats predict as.formula loess loess.control median
#' @export
normalize_sample <- function(targets, reference_pca) {
  # 1. Prepare Targets -------------------------------------------------------
  # Filter bins
  targets <- targets[bin %in% reference_pca$valid_bins]

  # Merge refmedian
  targets <- merge(targets, reference_pca$bins_meta, by = "bin", all.x = TRUE)

  # LogR
  min1 <- function(x) {
    x[x < 1] <- 1
    x[is.na(x)] <- 1
    x
  }
  targets[, rawLR := log2(min1(count))]
  targets[, rawLR_short := log2(min1(count_short))]

  # Median correct to backbone
  targets[, rawLR := rawLR - median(rawLR[is_backbone]), by = "is_target"]
  targets[, rawLR_short := rawLR_short - median(rawLR_short[is_backbone]),
    by = "is_target"
  ]

  # 2. Apply Reference Corrections -------------------------------------------
  # Correct by reference median
  targets[, rawLR := rawLR - refmedian]
  targets[, rawLR_short := rawLR_short - refmedian_short]

  # Impute missing (random noise near 0)
  targets[is.na(rawLR), rawLR := rnorm(.N, mean = 0, sd = 0.01)]
  targets[is.na(rawLR_short), rawLR_short := rnorm(.N, mean = 0, sd = 0.01)]

  # PCA Correction Function
  jcorrect <- function(temp, train_ix = NULL) {
    if (is.null(train_ix)) train_ix <- rep(TRUE, nrow(temp))
    if (is.logical(train_ix) && length(which(train_ix)) > 20000) {
      train_ix_idx <- which(train_ix)
      train_ix <- sample(train_ix_idx, 20000)
    } else if (length(train_ix) > 20000 && !is.logical(train_ix)) {
      train_ix <- sample(train_ix, 20000)
    }

    pcs <- sum(grepl("^PC", colnames(temp)))
    if (pcs == 0) {
      return(temp$lr)
    }

    formula_str <- paste("lr ~ ", paste0("PC", seq_len(min(12, pcs)),
      collapse = "+"
    ))

    # Robust Linear Model
    rlm_mod <- tryCatch(
      {
        rlm(as.formula(formula_str), data = temp, subset = train_ix)
      },
      error = function(e) NULL
    )

    if (!is.null(rlm_mod)) {
      temp[, lr := lr - predict(rlm_mod, temp)]
    }

    # GC Correct (Loess)
      # Clean potential NAs/Infs in lr before loess
      valid_rows <- is.finite(temp$lr) & is.finite(temp$gc)
      if (sum(valid_rows) > 10) {
          # Update train_ix to intersect with valid rows
          train_ix_clean <- train_ix & valid_rows
          
          # Adjust span for small datasets
          n_points <- sum(train_ix_clean)
          if (n_points > 5) { # Only run if we have enough points
              span_val <- 0.75
              if (n_points < 50) span_val <- 1.0 
              
              tryCatch({
                  loess_mod <- loess(lr ~ gc,
                    data = temp, subset = train_ix_clean, span = span_val,
                    family = "symmetric", control = loess.control(surface = "interpolate")
                  )
                  # Predict safely
                  preds <- predict(loess_mod, temp)
                  # Replace NA preds with 0 (no correction)
                  preds[is.na(preds)] <- 0
                  temp[, lr := lr - preds]
              }, error = function(e) warning("Loess failed: ", e$message))
          }
      }


    return(temp$lr)
  }

  # 3. PCA & GC Correction ---------------------------------------------------
  # Apply Correction (Standard)
  ix <- targets$is_target
  if (!is.null(reference_pca$tpca)) {
    # Merge PCA
    temp <- merge(targets[ix], reference_pca$tpca, by = "bin", all.x = TRUE)
    temp[, lr := rawLR]
    corrected <- jcorrect(temp, temp$is_backbone)

    # Assign back safely
    # We use a join update to ensure alignment
    temp[, corrected_lr := corrected]
    targets[temp, log2 := i.corrected_lr, on = "bin"]
  } else {
    targets[ix, log2 := rawLR]
  }

  # Apply Correction (Short)
  if (!is.null(reference_pca$tpca_short)) {
    temp <- merge(targets[ix], reference_pca$tpca_short,
      by = "bin",
      all.x = TRUE
    )
    temp[, lr := rawLR_short]
    corrected <- jcorrect(temp, temp$is_backbone)
    temp[, corrected_lr := corrected]
    targets[temp, log2_short := i.corrected_lr, on = "bin"]
  } else {
    targets[ix, log2_short := rawLR_short]
  }

  # Background correction
  if (any(!targets$is_target)) {
    ix_bg <- !targets$is_target

    if (!is.null(reference_pca$bgpca)) {
      temp <- merge(targets[ix_bg], reference_pca$bgpca,
        by = "bin",
        all.x = TRUE
      )
      temp[, lr := rawLR]
      corrected <- jcorrect(temp, temp$is_backbone)
      temp[, corrected_lr := corrected]
      targets[temp, log2 := i.corrected_lr, on = "bin"]
    } else {
      targets[ix_bg, log2 := rawLR]
    }

    if (!is.null(reference_pca$bgpca_short)) {
      temp <- merge(targets[ix_bg], reference_pca$bgpca_short,
        by = "bin",
        all.x = TRUE
      )
      temp[, lr := rawLR_short]
      corrected <- jcorrect(temp, temp$is_backbone)
      temp[, corrected_lr := corrected]
      targets[temp, log2_short := i.corrected_lr, on = "bin"]
    } else {
      targets[ix_bg, log2_short := rawLR_short]
    }
  }

  # 4. Clean Up --------------------------------------------------------------
  # Set min/max
  targets[log2 < -5, log2 := -5]
  targets[log2 > 7, log2 := 7]
  targets[log2_short < -4, log2_short := -4]
  targets[log2_short > 7, log2_short := 7]
  targets[chromosome == "Y", log2 := NA]
  targets[chromosome == "Y", log2_short := NA]

  # 5. Sort Genomically ------------------------------------------------------
  # Ensure targets are sorted by chromosome and start for correct plotting/downstream analysis
  
  # Standardize chrom names for sorting
  # Use clean_chrom_names helper (assumed available in package)
  # We handle standard 1..22, X, Y. Others sort at end.
  
  chrom_levels <- c(as.character(1:22), "X", "Y")
  
  # Use temporary columns for sorting
  targets[, sort_chr := as.character(chromosome)]
  targets[, sort_chr := stringr::str_remove(sort_chr, "^chr")]
  targets[, sort_fac := factor(sort_chr, levels = chrom_levels)]
  
  # Sort: Factor first (NA last), then Start
  data.table::setorder(targets, sort_fac, start, na.last = TRUE)
  
  # Clean temp
  targets[, c("sort_chr", "sort_fac") := NULL]

  return(targets)
}
