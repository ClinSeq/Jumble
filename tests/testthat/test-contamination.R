test_that("estimate_contamination returns NA for NULL input", {
  expect_identical(estimate_contamination(NULL), NA_real_)
})

test_that("estimate_contamination returns NA for empty data.table", {
  empty_dt <- data.table::data.table(
    chromosome = character(0), start = integer(0),
    allele_ratio = numeric(0), DP = integer(0),
    AD = integer(0), RD = integer(0)
  )
  expect_identical(estimate_contamination(empty_dt), NA_real_)
})

test_that("estimate_contamination returns NA when required columns are missing", {
  dt <- data.table::data.table(
    chromosome = "1", start = 1000L, allele_ratio = 0.99
    # Missing DP, AD, RD

  )
  expect_identical(estimate_contamination(dt), NA_real_)
})

test_that("estimate_contamination returns NA when no hom-alt SNPs", {
  # All SNPs are het (allele_ratio ~0.5)
  n <- 200
  dt <- data.table::data.table(
    chromosome = rep("1", n),
    start = seq(1e6, by = 1000, length.out = n),
    allele_ratio = runif(n, 0.4, 0.6),
    DP = rep(100L, n),
    AD = rep(50L, n),
    RD = rep(50L, n)
  )
  expect_identical(estimate_contamination(dt), NA_real_)
})

test_that("estimate_contamination returns NA when fewer than 50 hom-alt SNPs", {
  n <- 30  # Only 30 hom-alt SNPs
  dt <- data.table::data.table(
    chromosome = rep("1", n),
    start = seq(1e6, by = 1000, length.out = n),
    allele_ratio = rep(0.99, n),
    DP = rep(100L, n),
    AD = rep(99L, n),
    RD = rep(1L, n)
  )
  expect_identical(estimate_contamination(dt), NA_real_)
})

test_that("estimate_contamination returns numeric scalar for clean synthetic sample", {
  skip_if_not_installed("randomForest")

  # Simulate a clean sample: 500 hom-alt SNPs with allele_ratio ~1.0
  # Plus some het SNPs for balanced-region selection
  set.seed(42)
  n_hom <- 500
  n_het <- 200

  hom_snps <- data.table::data.table(
    chromosome = rep("1", n_hom),
    start = seq(1e6, by = 10000, length.out = n_hom),
    allele_ratio = pmin(1.0, rnorm(n_hom, mean = 0.998, sd = 0.003)),
    DP = as.integer(rpois(n_hom, 200)),
    AD = integer(n_hom),
    RD = integer(n_hom)
  )
  hom_snps[, AD := as.integer(round(allele_ratio * DP))]
  hom_snps[, RD := DP - AD]

  het_snps <- data.table::data.table(
    chromosome = rep("1", n_het),
    start = seq(1e6 + 5000, by = 10000, length.out = n_het),
    allele_ratio = rnorm(n_het, mean = 0.5, sd = 0.05),
    DP = as.integer(rpois(n_het, 200)),
    AD = integer(n_het),
    RD = integer(n_het)
  )
  het_snps[, AD := as.integer(round(allele_ratio * DP))]
  het_snps[, RD := DP - AD]

  dt <- rbind(hom_snps, het_snps)

  result <- estimate_contamination(dt)

  expect_true(is.numeric(result))
  expect_length(result, 1)
  expect_true(result >= 0 && result <= 1)
  # A clean sample should predict low contamination
  expect_true(result < 0.05)
})

test_that("estimate_contamination result is clamped to [0, 1]", {
  skip_if_not_installed("randomForest")

  # Any valid prediction should be in [0, 1]
  set.seed(123)
  n_hom <- 500
  n_het <- 200

  hom_snps <- data.table::data.table(
    chromosome = rep("2", n_hom),
    start = seq(1e6, by = 10000, length.out = n_hom),
    allele_ratio = pmin(1.0, rnorm(n_hom, mean = 0.95, sd = 0.02)),
    DP = as.integer(rpois(n_hom, 150)),
    AD = integer(n_hom),
    RD = integer(n_hom)
  )
  hom_snps[, AD := as.integer(round(allele_ratio * DP))]
  hom_snps[, RD := DP - AD]

  het_snps <- data.table::data.table(
    chromosome = rep("2", n_het),
    start = seq(1e6 + 5000, by = 10000, length.out = n_het),
    allele_ratio = rnorm(n_het, mean = 0.5, sd = 0.05),
    DP = as.integer(rpois(n_het, 150)),
    AD = integer(n_het),
    RD = integer(n_het)
  )
  het_snps[, AD := as.integer(round(allele_ratio * DP))]
  het_snps[, RD := DP - AD]

  dt <- rbind(hom_snps, het_snps)

  result <- estimate_contamination(dt)

  expect_true(is.numeric(result))
  expect_true(result >= 0)
  expect_true(result <= 1)
})

test_that("load_contamination_model caches the model", {
  # Reset cache
  .jumble_contam_env$model <- NULL

  model1 <- load_contamination_model()
  model2 <- load_contamination_model()

  # Both should return the same object (cached)
  if (!is.null(model1)) {
    expect_identical(model1, model2)
  }
})
