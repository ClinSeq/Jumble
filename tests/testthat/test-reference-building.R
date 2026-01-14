library(testthat)
library(Jumble)

test_that("Reference building works for gene panel", {
  # Use package test data
  testdata_dir <- system.file("testdata", package = "Jumble")
  if (testdata_dir == "") testdata_dir <- "inst/testdata" # For development

  ref_samples <- list.files(
    file.path(testdata_dir, "gene_panel/reference"),
    pattern = "\\.counts\\.RDS$",
    full.names = TRUE
  )

  skip_if(length(ref_samples) == 0, "No reference samples found")
  skip_if_offline()

  out_dir <- tempdir()
  ref_file <- file.path(out_dir, "test_gene_panel_ref.RDS")

  reference <- build_reference(
    count_files = ref_samples,
    annotation_source = "biomart",
    # genome = "hg19", # Let it default
    output_file = ref_file,
    cores = 1
  )

  expect_true(file.exists(ref_file))
  expect_true("target_template" %in% names(reference))
  expect_true("ranges" %in% names(reference))
  expect_true("chromlength" %in% names(reference))
  expect_equal(length(reference$samples), length(ref_samples))
  expect_gt(nrow(reference$target_template), 10000) # Should have many bins
})

test_that("Reference building works for WGS", {
  # Use package test data
  testdata_dir <- system.file("testdata", package = "Jumble")
  if (testdata_dir == "") testdata_dir <- "inst/testdata" # For development

  wgs_samples <- list.files(
    file.path(testdata_dir, "wgs/reference"),
    pattern = "\\.counts\\.RDS$",
    full.names = TRUE
  )

  skip_if(length(wgs_samples) == 0, "No WGS reference samples found")
  skip_if_offline()

  out_dir <- tempdir()
  ref_file <- file.path(out_dir, "test_wgs_ref.RDS")

  reference <- build_reference(
    count_files = wgs_samples,
    annotation_source = "biomart",
    # genome used to be "hg19", now defaulting to auto/hg19 since counts don't have genome yet
    output_file = ref_file,
    cores = 1
  )

  expect_true(file.exists(ref_file))
  expect_true("target_template" %in% names(reference))
  expect_equal(length(reference$samples), length(wgs_samples))
  expect_gt(nrow(reference$target_template), 100000) # WGS should have many more bins
})
