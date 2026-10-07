library(testthat)
library(JumbleCNV)

test_that("Reference building works for gene panel", {
  # Use package test data
  testdata_dir <- system.file("testdata", package = "JumbleCNV")
  if (testdata_dir == "") testdata_dir <- "inst/testdata" # For development

  ref_samples <- list.files(
    file.path(testdata_dir, "gene_panel/reference"),
    pattern = "\\.counts\\.RDS$",
    full.names = TRUE
  )

  skip_if(length(ref_samples) == 0, "No reference samples found")

  out_dir <- tempdir()
  ref_file <- file.path(out_dir, "test_gene_panel_ref.RDS")
  annotation_file <- file.path(out_dir, "empty_annotation.RDS")

  saveRDS(
    list(
      allgenes = data.table::data.table(),
      allexons = data.table::data.table(),
      cancergenes_clinseq = data.table::data.table()
    ),
    annotation_file
  )

  reference <- build_reference(
    count_files = ref_samples,
    annotation_source = annotation_file,
    genome = "hg19",
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
  bam_file <- system.file("extdata", "test_sample.bam", package = "JumbleCNV")
  if (bam_file == "") bam_file <- "inst/extdata/test_sample.bam"

  expect_true(file.exists(bam_file))

  out_dir <- tempdir()
  count_file <- file.path(out_dir, "test_sample.wgs.counts.RDS")
  ref_file <- file.path(out_dir, "test_wgs_ref.RDS")
  annotation_file <- file.path(out_dir, "empty_annotation.RDS")

  counts <- generate_counts(
    bam_file = bam_file,
    target_bed = NULL,
    wgs_bin_size = 1000000
  )
  saveRDS(counts, count_file)

  saveRDS(
    list(
      allgenes = data.table::data.table(),
      allexons = data.table::data.table(),
      cancergenes_clinseq = data.table::data.table()
    ),
    annotation_file
  )

  reference <- build_reference(
    count_files = count_file,
    annotation_source = annotation_file,
    genome = "hg19",
    output_file = ref_file,
    cores = 1
  )

  expect_true(file.exists(ref_file))
  expect_true("target_template" %in% names(reference))
  expect_equal(length(reference$samples), 1L)
  expect_gt(nrow(reference$target_template), 1000)
  expect_equal(reference$target_bed_file, "wgs")
})
