library(testthat)
library(Jumble)

test_that("SNP processing works with VCF file", {
  # Use package test data
  testdata_dir <- system.file("testdata", package = "Jumble")
  if (testdata_dir == "") testdata_dir <- "inst/testdata" # For development

  ref_file <- file.path(testdata_dir, "gene_panel/reference.RDS")
  sample_file <- file.path(testdata_dir, "gene_panel/samples/test_sample_1.counts.RDS")
  vcf_file <- file.path(testdata_dir, "gene_panel/samples/test_sample_1.vcf.gz")

  skip_if(!file.exists(ref_file), "Reference file not found")
  skip_if(!file.exists(sample_file), "Sample file not found")
  skip_if(!file.exists(vcf_file), "VCF file not found")

  out_dir <- tempdir()

  res <- run_jumble(
    bam_file = sample_file,
    reference_file = ref_file,
    output_dir = out_dir,
    snp_vcf = vcf_file
  )

  # Check SNP CSV file was created
  sample_name <- "test_sample_1"
  snp_file <- file.path(out_dir, paste0(sample_name, ".jumble_snps.csv"))

  expect_true(file.exists(snp_file))

  # Load and check SNP data
  snp_data <- data.table::fread(snp_file)
  expect_true(is.data.frame(snp_data) || is.data.table(snp_data))
  expect_true("allele_ratio" %in% names(snp_data))
  expect_true("chromosome" %in% names(snp_data))
  expect_gt(nrow(snp_data), 0)
})

test_that("Pipeline works without SNP VCF", {
  # Use package test data
  testdata_dir <- system.file("testdata", package = "Jumble")
  if (testdata_dir == "") testdata_dir <- "inst/testdata" # For development

  ref_file <- file.path(testdata_dir, "gene_panel/reference.RDS")
  sample_file <- file.path(testdata_dir, "gene_panel/samples/test_sample_2.counts.RDS")

  skip_if(!file.exists(ref_file), "Reference file not found")
  skip_if(!file.exists(sample_file), "Sample file not found")

  out_dir <- tempdir()

  res <- run_jumble(
    bam_file = sample_file,
    reference_file = ref_file,
    output_dir = out_dir,
    snp_vcf = NULL
  )

  # Check SNP CSV file was NOT created
  sample_name <- "test_sample_2"
  snp_file <- file.path(out_dir, paste0(sample_name, ".jumble_snps.csv"))

  expect_false(file.exists(snp_file))

  # But other files should exist
  expect_true(file.exists(file.path(out_dir, paste0(sample_name, ".jumble.csv"))))
  expect_true(file.exists(file.path(out_dir, paste0(sample_name, ".png"))))
})
