test_that("Full pipeline runs on test data", {
  # Paths
  bam_file <- system.file("extdata", "test_sample.bam", package = "Jumble")
  # If not installed yet, use local path
  if (bam_file == "") bam_file <- "../../inst/extdata/test_sample.bam"

  expect_true(file.exists(bam_file))

  # Output dir
  out_dir <- tempdir()

  # 1. Generate Counts
  # We use WGS mode for simplicity in this test as we don't have a BED file for the test sample
  # (unless we create a dummy one, but WGS is fine for testing pipeline mechanics)

  counts <- generate_counts(bam_file, target_bed = NULL, wgs_bin_size = 10000)

  expect_type(counts, "list")
  expect_true("count" %in% names(counts))
  expect_true("ranges" %in% names(counts))

  # Save counts for reference building
  count_file <- file.path(out_dir, "test_sample.counts.RDS")
  saveRDS(counts, count_file)

  # 2. Build Reference
  # We use the same sample as reference (self-reference) for testing mechanics
  # We need to mock annotation or use biomaRt?
  # BiomaRt might fail without internet or if Ensembl is down.
  # We should use the fallback if possible, or mock it.

  # For testing, let's try to use biomaRt but skip if offline?
  # Or better, use a small mocked annotation if possible.
  # But build_reference expects specific columns.

  # Let's try running with biomaRt. If it fails, we skip.
  skip_if_offline()

  # We need to handle the case where we don't want to download huge data during test.
  # Maybe we can mock generate_gene_annotation?
  # For now, let's assume it works or we skip.

  # Note: The test sample is likely hg19 or hg38. We should detect it.
  genome <- detect_genome(bam_file)
  if (is.null(genome)) genome <- "hg19" # Default

  # Create a dummy reference object manually to avoid biomaRt dependency in basic test?
  # No, we want to test build_reference.

  # Let's try to run it.
  ref_file <- file.path(out_dir, "reference.RDS")

  # We might need to mock generate_gene_annotation to return empty tables to avoid network.
  # But let's try the real thing first.

  # To avoid long wait, maybe we can mock it.
  # But I can't easily mock internal functions in testthat without 'mockery'.

  # Let's just run it.
  reference <- build_reference(
    count_files = c(count_file),
    annotation_source = "biomart",
    genome = genome,
    output_file = ref_file
  )

  expect_true(file.exists(ref_file))
  expect_true("target_template" %in% names(reference))

  # 3. Run Jumble
  # Run on the same sample
  res <- run_jumble(bam_file, ref_file, output_dir = out_dir)

  expect_type(res, "list")
  expect_true("targets" %in% names(res))
  expect_true("segments" %in% names(res))

  # Check output files
  # Check output files
  sample_name <- sub("\\.bam$", "", basename(bam_file), ignore.case = TRUE)
  expect_true(file.exists(file.path(out_dir, paste0(sample_name, ".jumble.csv"))))
  expect_true(file.exists(file.path(out_dir, paste0(sample_name, ".png"))))
})
