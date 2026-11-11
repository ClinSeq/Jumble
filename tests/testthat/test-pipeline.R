test_that("Full pipeline runs on test data", {
  # Paths
  bam_file <- system.file("extdata", "test_sample.bam", package = "JumbleCNV")
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

  # Save counts for reference-building mechanics. The generated WGS count object
  # is intentionally not used for the run_jumble() smoke test because a
  # single-sample self-reference is numerically unstable for PCA normalization.
  count_file <- file.path(out_dir, "test_sample.counts.RDS")
  saveRDS(counts, count_file)

  # 2. Build a small WGS reference with local empty annotation so this formal
  # package test remains deterministic and offline-safe.
  genome <- detect_genome(bam_file)
  if (is.null(genome)) genome <- "hg19"

  wgs_ref_file <- file.path(out_dir, "wgs_reference.RDS")
  annotation_file <- file.path(out_dir, "empty_annotation.RDS")
  saveRDS(
    list(
      allgenes = data.table::data.table(),
      allexons = data.table::data.table(),
      cancergenes_clinseq = data.table::data.table()
    ),
    annotation_file
  )

  wgs_reference <- build_reference(
    count_files = c(count_file),
    annotation_source = annotation_file,
    genome = genome,
    output_file = wgs_ref_file
  )

  expect_true(file.exists(wgs_ref_file))
  expect_true("target_template" %in% names(wgs_reference))
  expect_equal(wgs_reference$target_bed_file, "wgs")

  # 3. Run Jumble using stable bundled gene-panel fixtures. This keeps the
  # end-to-end run_jumble() smoke test independent of the synthetic WGS
  # self-reference above.
  testdata_dir <- system.file("testdata", package = "JumbleCNV")
  if (testdata_dir == "") testdata_dir <- "inst/testdata"
  sample_file <- file.path(testdata_dir, "gene_panel/samples/test_sample_1.counts.RDS")
  ref_file <- file.path(testdata_dir, "gene_panel/reference.RDS")

  res <- run_jumble(
    sample_file,
    ref_file,
    output_dir = out_dir,
    annotation_mode = "none"
  )

  expect_type(res, "list")
  expect_true("targets" %in% names(res))
  expect_true("segments" %in% names(res))

  # Check output files
  sample_name <- sub("\\.counts\\.RDS$", "", basename(sample_file), ignore.case = TRUE)
  expect_true(file.exists(file.path(out_dir, paste0(sample_name, ".jumble.csv"))))
  expect_true(file.exists(file.path(out_dir, paste0(sample_name, ".png"))))
})

test_that("generate_counts works with targeted mode (BED file)", {
  bam_file <- system.file("extdata", "test_sample.bam", package = "JumbleCNV")
  if (bam_file == "") bam_file <- "../../inst/extdata/test_sample.bam"
  
  expect_true(file.exists(bam_file))
  
  # Create a mock BED file
  mock_bed_path <- tempfile(fileext = ".bed")
  bed_content <- data.frame(
    chromosome = c("1", "2"),
    start = c(10000, 50000),
    end = c(20000, 60000)
  )
  write.table(bed_content, mock_bed_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  # Ensure the mock bed file exists
  expect_true(file.exists(mock_bed_path))
  
  # Test generate_counts in targeted mode
  counts <- generate_counts(bam_file, target_bed = mock_bed_path)
  
  expect_type(counts, "list")
  expect_true("count" %in% names(counts))
  expect_true("ranges" %in% names(counts))
  expect_true("bed" %in% names(counts))
  expect_false(is.null(counts$bed))
})
