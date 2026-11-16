library(testthat)
library(Jumble)

test_that("Regression test with anonymized test data", {
  # Use package test data
  testdata_dir <- system.file("testdata", package = "Jumble")
  if (testdata_dir == "") testdata_dir <- "inst/testdata" # For development

  ref_file <- file.path(testdata_dir, "gene_panel/reference.RDS")
  samples_dir <- file.path(testdata_dir, "gene_panel/samples")

  # Check if files exist
  if (!file.exists(ref_file)) {
    skip("Test data not found, skipping regression test")
  }

  # Get list of count files
  count_files <- list.files(samples_dir, pattern = "\\.counts\\.RDS$", full.names = TRUE)

  if (length(count_files) == 0) {
    skip("No count files found in test data")
  }

  # Output dir for this test
  out_dir <- tempfile()
  dir.create(out_dir)

  for (count_file in count_files) {
    sample_name <- basename(count_file)
    sample_id <- sub("\\.counts\\.RDS$", "", sample_name)

    # Check for VCF
    vcf_pattern <- paste0(sample_id, "\\.vcf")
    vcf_files <- list.files(samples_dir, pattern = vcf_pattern, full.names = TRUE)
    vcf_file <- if (length(vcf_files) > 0) vcf_files[1] else NULL

    message(sprintf("Testing sample: %s", sample_id))

    # Run Jumble
    res <- run_jumble(
      bam_file = count_file,
      reference_file = ref_file,
      output_dir = out_dir,
      snp_vcf = vcf_file
    )

    # Check output exists
    out_csv <- file.path(out_dir, paste0(sample_id, ".jumble.csv"))
    expect_true(file.exists(out_csv))

    # Load results
    actual <- data.table::fread(out_csv)

    # --- Check Output Format ---
    expect_true("log2" %in% names(actual),
      info = paste(sample_id, "missing log2 column")
    )
    expect_true("bin" %in% names(actual))
    expect_true("chromosome" %in% names(actual))
    expect_true(nrow(actual) > 0)

    # Check for segmentation output
    has_segment <- "segment" %in% names(actual)
    has_smooth <- "smooth_log2" %in% names(actual)
    expect_true(has_segment || has_smooth,
      info = paste(sample_id, "missing segment/smooth_log2 column")
    )

    # Check segments were generated
    # Check segments were generated
    expect_true(nrow(actual[!is.na(get(if (has_segment) "segment" else "smooth_log2"))]) > 0,
      info = paste(sample_id, "no segments generated")
    )

    # If VCF was provided, check SNP file
    if (!is.null(vcf_file)) {
      snp_file <- file.path(out_dir, paste0(sample_id, ".jumble_snps.csv"))
      expect_true(file.exists(snp_file),
        info = paste(sample_id, "SNP file not created despite VCF input")
      )
    }

    # Check GIS output
    if (!is.null(vcf_file)) {
      gis_file <- file.path(out_dir, paste0(sample_id, ".jumble_gis.csv"))
      expect_true(file.exists(gis_file),
        info = paste(sample_id, "GIS file not created")
      )
  
      if (file.exists(gis_file)) {
        gis_dt <- data.table::fread(gis_file)
        expect_true(all(c("fraction", "predicted_gis") %in% names(gis_dt)),
          info = "GIS table missing required columns"
        )
        expect_true(nrow(gis_dt) > 0)
      }
    }
  }
})
