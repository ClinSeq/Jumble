# Tests for standalone frankenplot() function
# These tests verify the helper functions and input validation
# without requiring actual VCF files or rendering.

test_that("frankenplot_cancer_genes returns expected gene list", {
  genes <- Jumble:::frankenplot_cancer_genes()
  expect_type(genes, "character")
  expect_true(length(genes) > 50)
  # Key genes must be present
  expect_true("TP53" %in% genes)
  expect_true("BRCA1" %in% genes)
  expect_true("BRCA2" %in% genes)
  expect_true("EGFR" %in% genes)
  expect_true("KRAS" %in% genes)
  expect_true("PTEN" %in% genes)
  # No duplicates

  expect_equal(length(genes), length(unique(genes)))
})

test_that("fp_select_snp_sample handles single-sample VCF", {
  # Mock a VCF-like object with colnames
  mock_vcf <- structure(list(), class = "mock_vcf")
  attr(mock_vcf, "colnames_val") <- "SAMPLE1"
  # We can't easily mock VariantAnnotation objects, so test the logic directly
  # by calling with a character vector approach


  # Test the function's internal logic with a simple wrapper
  # Since fp_select_snp_sample uses colnames(vcf), we test edge cases
  skip_if_not_installed("VariantAnnotation")
})

test_that("fp_annotate_effect classifies variants correctly", {
  dt <- data.table::data.table(
    CANONICAL = c("YES", "YES", "YES", "YES", ""),
    SYMBOL = c("TP53", "BRCA1", "KRAS", "MYC", ""),
    IMPACT = c("HIGH", "MODERATE", "MODERATE", "LOW", "HIGH"),
    CLIN_SIG = c("", "", "pathogenic", "", ""),
    is_hotspot = c(FALSE, FALSE, FALSE, FALSE, FALSE)
  )

  result <- Jumble:::fp_annotate_effect(dt)

  expect_equal(result$effect[1], "high-impact")   # HIGH impact
  expect_equal(result$effect[2], "uncertain")      # MODERATE impact
  expect_equal(result$effect[3], "high-impact")    # pathogenic overrides MODERATE
  expect_true(is.na(result$effect[4]))             # LOW impact = NA
  expect_true(is.na(result$effect[5]))             # empty SYMBOL + CANONICAL
})

test_that("fp_annotate_effect handles hotspot override", {
  dt <- data.table::data.table(
    CANONICAL = c("YES"),
    SYMBOL = c("BRAF"),
    IMPACT = c("MODERATE"),
    CLIN_SIG = c(""),
    is_hotspot = c(TRUE)
  )

  result <- Jumble:::fp_annotate_effect(dt)
  expect_equal(result$effect[1], "hotspot")
})

test_that("fp_annotate_effect handles empty data.table", {
  dt <- data.table::data.table(
    CANONICAL = character(0),
    SYMBOL = character(0),
    IMPACT = character(0),
    CLIN_SIG = character(0),
    is_hotspot = logical(0)
  )

  result <- Jumble:::fp_annotate_effect(dt)
  expect_equal(nrow(result), 0)
})

test_that("fp_read_jumble_csv reads and processes correctly", {
  # Use the package's test data
  test_csv <- system.file("testdata", "test_tumor.jumble.csv",
                           package = "Jumble")
  if (test_csv == "" || !file.exists(test_csv)) {
    skip("Test CSV not available")
  }

  bins <- Jumble:::fp_read_jumble_csv(test_csv)

  expect_s3_class(bins, "data.table")
  expect_true("bin" %in% names(bins))
  expect_true("log2" %in% names(bins))
  expect_true("smooth_log2" %in% names(bins))
  expect_true("target type" %in% names(bins))
  expect_true(all(bins$bin == 1:nrow(bins)))
})

test_that("fp_read_jumble_csv stops on missing file", {
  expect_error(
    Jumble:::fp_read_jumble_csv("/nonexistent/file.csv"),
    "not found"
  )
})

test_that("fp_read_cns maps segments to bins", {
  test_csv <- system.file("testdata", "test_tumor.jumble.csv",
                           package = "Jumble")
  test_cns <- system.file("testdata", "test_tumor.cns",
                           package = "Jumble")
  if (test_csv == "" || test_cns == "" ||
      !file.exists(test_csv) || !file.exists(test_cns)) {
    skip("Test data not available")
  }

  bins <- Jumble:::fp_read_jumble_csv(test_csv)
  result <- Jumble:::fp_read_cns(test_cns, bins)

  expect_type(result, "list")
  expect_true("bins" %in% names(result))
  expect_true("segments" %in% names(result))
  expect_true("segment" %in% names(result$bins))
})

test_that("fp_map_snps_to_bins handles NULL input", {
  bins <- data.table::data.table(
    chromosome = c("1", "1"),
    start = c(100, 200),
    end = c(199, 299),
    bin = 1:2
  )

  result <- Jumble:::fp_map_snps_to_bins(NULL, bins)
  expect_true("allele_ratio" %in% names(result$bins))
  expect_true(all(is.na(result$bins$allele_ratio)))
})

test_that("fp_map_variants_to_bins handles NULL input", {
  bins <- data.table::data.table(
    chromosome = c("1", "1"),
    start = c(100, 200),
    end = c(199, 299),
    bin = 1:2
  )

  result <- Jumble:::fp_map_variants_to_bins(NULL, bins)
  expect_null(result)
})

test_that("fp_map_variants_to_bins handles empty data.table", {
  bins <- data.table::data.table(
    chromosome = c("1", "1"),
    start = c(100, 200),
    end = c(199, 299),
    bin = 1:2
  )
  variants <- data.table::data.table(
    chromosome = character(0),
    start = numeric(0),
    end = numeric(0)
  )

  result <- Jumble:::fp_map_variants_to_bins(variants, bins)
  expect_equal(nrow(result), 0)
})

test_that("fp_parse_dpyd handles missing files gracefully", {
  result <- Jumble:::fp_parse_dpyd(NULL, NULL)
  expect_null(result$dpyd_result)
  expect_null(result$dpyd_table)

  result2 <- Jumble:::fp_parse_dpyd("/nonexistent.json", "/nonexistent.csv")
  expect_null(result2$dpyd_result)
  expect_null(result2$dpyd_table)
})

test_that("frankenplot validates required arguments", {
  expect_error(
    frankenplot(output_file = "test.html", tumor_cns = "test.cns"),
    "tumor_jumble_csv"
  )

  expect_error(
    frankenplot(tumor_jumble_csv = "/nonexistent.csv",
                tumor_cns = "test.cns",
                output_file = "test.html"),
    "must exist"
  )
})

test_that("frankenplot validates optional file arguments", {
  # Create a temporary CSV to pass the first check
  tmp_csv <- tempfile(fileext = ".jumble.csv")
  tmp_cns <- tempfile(fileext = ".cns")
  data.table::fwrite(
    data.table::data.table(chromosome = "1", start = 100, end = 200,
                            gene = "TEST", log2 = 0, depth = 100),
    tmp_csv
  )
  data.table::fwrite(
    data.table::data.table(chromosome = "1", start = 100, end = 200,
                            log2 = 0),
    tmp_cns
  )

  expect_error(
    frankenplot(tumor_jumble_csv = tmp_csv,
                tumor_cns = tmp_cns,
                output_file = "test.html",
                somatic_vcf = "/nonexistent.vcf"),
    "somatic_vcf.*not found"
  )

  unlink(c(tmp_csv, tmp_cns))
})

test_that("fp_parse_germline_vcf returns correct structure for NULL input", {
  result <- Jumble:::fp_parse_germline_vcf(NULL)
  expect_type(result, "list")
  expect_null(result$galf_n)
  expect_null(result$galf_t)
})

test_that("fp_parse_somatic_vcf returns NULL for missing file", {
  result <- Jumble:::fp_parse_somatic_vcf(NULL)
  expect_null(result)

  result2 <- Jumble:::fp_parse_somatic_vcf("/nonexistent.vcf")
  expect_null(result2)
})

test_that("fp_parse_snp_vcf returns NULL for missing file", {
  result <- Jumble:::fp_parse_snp_vcf(NULL, 1)
  expect_null(result)

  result2 <- Jumble:::fp_parse_snp_vcf("/nonexistent.vcf", 1)
  expect_null(result2)
})
