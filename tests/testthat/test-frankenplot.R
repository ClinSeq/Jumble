# Tests for standalone frankenplot() function
# These tests verify the helper functions and input validation
# without requiring actual VCF files or rendering.

test_that("frankenplot_cancer_genes returns expected gene list", {
  genes <- JumbleCNV:::frankenplot_cancer_genes()
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

test_that("fp_select_snp_sample handles single-sample VCF-like input", {
  mock_vcf <- matrix(nrow = 0, ncol = 1)
  colnames(mock_vcf) <- "SAMPLE1"

  expect_identical(JumbleCNV:::fp_select_snp_sample(mock_vcf, role = "tumor"), 1L)
  expect_identical(JumbleCNV:::fp_select_snp_sample(mock_vcf, role = "normal"), 1L)
})

test_that("fp_annotate_effect classifies variants correctly", {
  dt <- data.table::data.table(
    CANONICAL = c("YES", "YES", "YES", "YES", ""),
    SYMBOL = c("TP53", "BRCA1", "KRAS", "MYC", ""),
    IMPACT = c("HIGH", "MODERATE", "MODERATE", "LOW", "HIGH"),
    CLIN_SIG = c("", "", "pathogenic", "", ""),
    is_hotspot = c(FALSE, FALSE, FALSE, FALSE, FALSE)
  )

  result <- JumbleCNV:::fp_annotate_effect(dt)

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

  result <- JumbleCNV:::fp_annotate_effect(dt)
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

  result <- JumbleCNV:::fp_annotate_effect(dt)
  expect_equal(nrow(result), 0)
})

test_that("fp_read_jumble_csv reads and processes correctly", {
  test_csv <- tempfile(fileext = ".jumble.csv")
  on.exit(unlink(test_csv), add = TRUE)

  data.table::fwrite(
    data.table::data.table(
      chromosome = rep("1", 30),
      start = seq(1000, by = 1000, length.out = 30),
      end = seq(1100, by = 1000, length.out = 30),
      gene = c(rep("GENE1", 10), rep("Antitarget", 20)),
      log2 = rep(0, 30),
      count = rep(100L, 30)
    ),
    test_csv
  )

  bins <- JumbleCNV:::fp_read_jumble_csv(test_csv)

  expect_s3_class(bins, "data.table")
  expect_true("bin" %in% names(bins))
  expect_true("log2" %in% names(bins))
  expect_true("smooth_log2" %in% names(bins))
  expect_true("target type" %in% names(bins))
  expect_true("depth" %in% names(bins))
  expect_true(all(bins$bin == seq_len(nrow(bins))))
  expect_true("Background" %in% bins$gene)
})

test_that("fp_read_jumble_csv stops on missing file", {
  expect_error(
    JumbleCNV:::fp_read_jumble_csv("/nonexistent/file.csv"),
    "not found"
  )
})

test_that("fp_read_cns maps segments to bins", {
  test_csv <- tempfile(fileext = ".jumble.csv")
  test_cns <- tempfile(fileext = ".cns")
  on.exit(unlink(c(test_csv, test_cns)), add = TRUE)

  data.table::fwrite(
    data.table::data.table(
      chromosome = rep("1", 30),
      start = seq(1000, by = 1000, length.out = 30),
      end = seq(1100, by = 1000, length.out = 30),
      gene = rep("GENE1", 30),
      log2 = rep(0, 30),
      count = rep(100L, 30)
    ),
    test_csv
  )
  data.table::fwrite(
    data.table::data.table(
      chromosome = "1",
      start = 1000,
      end = 10100,
      log2 = 0
    ),
    test_cns
  )

  bins <- JumbleCNV:::fp_read_jumble_csv(test_csv)
  result <- JumbleCNV:::fp_read_cns(test_cns, bins)

  expect_type(result, "list")
  expect_true("bins" %in% names(result))
  expect_true("segments" %in% names(result))
  expect_true("segment" %in% names(result$bins))
  expect_gt(sum(!is.na(result$bins$segment)), 0)
})

test_that("fp_map_snps_to_bins handles NULL input", {
  bins <- data.table::data.table(
    chromosome = c("1", "1"),
    start = c(100, 200),
    end = c(199, 299),
    bin = 1:2
  )

  result <- JumbleCNV:::fp_map_snps_to_bins(NULL, bins)
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

  result <- JumbleCNV:::fp_map_variants_to_bins(NULL, bins)
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

  result <- JumbleCNV:::fp_map_variants_to_bins(variants, bins)
  expect_equal(nrow(result), 0)
})

test_that("fp_parse_dpyd handles missing files gracefully", {
  result <- JumbleCNV:::fp_parse_dpyd(NULL, NULL)
  expect_null(result$dpyd_result)
  expect_null(result$dpyd_table)

  result2 <- JumbleCNV:::fp_parse_dpyd("/nonexistent.json", "/nonexistent.csv")
  expect_null(result2$dpyd_result)
  expect_null(result2$dpyd_table)
})

test_that("frankenplot validates required arguments", {
  expect_error(
    frankenplot(output_file = tempfile(fileext = ".html"), tumor_cns = "test.cns"),
    "tumor_jumble_csv"
  )

  expect_error(
    frankenplot(tumor_jumble_csv = "/nonexistent.csv",
                tumor_cns = "test.cns",
                output_file = tempfile(fileext = ".html")),
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
                output_file = tempfile(fileext = ".html"),
                somatic_vcf = "/nonexistent.vcf"),
    "somatic_vcf.*not found"
  )

  unlink(c(tmp_csv, tmp_cns))
})

test_that("frankenplot validates output_png length", {
  tmp_csv <- tempfile(fileext = ".jumble.csv")
  tmp_cns <- tempfile(fileext = ".cns")
  on.exit(unlink(c(tmp_csv, tmp_cns)), add = TRUE)

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
                output_file = tempfile(fileext = ".html"),
                output_png = c("a.png", "b.png")),
    "output_png must be a single file path"
  )
})

test_that("fp_render_report passes static PNG path to the template", {
  skip_if_not_installed("rmarkdown")

  tmp_dir <- tempfile("frankenplot-render-")
  dir.create(tmp_dir)
  html_file <- file.path(tmp_dir, "sample.frankenplot.html")
  png_file <- file.path(tmp_dir, "sample.frankenplot.png")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  bins <- data.table::data.table(
    chromosome = rep(c("1", "2"), each = 20),
    start = rep(seq(1, by = 1000, length.out = 20), 2),
    end = rep(seq(500, by = 1000, length.out = 20), 2),
    gene = rep(c("TP53", "Background"), 20),
    log2 = rep(0, 40),
    depth = rep(100, 40),
    bin = seq_len(40),
    smooth_log2 = rep(0, 40),
    `target type` = rep(c("target", "background"), 20),
    type = rep(c("target", "background"), 20),
    allele_ratio = as.numeric(NA)
  )
  segments <- data.table::data.table(
    chromosome = c("1", "2"),
    start = c(1L, 21L),
    end = c(20L, 40L),
    start_pos = c(1, 1),
    end_pos = c(20000, 20000),
    log2 = c(0, 0),
    gstart = c(1, 20001),
    gstop = c(20000, 40000)
  )

  JumbleCNV:::fp_render_report(
    bins_t = data.table::copy(bins),
    segments_t = data.table::copy(segments),
    bins_n = NULL,
    segments_n = NULL,
    alf = NULL,
    alf_n = NULL,
    salf = NULL,
    galf_t = NULL,
    galf_n = NULL,
    hrdtable = NULL,
    qc_metrics = NULL,
    dpyd_result = NULL,
    dpyd_table = NULL,
    genome = "hg19",
    output_file = html_file,
    tumor_jumble_csv = "sample.jumble.csv",
    output_png = png_file
  )

  expect_true(file.exists(html_file))
  expect_true(file.exists(png_file))
  expect_gt(file.info(png_file)$size, 0)
})

test_that("fp_parse_germline_vcf returns correct structure for NULL input", {
  result <- JumbleCNV:::fp_parse_germline_vcf(NULL)
  expect_type(result, "list")
  expect_null(result$galf_n)
  expect_null(result$galf_t)
})

test_that("fp_parse_somatic_vcf returns NULL for missing file", {
  result <- JumbleCNV:::fp_parse_somatic_vcf(NULL)
  expect_null(result)

  result2 <- JumbleCNV:::fp_parse_somatic_vcf("/nonexistent.vcf")
  expect_null(result2)
})

test_that("fp_parse_snp_vcf returns NULL for missing file", {
  result <- JumbleCNV:::fp_parse_snp_vcf(NULL, 1)
  expect_null(result)

  result2 <- JumbleCNV:::fp_parse_snp_vcf("/nonexistent.vcf", 1)
  expect_null(result2)
})
