# Tests for Frankenplot command-line wrapper
# These tests exercise wrapper availability and lightweight CLI validation paths
# without requiring full HTML rendering.

frankenplot_wrapper_path <- function() {
  script_path <- system.file("scripts", "jumble-frankenplot.R", package = "Jumble")
  if (nzchar(script_path) && file.exists(script_path)) {
    return(script_path)
  }

  candidate_paths <- c(
    file.path("..", "..", "inst", "scripts", "jumble-frankenplot.R"),
    file.path("inst", "scripts", "jumble-frankenplot.R")
  )
  candidate_paths <- normalizePath(candidate_paths, mustWork = FALSE)
  candidate_paths[file.exists(candidate_paths)][1]
}

run_frankenplot_wrapper <- function(args) {
  script_path <- frankenplot_wrapper_path()
  skip_if(is.na(script_path) || !file.exists(script_path),
          "Frankenplot wrapper script is not available")

  rscript <- file.path(R.home("bin"), "Rscript")
  out <- tempfile("frankenplot-wrapper-stdout-", fileext = ".txt")
  err <- tempfile("frankenplot-wrapper-stderr-", fileext = ".txt")
  on.exit(unlink(c(out, err)), add = TRUE)

  status <- system2(
    rscript,
    args = c(script_path, args),
    stdout = out,
    stderr = err
  )

  list(
    status = status,
    stdout = paste(readLines(out, warn = FALSE), collapse = "\n"),
    stderr = paste(readLines(err, warn = FALSE), collapse = "\n")
  )
}

test_that("jumble-frankenplot wrapper script is present and parseable", {
  script_path <- frankenplot_wrapper_path()
  expect_false(is.na(script_path))
  expect_true(file.exists(script_path))
  expect_silent(parse(script_path))
})

test_that("jumble-frankenplot wrapper exposes help text", {
  res <- run_frankenplot_wrapper("--help")

  if (!identical(res$status, 0L) && grepl("there is no package called", res$stderr)) {
    skip("Child Rscript process cannot load local Jumble package")
  }

  expect_identical(res$status, 0L)
  expect_match(res$stdout, "Generate a Frankenplot HTML genome report")
  expect_match(res$stdout, "--jumble-csv")
  expect_match(res$stdout, "--cns")
  expect_match(res$stdout, "--output")
  expect_match(res$stdout, "--output-png")
  expect_match(res$stdout, "--somatic-vcf")
  expect_match(res$stdout, "--compute-gis")
  expect_match(res$stdout, "--gis-output")
})

test_that("jumble-frankenplot wrapper validates missing required arguments", {
  res <- run_frankenplot_wrapper(character())

  if (identical(res$status, 0L)) {
    fail("Wrapper unexpectedly succeeded without required arguments")
  }
  if (grepl("there is no package called", res$stderr)) {
    skip("Child Rscript process cannot load local Jumble package")
  }

  expect_match(res$stderr, "--jumble-csv")
})

test_that("jumble-frankenplot wrapper validates optional file arguments before rendering", {
  tmp_csv <- tempfile(fileext = ".jumble.csv")
  tmp_cns <- tempfile(fileext = ".cns")
  tmp_html <- tempfile(fileext = ".html")
  on.exit(unlink(c(tmp_csv, tmp_cns, tmp_html)), add = TRUE)

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

  missing_vcf <- tempfile(fileext = ".vcf")
  res <- run_frankenplot_wrapper(c(
    "--jumble-csv", tmp_csv,
    "--cns", tmp_cns,
    "--output", tmp_html,
    "--somatic-vcf", missing_vcf
  ))

  if (identical(res$status, 0L)) {
    fail("Wrapper unexpectedly succeeded with a missing optional VCF")
  }
  if (grepl("there is no package called", res$stderr)) {
    skip("Child Rscript process cannot load local Jumble package")
  }

  expect_match(res$stderr, "--somatic-vcf file not found")
})

# ---------------------------------------------------------------------------
# GIS computation argument validation (--compute-gis; fails before any GIS
# work or rendering)
# ---------------------------------------------------------------------------

write_wrapper_inputs <- function() {
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
  list(csv = tmp_csv, cns = tmp_cns, html = tempfile(fileext = ".html"))
}

write_snp_vcf <- function(samples) {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=1,length=249250621>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples),
          collapse = "\t"),
    paste(c("1", "150", "rs1", "A", "G", ".", "PASS", ".", "GT:AD:DP",
            rep("0/1:50,50:100", length(samples))), collapse = "\t")
  ), vcf)
  vcf
}

expect_wrapper_error <- function(res, pattern) {
  if (identical(res$status, 0L)) {
    fail("Wrapper unexpectedly succeeded")
  }
  if (grepl("there is no package called", res$stderr)) {
    skip("Child Rscript process cannot load local Jumble package")
  }
  expect_match(res$stderr, pattern, fixed = TRUE)
}

test_that("GIS computation requires --genome", {
  inp <- write_wrapper_inputs()
  vcf <- write_snp_vcf("TUMOR")
  on.exit(unlink(c(unlist(inp), vcf)), add = TRUE)

  res <- run_frankenplot_wrapper(c(
    "--jumble-csv", inp$csv, "--cns", inp$cns, "--output", inp$html,
    "--tumor-snp-vcf", vcf, "--compute-gis"
  ))
  expect_wrapper_error(res, "--genome is required with --compute-gis")
})

test_that("GIS computation requires --tumor-sample-name for multi-sample VCFs", {
  inp <- write_wrapper_inputs()
  vcf <- write_snp_vcf(c("NORMAL", "TUMOR"))
  on.exit(unlink(c(unlist(inp), vcf)), add = TRUE)

  res <- run_frankenplot_wrapper(c(
    "--jumble-csv", inp$csv, "--cns", inp$cns, "--output", inp$html,
    "--tumor-snp-vcf", vcf, "--genome", "hg19", "--compute-gis"
  ))
  expect_wrapper_error(res, "--tumor-sample-name is required")
})

test_that("GIS computation requires an exact --tumor-sample-name match", {
  inp <- write_wrapper_inputs()
  vcf <- write_snp_vcf(c("P1-N-DNA", "P1-T-DNA"))
  on.exit(unlink(c(unlist(inp), vcf)), add = TRUE)

  # Partial name: would fuzzy-match in match_vcf_sample(), must be rejected here.
  res <- run_frankenplot_wrapper(c(
    "--jumble-csv", inp$csv, "--cns", inp$cns, "--output", inp$html,
    "--tumor-snp-vcf", vcf, "--genome", "hg19", "--compute-gis", "--tumor-sample-name", "P1-T"
  ))
  expect_wrapper_error(res, "must exactly match one sample")
})

test_that("GIS computation rejects a non-matching name for single-sample VCFs", {
  inp <- write_wrapper_inputs()
  vcf <- write_snp_vcf("TUMOR")
  on.exit(unlink(c(unlist(inp), vcf)), add = TRUE)

  res <- run_frankenplot_wrapper(c(
    "--jumble-csv", inp$csv, "--cns", inp$cns, "--output", inp$html,
    "--tumor-snp-vcf", vcf, "--genome", "hg19", "--compute-gis", "--tumor-sample-name", "OTHER"
  ))
  expect_wrapper_error(res, "must exactly match one sample")
})

test_that("--compute-gis requires --tumor-snp-vcf", {
  inp <- write_wrapper_inputs()
  on.exit(unlink(unlist(inp)), add = TRUE)

  res <- run_frankenplot_wrapper(c(
    "--jumble-csv", inp$csv, "--cns", inp$cns, "--output", inp$html,
    "--genome", "hg19", "--compute-gis"
  ))
  expect_wrapper_error(res, "--tumor-snp-vcf is required with --compute-gis")
})

test_that("--compute-gis cannot be combined with --hrdtable", {
  inp <- write_wrapper_inputs()
  vcf <- write_snp_vcf("TUMOR")
  hrd <- tempfile(fileext = ".csv")
  writeLines("Measure,Value", hrd)
  on.exit(unlink(c(unlist(inp), vcf, hrd)), add = TRUE)

  res <- run_frankenplot_wrapper(c(
    "--jumble-csv", inp$csv, "--cns", inp$cns, "--output", inp$html,
    "--tumor-snp-vcf", vcf, "--genome", "hg19", "--hrdtable", hrd, "--compute-gis"
  ))
  expect_wrapper_error(res, "--compute-gis cannot be combined with --hrdtable")
})

test_that("--gis-output requires --compute-gis", {
  inp <- write_wrapper_inputs()
  on.exit(unlink(unlist(inp)), add = TRUE)

  res <- run_frankenplot_wrapper(c(
    "--jumble-csv", inp$csv, "--cns", inp$cns, "--output", inp$html,
    "--gis-output", tempfile(fileext = ".csv")
  ))
  expect_wrapper_error(res, "--gis-output requires --compute-gis")
})
