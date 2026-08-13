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
  expect_match(res$stdout, "--somatic-vcf")
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
