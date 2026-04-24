test_that("load_cytobands works", {
  cb19 <- load_cytobands("hg19")
  expect_s3_class(cb19, "data.table")
  expect_true(nrow(cb19) > 400)
  expect_true("chromosome" %in% names(cb19))
  expect_true("band" %in% names(cb19))

  cb38 <- load_cytobands("hg38")
  expect_s3_class(cb38, "data.table")
  expect_true(nrow(cb38) > 400)
})

test_that("load_target_bed preserves bed_name", {
  mock_bed_path <- tempfile(fileext = ".bed")
  bed_content <- data.frame(
    chromosome = c("1", "2"),
    start = c(10000, 50000),
    end = c(20000, 60000),
    name = c("RegionA", "RegionB")
  )
  write.table(bed_content, mock_bed_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  bed <- load_target_bed(mock_bed_path)
  
  expect_true("bed_name" %in% names(bed))
  expect_equal(bed$bed_name, c("RegionA", "RegionB"))
})
