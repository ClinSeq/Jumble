test_that("get_chrom_arms returns correct structure for hg19 and hg38", {
    library(data.table)
    arms_hg19 <- Jumble:::get_chrom_arms("hg19")
    expect_true(is.data.table(arms_hg19))
    expect_equal(nrow(arms_hg19), 46) # 22 autosomes * 2 + X*2 = 46

    arms_hg38 <- get_chrom_arms("hg38")
    expect_true(is.data.table(arms_hg38))
    # Check an hg38 specific coordinate (e.g. chr1 p-arm end)
    # In our code: 123.4Mb
    chr1_p_end <- arms_hg38[chromosome == "1" & arm == "1p"]$end
    expect_equal(chr1_p_end, 123.4e6)
})
