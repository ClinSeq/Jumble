test_that("classify_msi correctly classifies mononucleotide repeats", {
  # poly-A repeat: variant is "A" in a tract of AAAAAA... (>= 6 repeats total)
  # Left flank ends with 5 A's, variant is A, right flank starts with 2 A's
  # Total tract = 5 + 1 + 2 = 8 >= 6 threshold -> period 1
  result <- classify_msi(
    nonanchored_alt = "A",
    left_flank      = "TTTTTGCGCGAAAAAA",  # last 5 chars are A, but we need 6 contiguous from end
    right_flank     = "AATTTTT"
  )
  # Left: scanning from end of "TTTTTGCGCGAAAAAA" (16 chars), pos=16: A, 15: A, 14: A, 13: A, 12: A, 11: A = 6 A's
  # Var: 1 A

  # Right: pos=1: A, pos=2: A, pos=3: T = 2 A's
  # Total = 6 + 1 + 2 = 9 >= 6 -> period 1
  expect_equal(result, 1L)
})

test_that("classify_msi returns 0 for non-repeat context", {
  # Single A but no repeat tract around it
  result <- classify_msi(
    nonanchored_alt = "A",
    left_flank      = "TTTTTGCGCGTTTTTTT",
    right_flank     = "GCCCCCC"
  )
  expect_equal(result, 0L)
})

test_that("classify_msi correctly classifies dinucleotide repeats", {
  # CA repeat: variant is "CA" in a CA-repeat tract
  # Left flank ends with CACACACA (4 CA repeats)
  # Variant is CA (1 repeat)
  # Right flank starts with CACACA (3 CA repeats)
  # Total = 4 + 1 + 3 = 8 >= 4 threshold -> period 2
  result <- classify_msi(
    nonanchored_alt = "CA",
    left_flank      = "TTTTTCACACACA",
    right_flank     = "CACACACTTTT"
  )
  # Left: from end, pos=12-13: CA, 10-11: CA, 8-9: CA, 6-7: CA = 4
  # Wait, need to check: left_flank is "TTTTTCACACACA" (13 chars)
  # Scanning backward from end: pos=12-13: CA, 10-11: CA, 8-9: CA, 6-7: CA = 4
  # Right: "CACACACTTTT" -> pos 1-2: CA, 3-4: CA, 5-6: CA, 7-8: CT -> 3

  # Total = 4 + 1 + 3 = 8 >= 4 -> period 2
  expect_equal(result, 2L)
})

test_that("classify_msi correctly classifies trinucleotide repeats", {
  # CAG repeat
  result <- classify_msi(
    nonanchored_alt = "CAG",
    left_flank      = "TTTTCAGCAGCAGCAG",
    right_flank     = "CAGCAGCAGTTTT"
  )
  # Left: 16 chars, scanning from end: 14-16: CAG, 11-13: CAG, 8-10: CAG, 5-7: CAG = 4
  # Var: 1
  # Right: 1-3: CAG, 4-6: CAG, 7-9: CAG = 3
  # Total = 4 + 1 + 3 = 8 >= 4 -> period 3
  expect_equal(result, 3L)
})

test_that("classify_msi handles vectorized input", {
  result <- classify_msi(
    nonanchored_alt = c("A",                "CA",              "G"),
    left_flank      = c("AAAAAAAAAAAAAAAA", "TTTTTCACACACA",   "TTTTTTTTT"),
    right_flank     = c("AAAAAAAAAAAAAAAA", "CACACACTTTT",     "CCCCCCCC")
  )
  expect_equal(length(result), 3)
  expect_equal(result[1], 1L)  # mono repeat
  expect_equal(result[2], 2L)  # di repeat
  expect_equal(result[3], 0L)  # not a repeat
})

test_that("classify_msi handles NA inputs gracefully", {
  result <- classify_msi(
    nonanchored_alt = c("A", NA, "T"),
    left_flank      = c("AAAAAAAAAA", "AAAAAAA", NA),
    right_flank     = c("AAAAAAAAAA", "AAAAAAA", "AAAAAAA")
  )
  expect_equal(length(result), 3)
  expect_equal(result[2], 0L)  # NA alt -> 0
  expect_equal(result[3], 0L)  # NA flank -> 0
})

test_that("classify_msi shortest period wins", {
  # "AA" could be period 1 (A+A) or period 2 (AA)
  # If both qualify, period 1 should win
  result <- classify_msi(
    nonanchored_alt = "AA",
    left_flank      = "AAAAAAAAAAAA",
    right_flank     = "AAAAAAAAAAAA"
  )
  expect_equal(result, 1L)  # period 1 wins over period 2
})

test_that("classify_somatic_msi adds MSI column to somatic table", {
  # Create a minimal somatic-like data.table
  somatic <- data.table::data.table(
    chromosome = c("1", "1", "1"),
    start      = c(100, 200, 300),
    end        = c(100, 200, 300),
    REF        = c("A",  "AT",  "A"),
    ALT        = c("G",  "A",   "AG"),
    AF         = c(0.3,  0.2,   0.4),
    AO         = c(10,   8,     12),
    DP         = c(30,   40,    30)
  )

  # Run classification — skip if BSgenome not available
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")

  result <- classify_somatic_msi(somatic, "hg19")

  expect_true("MSI" %in% names(result))
  expect_equal(nrow(result), 3)

  # Row 1 is SNV -> MSI should be NA
  expect_true(is.na(result$MSI[1]))

  # Rows 2 and 3 are indels -> MSI should be integer (0, 1, 2, or 3)
  expect_true(!is.na(result$MSI[2]))
  expect_true(!is.na(result$MSI[3]))
  expect_true(result$MSI[2] %in% 0:3)
  expect_true(result$MSI[3] %in% 0:3)
})

test_that("classify_somatic_msi handles empty somatic table", {
  somatic <- data.table::data.table(
    chromosome = character(0),
    start      = integer(0),
    end        = integer(0),
    REF        = character(0),
    ALT        = character(0),
    AF         = numeric(0),
    AO         = integer(0),
    DP         = integer(0)
  )

  result <- classify_somatic_msi(somatic, "hg19")
  expect_true("MSI" %in% names(result))
  expect_equal(nrow(result), 0)
})

test_that("classify_somatic_msi handles SNV-only table", {
  somatic <- data.table::data.table(
    chromosome = c("1", "2"),
    start      = c(100, 200),
    end        = c(100, 200),
    REF        = c("A", "G"),
    ALT        = c("T", "C"),
    AF         = c(0.3, 0.5),
    AO         = c(10, 15),
    DP         = c(30, 30)
  )

  result <- classify_somatic_msi(somatic, "hg19")
  expect_true("MSI" %in% names(result))
  expect_true(all(is.na(result$MSI)))  # No indels -> all NA
})

test_that("compute_qc_metrics includes MSI columns when somatic provided", {
  # Minimal targets
  targets <- data.table::data.table(
    chromosome = c("1", "1"),
    start = c(100, 200),
    end   = c(200, 300),
    mid   = c(150, 250),
    count = c(100, 100),
    gc    = c(0.4, 0.5),
    log2  = c(0, 0),
    is_target = c(TRUE, TRUE)
  )

  # Somatic with MSI column and REF/ALT for total counts
  # AF column required by compute_qc_metrics (filters somatic[AF >= 0.05])
  somatic <- data.table::data.table(
    REF = c("A", "AT", "A", "A", "AC", "ACG"),
    ALT = c("G", "A",  "AT", "AT", "A",  "A"),
    MSI = c(NA, 0L, 1L, 1L, 2L, 3L),
    AF  = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
  )

  qc <- compute_qc_metrics(
    targets = targets,
    bam_file = "test.bam",
    reference_file = "ref.RDS",
    somatic = somatic
  )

  expect_true("somatic_snvs" %in% names(qc))
  expect_true("somatic_indels" %in% names(qc))
  expect_equal(qc$somatic_snvs, 1L)    # A>G is the only SNV
  expect_equal(qc$somatic_indels, 5L)   # 5 indels (different REF/ALT lengths)
  expect_equal(qc$MSI_mono, 2L)
  expect_equal(qc$MSI_di, 1L)
  expect_equal(qc$MSI_tri, 1L)
})

test_that("compute_qc_metrics always returns all 21 columns", {
  targets <- data.table::data.table(
    chromosome = c("1"),
    start = c(100),
    end   = c(200),
    mid   = c(150),
    count = c(100),
    gc    = c(0.4),
    log2  = c(0),
    is_target = c(TRUE)
  )

  qc <- compute_qc_metrics(
    targets = targets,
    bam_file = "test.bam",
    reference_file = "ref.RDS"
  )

  # All 21 columns should be present, even without VCFs
  expected_cols <- c("sample", "bam_file", "reference_file", "snp_vcf", "somatic_vcf",
                     "median_target_count", "gc_bias", "noise", "waviness",
                     "het_snps", "hom_snps", "sex", "contamination",
                     "somatic_snvs", "somatic_indels", "MSI_mono", "MSI_di", "MSI_tri",
                     "TMB_snv", "TMB_indel", "TMB_score")
  expect_equal(names(qc), expected_cols)

  # Somatic/SNP columns should be NA when not provided
  expect_true(is.na(qc$somatic_snvs))
  expect_true(is.na(qc$somatic_indels))
  expect_true(is.na(qc$het_snps))
  expect_true(is.na(qc$sex))
})

test_that("compute_snp_stats infers sex correctly", {
  # Female: normal het density on X (non-PAR positions: 50M+)
  snp_table <- data.table::data.table(
    chromosome = c(rep("1", 100), rep("X", 40)),
    start = c(seq(1e6, by = 1000, length.out = 100), seq(50e6, by = 1000, length.out = 40)),
    AD = c(rep(15, 140)),
    RD = c(rep(15, 140)),
    DP = c(rep(30, 140))
  )
  targets <- data.table::data.table(
    chromosome = c(rep("1", 500), rep("X", 200)),
    start = c(seq(1e6, by = 500, length.out = 500), seq(50e6, by = 500, length.out = 200)),
    end   = c(seq(1e6, by = 500, length.out = 500) + 499, seq(50e6, by = 500, length.out = 200) + 499),
    is_target = TRUE
  )
  result <- compute_snp_stats(snp_table, targets)
  expect_equal(result$sex, "female")

  # Male: no hets on X — only autosomal SNPs
  snp_table_male <- data.table::data.table(
    chromosome = c(rep("1", 100)),
    start = seq(1e6, by = 1000, length.out = 100),
    AD = c(rep(15, 100)),
    RD = c(rep(15, 100)),
    DP = c(rep(30, 100))
  )
  result_male <- compute_snp_stats(snp_table_male, targets)
  expect_equal(result_male$sex, "male")
})
