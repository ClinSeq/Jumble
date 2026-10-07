test_that("reference annotation modes resolve legacy and slim references", {
  legacy_ref <- list(
    target_bed_file = "mock.bed",
    genome = "hg19",
    allgenes = data.table::data.table(
      `Gene name` = "GENE1",
      `Gene stable ID` = "ENSG000001",
      `Chromosome/scaffold name` = "1",
      `Gene start (bp)` = 100,
      `Gene end (bp)` = 200
    ),
    allexons = data.table::data.table(
      `Gene stable ID` = "ENSG000001",
      `Chromosome/scaffold name` = "1",
      `Exon region start (bp)` = 120,
      `Exon region end (bp)` = 180
    ),
    cancergenes_clinseq = data.table::data.table(
      hugo_symbol = "GENE1",
      ensembl_gene_id_version = "ENSG000001",
      ANNOT = "ONCO",
      chromosome = "1",
      start = 100,
      end = 200
    ),
    cytobands = data.table::data.table(
      chromosome = "1",
      start = 1,
      end = 300,
      band = "p36.33"
    )
  )

  normalized_legacy <- JumbleCNV:::normalize_reference_object(legacy_ref)
  expect_identical(normalized_legacy$reference_schema_version, "legacy-embedded")
  expect_identical(normalized_legacy$annotation_resources$mode, "embedded")
  expect_true(normalized_legacy$annotation_resources$has_gene_annotation)
  expect_true(normalized_legacy$annotation_resources$has_exon_annotation)
  expect_true(normalized_legacy$annotation_resources$has_cancer_gene_annotation)
  expect_true(normalized_legacy$annotation_resources$has_cytoband_annotation)

  slim_ref <- legacy_ref
  slim_ref$allgenes <- NULL
  slim_ref$allexons <- NULL
  slim_ref$cancergenes_clinseq <- NULL
  slim_ref$cytobands <- NULL

  normalized_slim <- JumbleCNV:::normalize_reference_object(slim_ref)
  expect_identical(normalized_slim$reference_schema_version, "slim")
  expect_identical(normalized_slim$annotation_resources$mode, "none")
  expect_false(normalized_slim$annotation_resources$has_gene_annotation)

  explicit_none <- JumbleCNV:::normalize_reference_object(legacy_ref, annotation_mode = "none")
  expect_identical(explicit_none$annotation_resources$mode, "none")
  expect_false(explicit_none$annotation_resources$has_gene_annotation)

  expect_error(
    JumbleCNV:::normalize_reference_object(slim_ref, annotation_mode = "embedded"),
    "contains no embedded annotation"
  )
})

test_that("target and segment annotation helpers are silent without annotation", {
  reference <- JumbleCNV:::normalize_reference_object(list(target_bed_file = "mock.bed"), annotation_mode = "none")

  targets <- data.table::data.table(
    chromosome = c("1", "1", "1"),
    start = c(1, 101, 201),
    end = c(100, 200, 300),
    mid = c(50, 150, 250),
    is_target = c(TRUE, TRUE, FALSE),
    type = c("target", "target", "bin")
  )

  expect_silent({
    annotated_targets <- JumbleCNV:::add_gene_annotations(data.table::copy(targets), reference)
    annotated_targets <- JumbleCNV:::add_cytoband_annotations(annotated_targets, reference)
    annotated_targets <- JumbleCNV:::add_exon_annotations(annotated_targets, reference)
  })

  expect_true("gene" %in% names(annotated_targets))
  expect_true("band" %in% names(annotated_targets))
  expect_true(all(annotated_targets$gene == ""))
  expect_true(all(annotated_targets$band == ""))
  expect_identical(annotated_targets[is_target == FALSE]$type, "background")

  segments <- data.table::data.table(
    chromosome = "1",
    start = 1,
    end = 3,
    nbrOfLoci = 3,
    mean = 0,
    start_pos = 1,
    end_pos = 300,
    genes = "",
    relevance = ""
  )

  expect_silent({
    annotated_segments <- JumbleCNV:::annotate_segments(
      data.table::copy(segments),
      cancergenes = NULL,
      cancerexons = NULL,
      allgenes = NULL,
      cytobands = NULL
    )
  })
  expect_identical(annotated_segments$genes, "")
  expect_identical(annotated_segments$relevance, "")
})

test_that("run_jumble accepts explicit no-annotation mode with existing fixtures", {
  testdata_dir <- system.file("testdata", package = "JumbleCNV")
  if (testdata_dir == "") testdata_dir <- "inst/testdata"

  ref_file <- file.path(testdata_dir, "gene_panel/reference.RDS")
  sample_file <- file.path(testdata_dir, "gene_panel/samples/test_sample_2.counts.RDS")

  skip_if(!file.exists(ref_file), "Reference file not found")
  skip_if(!file.exists(sample_file), "Sample file not found")

  out_dir <- tempfile()
  dir.create(out_dir)

  expect_warning(
    res <- run_jumble(
      bam_file = sample_file,
      reference_file = ref_file,
      output_dir = out_dir,
      snp_vcf = NULL,
      annotation_mode = "none"
    ),
    NA
  )

  expect_type(res, "list")
  expect_true("targets" %in% names(res))
  expect_true("segments" %in% names(res))
  expect_true(file.exists(file.path(out_dir, "test_sample_2.jumble.csv")))
})

test_that("run_jumble supports slim references without embedded annotation", {
  testdata_dir <- system.file("testdata", package = "JumbleCNV")
  if (testdata_dir == "") testdata_dir <- "inst/testdata"

  ref_file <- file.path(testdata_dir, "gene_panel/reference.RDS")
  sample_file <- file.path(testdata_dir, "gene_panel/samples/test_sample_2.counts.RDS")

  skip_if(!file.exists(ref_file), "Reference file not found")
  skip_if(!file.exists(sample_file), "Sample file not found")

  slim_ref <- readRDS(ref_file)
  slim_ref$allgenes <- NULL
  slim_ref$allexons <- NULL
  slim_ref$cancergenes_clinseq <- NULL
  slim_ref$cytobands <- NULL
  slim_ref$reference_schema_version <- "slim"

  out_dir <- tempfile()
  dir.create(out_dir)
  slim_ref_file <- file.path(out_dir, "slim.reference.RDS")
  saveRDS(slim_ref, slim_ref_file)

  expect_warning(
    res <- run_jumble(
      bam_file = sample_file,
      reference_file = slim_ref_file,
      output_dir = out_dir,
      snp_vcf = NULL,
      annotation_mode = "auto"
    ),
    NA
  )

  expect_type(res, "list")
  expect_true("targets" %in% names(res))
  expect_true("segments" %in% names(res))
  expect_true("gene" %in% names(res$targets))
  expect_true(all(res$targets$gene == "" | is.na(res$targets$gene)))
  expect_true(file.exists(file.path(out_dir, "test_sample_2.jumble.csv")))
})
