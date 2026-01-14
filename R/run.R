#' Run Jumble Analysis
#'
#' Runs the complete Jumble analysis pipeline on a query sample.
#'
#' @param bam_file Path to the input BAM file.
#' @param reference_file Path to the reference RDS file.
#' @param output_dir Directory to save results.
#' @param snp_vcf Optional path to a VCF file with SNPs.
#' @return A list containing results (targets, segments, etc.).
#' @importFrom data.table fread fwrite
#' @importFrom ggplot2 ggsave
#' @export
run_jumble <- function(bam_file, reference_file, output_dir = ".",
                       snp_vcf = NULL, somatic_vcf = NULL, cores = 1, genome = NULL, ...) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # 1. Load Reference --------------------------------------------------------
  # 1. Load Reference --------------------------------------------------------
  if (is.character(reference_file)) {
    message("Loading reference: ", reference_file)
    reference <- readRDS(reference_file)
  } else {
    reference <- reference_file
  }

  # 2. Reference PCA ---------------------------------------------------------
  message("Computing reference PCA...")
  ref_pca <- compute_reference_pca(reference)

  # 3. Counts ----------------------------------------------------------------
  message("Generating counts for query sample...")
  # Check if input is BAM or Counts
  if (grepl("\\.RDS$", bam_file, ignore.case = TRUE)) {
    counts <- readRDS(bam_file)
    bam_path <- counts$input_bam_file

    # Ensure counts match reference order/template
    # If the counts RDS was generated with a different order (e.g. unsorted), we must reorder it
    if (!is.null(counts$ranges) && !is.null(reference$ranges)) {
      # Create GRanges for matching
      # Assuming identical binning, just different order

      # Check if identical first (quick check)
      if (!identical(counts$ranges, reference$ranges)) {
        message("Reordering input counts to match reference template...")

        # Use findOverlaps to map input ranges to reference ranges
        # We expect 1:1 mapping
        ol <- GenomicRanges::findOverlaps(reference$ranges, counts$ranges, type = "equal")

        # Check for completeness
        if (length(ol) != length(reference$ranges)) {
          warning("Input counts RDS does not fully match reference ranges. Some bins may be missing or mismatched.")
        }

        # Get the index of the matching input bin for each reference bin
        # subjectHits(ol) is indices in counts$ranges
        # queryHits(ol) is indices in reference$ranges

        # We want reference 1..N. map[i] = input_index_for_ref_i
        map <- rep(NA, length(reference$ranges))
        map[S4Vectors::queryHits(ol)] <- S4Vectors::subjectHits(ol)

        if (any(is.na(map))) {
          stop("Fatal: Input counts RDS ranges do not cover all reference ranges.")
        }

        # Reorder counts
        if (!is.null(counts$count)) counts$count <- counts$count[map]
        if (!is.null(counts$count_short)) counts$count_short <- counts$count_short[map]
        if (!is.null(counts$count_medium)) counts$count_medium <- counts$count_medium[map]
        counts$ranges <- counts$ranges[map] # Should now be identical to reference
      }
    }
  } else {
    # Use reference target bed
    target_bed <- reference$target_bed_file
    if (target_bed == "wgs") target_bed <- NULL

    counts <- list()
    counts$input_bam_file <- bam_file
    counts$ranges <- reference$ranges
    counts$chromlength <- reference$chromlength

    flag <- if (!is.null(reference$flag)) reference$flag else 2816

    message("Counting reads...")
    counts$count <- bamsignals::bamCount(bam_file, reference$ranges,
      paired.end = "midpoint",
      mapq = 20, filteredFlag = flag, verbose = FALSE
    )
    counts$count_short <- bamsignals::bamCount(bam_file, reference$ranges,
      paired.end = "midpoint",
      mapq = 20, filteredFlag = flag, tlenFilter = c(0, 150),
      verbose = FALSE
    )

    # Save counts
    saveRDS(counts, file.path(output_dir, paste0(
      basename(bam_file),
      ".counts.RDS"
    )))
  }

  # 4. Prepare Targets -------------------------------------------------------
  targets <- data.table::copy(reference$target_template)
  targets[, chromosome := clean_chrom_names(chromosome)]
  targets[, chromosome := clean_chrom_names(chromosome)]
  sample_name <- basename(bam_file)
  targets[, sample := sample_name]
  targets[, count := counts$count]
  targets[, count_short := counts$count_short]

  # Add count_all column if exists in counts
  if (!is.null(counts$count_all)) {
    targets[, count_all := counts$count_all]
  }

  # Annotate tiled bins (bins that are close together)
  targets[, left := c(Inf, abs(diff(mid)))]
  targets[, right := c(abs(diff(mid)), Inf)]
  targets[, farthest := left]
  targets[right > left, farthest := right]
  targets[, is_tiled := FALSE]
  targets[farthest < 250, is_tiled := TRUE]
  targets[, c("left", "right", "farthest") := NULL]

  # Set tiled type if enough tiled bins
  if (sum(targets$is_tiled == TRUE) > 500) {
    targets[is_tiled == TRUE, type := "tiled"]
  }

  # 5. Annotations (Genes/Exons) ---------------------------------------------
  if (!is.null(reference$allgenes) && nrow(reference$allgenes) > 0) {
    allgenes <- reference$allgenes
    # Use columns: 'Gene name', 'Chromosome/scaffold name', 'Gene start (bp)', 'Gene end (bp)'

    # makeGRanges for bins
    binranges <- GenomicRanges::makeGRangesFromDataFrame(targets,
      seqnames.field = "chromosome",
      start.field = "start", end.field = "end"
    )

    # makeGRanges for genes
    generanges <- GenomicRanges::makeGRangesFromDataFrame(allgenes,
      seqnames.field = "Chromosome/scaffold name",
      start.field = "Gene start (bp)", end.field = "Gene end (bp)",
      keep.extra.columns = TRUE # Keep Gene name
    )

    # Find overlaps
    ov <- GenomicRanges::findOverlaps(binranges, generanges)

    # Map genes to bins
    # Logic: A bin might overlap multiple genes. We want to list them.

    # Create mapping table
    mapping <- data.table::data.table(
      bin_id = S4Vectors::queryHits(ov),
      gene_symbol = allgenes$`Gene name`[S4Vectors::subjectHits(ov)]
    )

    # Aggregate genes per bin (comma separated)
    # Use efficient data.table aggregation
    mapping_agg <- mapping[, .(gene = paste(unique(gene_symbol), collapse = ",")), by = bin_id]

    # Assign to targets
    targets[mapping_agg$bin_id, gene := mapping_agg$gene]
  }

  # Annotate exonic bins (overlap with exons)
  if (!is.null(reference$allexons) && nrow(reference$allexons) > 0) {
    allexons <- reference$allexons
    # Prepare exon ranges
    exons_dt <- unique(allexons[, .(
      id = `Gene stable ID`,
      chromosome = `Chromosome/scaffold name`,
      start = `Exon region start (bp)`,
      end = `Exon region end (bp)`
    )])

    # Only proceed if chromosomes match format
    if (nrow(exons_dt) > 0) {
      binranges <- GenomicRanges::makeGRangesFromDataFrame(targets,
        seqnames.field = "chromosome",
        start.field = "start",
        end.field = "end"
      )
      exonranges <- GenomicRanges::makeGRangesFromDataFrame(exons_dt,
        seqnames.field = "chromosome",
        start.field = "start",
        end.field = "end"
      )

      exon_overlap <- data.table::as.data.table(
        GenomicRanges::findOverlaps(binranges, exonranges)
      )
      targets[unique(exon_overlap$queryHits), type := "exonic"]
    }
  }

  # Ensure background bins are labeled correctly
  targets[is_target == FALSE, type := "background"]

  # Define backbone
  targets <- define_backbone(targets, reference)

  # Save all targets before any processing
  alltargets <- data.table::copy(targets)

  # 6. Process SNPs ----------------------------------------------------------
  snp_table <- NULL
  if (!is.null(snp_vcf) && file.exists(snp_vcf)) {
    message("Processing SNPs from VCF...")
    snp_table <- process_snps(snp_vcf, targets, sample_name = sample_name)

    if (!is.null(snp_table)) {
      # SNP table will be saved at the end with correct sample name
    }
  }

  # 7. Normalize -------------------------------------------------------------
  message("Normalizing...")
  targets <- normalize_sample(targets, ref_pca)

  # Check if WGS
  is_wgs <- FALSE
  if (!is.null(reference$target_bed_file)) {
    if (grepl("wgs", reference$target_bed_file, ignore.case = TRUE)) {
      is_wgs <- TRUE
    }
  }

  # Adjust X chromosome log2 ratios to align with background
  # This correction matches the original script's logic (lines 744-760).
  # It calculates the median difference between background (bg) and target (tg)
  # bins in 1MB windows on the X chromosome and shifts the target bins by this
  # difference.
  # This helps normalize X chromosome coverage relative to the rest of the genome.
  if (!is_wgs) {
    try(
      {
        # Filter for X chromosome bins with valid log2 values
        temp <- targets[chromosome == "X" & !is.na(log2)]

        # Define 1MB windows
        temp[, pos_1M := 1 * round(start / 1e6)]

        # Calculate median log2 for background and target bins within each window
        temp[, bg := as.numeric(NA)][, bg := median(log2[is_target == FALSE],
          na.rm = TRUE
        ), by = pos_1M]
        temp[, tg := as.numeric(NA)][, tg := median(log2[is_target == TRUE &
          is_tiled == FALSE], na.rm = TRUE), by = pos_1M]

        # Reduce to unique windows
        temp <- unique(temp[, .(pos_1M, bg, tg)])

        # Smooth the medians using a running median of width 11
        temp[, bg_median := stats::runmed(bg, 11, na.action = "na.omit")]
        temp[, tg_median := stats::runmed(tg, 11, na.action = "na.omit")]

        # Calculate the correction factor (median difference)
        temp[, dif := bg_median - tg_median]
        x_correct <- median(temp$dif, na.rm = TRUE)

        # Apply correction to X chromosome target bins
        targets[chromosome == "X" & is_target == TRUE, log2 := log2 + x_correct]
      },
      silent = TRUE
    )
  }

  # 8. Segment ---------------------------------------------------------------
  message("Segmenting...")

  # Determine alpha parameter for segmentation
  # Alpha controls the significance level for change point detection.
  # WGS data uses a stricter alpha (1e-5) to avoid over-segmentation due to
  # higher density.
  alpha <- 0.02
  if (reference$target_bed_file == "wgs") alpha <- 1e-5

  # Prepare cancer genes and exons from reference
  cancergenes <- if (!is.null(reference$cancergenes_clinseq)) {
    reference$cancergenes_clinseq
  } else {
    NULL
  }

  cancerexons <- if (!is.null(reference$allexons)) {
    reference$allexons
  } else {
    NULL
  }

  seg_res <- segment_data(targets,
    alpha = alpha,
    cancergenes = cancergenes,
    cancerexons = cancerexons
  )
  targets <- seg_res$targets
  segments <- seg_res$segments

  # Merge back filtered bins (matching original script line 852)
  # Segmentation may filter out some bins (e.g., those with missing data).
  # We merge back with 'alltargets' to ensure the final output contains all
  # original bins, preserving the complete genomic structure.
  common_cols <- intersect(names(alltargets), names(targets))
  targets <- merge(alltargets, targets, by = common_cols, all = TRUE)
  targets <- targets[order(bin)]

  # Determine sample name for output
  sample_name <- basename(bam_file)
  sample_name <- sub("\\.counts\\.RDS$", "", sample_name, ignore.case = TRUE)
  sample_name <- sub("\\.bam$", "", sample_name, ignore.case = TRUE)

  # Add smooth_log2 for plotting
  targets[, smooth_log2 := stats::runmed(log2, k = 7), by = chromosome]

  # 9. Compute GIS -----------------------------------------------------------
  # Compute GIS (only if SNPs are present)
  gis_table <- NULL
  if (!is.null(snp_table)) {
    message("Computing GIS...")
    genome_version <- if (!is.null(reference$genome)) reference$genome else "hg19"
    gis_table <- compute_gis_table(targets, snp_table, genome = genome_version)

    # Map MAF to targets for plotting (plot_gis_score requires targets$maf)
    # Filter SNPs to valid allele ratio range
    snps_sub <- snp_table[allele_ratio > .02 & allele_ratio < .98]
    snps_sub[, maf := 0.5 + abs(allele_ratio - 0.5)]

    # Create GRanges
    targetranges <- GenomicRanges::makeGRangesFromDataFrame(targets,
      seqnames.field = "chromosome",
      start.field = "start",
      end.field = "end",
      ignore.strand = TRUE
    )
    snpranges <- GenomicRanges::makeGRangesFromDataFrame(snps_sub,
      seqnames.field = "chromosome",
      start.field = "start",
      end.field = "end",
      ignore.strand = TRUE
    )

    # Overlap
    hits <- GenomicRanges::findOverlaps(snpranges, targetranges)

    # Assign (if multiple SNPs map to one bin, this takes the last one, or random.
    # ideally we would average, but matching gis.R logic for consistency)
    targets[, maf := as.numeric(NA)]
    targets[S4Vectors::subjectHits(hits), maf := snps_sub[S4Vectors::queryHits(hits)]$maf]
  } else {
    message("Skipping GIS (no SNPs)...")
    targets[, maf := as.numeric(NA)] # Initialize NA for consistency
  }


  # 9b. Process Somatic VCF ----------------------------------------------------
  somatic <- NULL
  # 9.5 Process Somatic VCF --------------------------------------------------
  somatic <- NULL
  if (!is.null(somatic_vcf)) {
      # Determine genome version: Argument > Reference > Default
      use_genome <- if (!is.null(genome)) genome else if (!is.null(reference$genome)) reference$genome else "hg19"
      somatic <- process_somatic_vcf(somatic_vcf, reference, genome = use_genome)
  }

  # 10. Plot -----------------------------------------------------------------
  message("Plotting...")

  # Ensure targets are fully annotated before plotting (label, selected_genes)
  annotate_targets(targets)

  plot_file <- file.path(output_dir, paste0(sample_name, ".png"))
  plot_results(targets, segments,
    reference = reference,
    output_file = plot_file, title = sample_name,
    snp_table = snp_table,
    gis_table = gis_table,
    somatic_table = somatic
  )

  # Plot GIS separately if available
  if (!is.null(gis_table)) {
    gis_plot_file <- file.path(output_dir, paste0(sample_name, ".gis.png"))
    tryCatch(
      {
        plot_gis_score(gis_table, targets, gis_plot_file, title = sample_name)
      },
      error = function(e) {
        warning("Failed to create GIS plot: ", conditionMessage(e))
      }
    )
  }

  # 11. Save Results ---------------------------------------------------------
  message("Saving results...")

  # Main CSV file
  data.table::fwrite(targets, file.path(output_dir, paste0(
    sample_name,
    ".jumble.csv"
  )), sep = ",")

  # Save SNP table if present
  if (!is.null(snp_table)) {
    snp_output_file <- file.path(output_dir, paste0(
      sample_name,
      ".jumble_snps.csv"
    ))
    data.table::fwrite(snp_table, snp_output_file, sep = ",")
  }

  # Estimate Contamination
  contamination <- if (!is.null(snp_table)) {
    estimate_contamination(snp_table)
  } else {
    NA_real_
  }

  # Compute and save QC metrics
  qc_metrics <- compute_qc_metrics(
    targets = targets,
    bam_file = bam_file,
    reference_file = reference_file,
    snp_vcf = snp_vcf,
    sample_name = sample_name,
    contamination = contamination
  )

  qc_output_file <- file.path(output_dir, paste0(sample_name, ".qc.csv"))
  write_qc_metrics(qc_metrics, qc_output_file)

  # Save GIS Table
  if (!is.null(gis_table)) {
    gis_output_file <- file.path(output_dir, paste0(sample_name, ".jumble_gis.csv"))
    data.table::fwrite(gis_table, gis_output_file, sep = ",")
  }


  # Export CNR (copy number ratios)
  # Format: chromosome, start, end, gene, depth, log2, weight, gc, count, type
  cnr <- targets[!is.na(log2), .(
    chromosome = as.character(chromosome),
    start, end, gene,
    depth = count,
    log2,
    weight = 1,
    gc, count, type
  )]
  cnr[gene == "", gene := "-"]
  data.table::fwrite(cnr, file.path(output_dir, paste0(sample_name, ".cnr")),
    sep = "\t"
  )

  # Export CNS (copy number segments)
  # Format: chromosome, start, end, gene, log2, depth, probes, relevance
  if (nrow(segments) > 0) {
    # Add nbrOfLoci (number of bins in segment)
    segments[, nbrOfLoci := as.integer(NA)]
    for (i in seq_len(nrow(segments))) {
      segments$nbrOfLoci[i] <- sum(targets$segment == i, na.rm = TRUE)
    }

    cns <- segments[, .(
      chromosome,
      start = start_pos,
      end = end_pos,
      gene = ifelse(is.na(genes) | genes == "", "-", genes),
      log2 = round(mean, 3),
      depth = round(2^mean, 3),
      probes = nbrOfLoci,
      relevance = ifelse(is.na(relevance), "", relevance)
    )]
    data.table::fwrite(cns, file.path(output_dir, paste0(sample_name, ".cns")),
      sep = "\t"
    )
  }

  # Export SEG (DNAcopy format)
  # Format: ID, chrom, loc.start, loc.end, num.mark, seg.mean, C
  if (nrow(segments) > 0) {
    seg <- segments[, .(
      ID = sample_name,
      chrom = chromosome,
      loc.start = start_pos,
      loc.end = end_pos,
      num.mark = nbrOfLoci,
      seg.mean = mean,
      C = as.numeric(NA)
    )]
    # Convert X/Y to numeric for DNAcopy
    seg[, chrom := stringr::str_replace(chrom, "Y", "24")]
    seg[, chrom := stringr::str_replace(chrom, "X", "23")]
    seg[, chrom := as.numeric(chrom)]
    data.table::fwrite(seg, file.path(output_dir, paste0(
      sample_name,
      "_dnacopy.seg"
    )), sep = "\t")
  }

  message("Done.")
  invisible(list(targets = targets, segments = segments))
}
