#' Load Reference Data
#'
#' @param reference_file Path to reference RDS file or reference object
#' @return Reference object
#' @keywords internal
load_reference_data <- function(reference_file) {
  if (is.character(reference_file)) {
    message("Loading reference: ", reference_file)
    reference <- readRDS(reference_file)
  } else {
    reference <- reference_file
  }
  return(reference)
}

#' Load or Generate Counts
#'
#' @param bam_file Path to BAM file or counts RDS
#' @param reference Reference object
#' @param output_dir Output directory
#' @return Counts list
#' @keywords internal
load_or_generate_counts <- function(bam_file, reference, output_dir) {
  # Check if input is pre-computed counts
  if (grepl("\\.RDS$", bam_file, ignore.case = TRUE)) {
    counts <- readRDS(bam_file)
    counts <- reorder_counts_to_reference(counts, reference)
  } else {
    counts <- generate_counts_from_bam(bam_file, reference, output_dir)
  }
  return(counts)
}

#' Reorder Counts to Match Reference
#'
#' @param counts Counts object
#' @param reference Reference object
#' @return Reordered counts
#' @keywords internal
reorder_counts_to_reference <- function(counts, reference) {
  if (!is.null(counts$ranges) && !is.null(reference$ranges)) {
    if (!identical(counts$ranges, reference$ranges)) {
      message("Reordering input counts to match reference template...")
      
      ol <- GenomicRanges::findOverlaps(reference$ranges, counts$ranges, type = "equal")
      
      if (length(ol) != length(reference$ranges)) {
        warning("Input counts RDS does not fully match reference ranges. Some bins may be missing or mismatched.")
      }
      
      map <- rep(NA, length(reference$ranges))
      map[S4Vectors::queryHits(ol)] <- S4Vectors::subjectHits(ol)
      
      if (any(is.na(map))) {
        stop("Fatal: Input counts RDS ranges do not cover all reference ranges.")
      }
      
      # Reorder counts
      if (!is.null(counts$count)) counts$count <- counts$count[map]
      if (!is.null(counts$count_short)) counts$count_short <- counts$count_short[map]
      if (!is.null(counts$count_medium)) counts$count_medium <- counts$count_medium[map]
      counts$ranges <- counts$ranges[map]
    }
  }
  return(counts)
}

#' Generate Counts from BAM File
#'
#' @param bam_file Path to BAM file
#' @param reference Reference object
#' @param output_dir Output directory
#' @return Counts list
#' @keywords internal
generate_counts_from_bam <- function(bam_file, reference, output_dir) {
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
    basename(bam_file), ".counts.RDS"
  )))
  
  return(counts)
}

#' Prepare Targets Data Table
#'
#' @param reference Reference object
#' @param counts Counts object
#' @param bam_file BAM file path
#' @return Targets data.table
#' @keywords internal
prepare_targets <- function(reference, counts, bam_file) {
  targets <- data.table::copy(reference$target_template)
  targets[, chromosome := clean_chrom_names(chromosome)]
  
  sample_name <- basename(bam_file)
  targets[, sample := sample_name]
  targets[, count := counts$count]
  targets[, count_short := counts$count_short]
  
  if (!is.null(counts$count_all)) {
    targets[, count_all := counts$count_all]
  }
  
  targets <- annotate_tiled_bins(targets)
  
  return(targets)
}

#' Annotate Tiled Bins
#'
#' @param targets Targets data.table
#' @return Targets with tiled annotations
#' @keywords internal
annotate_tiled_bins <- function(targets) {
  targets[, left := c(Inf, abs(diff(mid)))]
  targets[, right := c(abs(diff(mid)), Inf)]
  targets[, farthest := left]
  targets[right > left, farthest := right]
  targets[, is_tiled := FALSE]
  targets[farthest < 250, is_tiled := TRUE]
  targets[, c("left", "right", "farthest") := NULL]
  
  if (sum(targets$is_tiled == TRUE) > 500) {
    targets[is_tiled == TRUE, type := "tiled"]
  }
  
  return(targets)
}

#' Add Gene Annotations to Targets
#'
#' @param targets Targets data.table
#' @param reference Reference object
#' @return Targets with gene annotations
#' @keywords internal
add_gene_annotations <- function(targets, reference) {
  if (!is.null(reference$allgenes) && nrow(reference$allgenes) > 0) {
    allgenes <- reference$allgenes
    
    binranges <- GenomicRanges::makeGRangesFromDataFrame(targets,
                                                         seqnames.field = "chromosome",
                                                         start.field = "start", end.field = "end"
    )
    
    generanges <- GenomicRanges::makeGRangesFromDataFrame(allgenes,
                                                          seqnames.field = "Chromosome/scaffold name",
                                                          start.field = "Gene start (bp)", end.field = "Gene end (bp)",
                                                          keep.extra.columns = TRUE
    )
    
    ov <- GenomicRanges::findOverlaps(binranges, generanges)
    
    mapping <- data.table::data.table(
      bin_id = S4Vectors::queryHits(ov),
      gene_symbol = allgenes$`Gene name`[S4Vectors::subjectHits(ov)]
    )
    
    mapping_agg <- mapping[, .(gene = paste(unique(gene_symbol), collapse = ",")), by = bin_id]
    targets[mapping_agg$bin_id, gene := mapping_agg$gene]
  }
  
  return(targets)
}

#' Add Exon Annotations to Targets
#'
#' @param targets Targets data.table
#' @param reference Reference object
#' @return Targets with exon annotations
#' @keywords internal
add_exon_annotations <- function(targets, reference) {
  if (!is.null(reference$allexons) && nrow(reference$allexons) > 0) {
    allexons <- reference$allexons
    
    exons_dt <- unique(allexons[, .(
      id = `Gene stable ID`,
      chromosome = `Chromosome/scaffold name`,
      start = `Exon region start (bp)`,
      end = `Exon region end (bp)`
    )])
    
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
  
  targets[is_target == FALSE, type := "background"]
  
  return(targets)
}

#' Apply X Chromosome Correction
#'
#' @param targets Targets data.table
#' @param is_wgs Logical indicating if WGS data
#' @return Corrected targets
#' @keywords internal
apply_x_chromosome_correction <- function(targets, is_wgs) {
  if (!is_wgs) {
    try(
      {
        temp <- targets[chromosome == "X" & !is.na(log2)]
        temp[, pos_1M := 1 * round(start / 1e6)]
        temp[, bg := as.numeric(NA)][, bg := median(log2[is_target == FALSE],
                                                    na.rm = TRUE
        ), by = pos_1M]
        temp[, tg := as.numeric(NA)][, tg := median(log2[is_target == TRUE &
                                                           is_tiled == FALSE], na.rm = TRUE), by = pos_1M]
        
        temp <- unique(temp[, .(pos_1M, bg, tg)])
        temp[, bg_median := stats::runmed(bg, 11, na.action = "na.omit")]
        temp[, tg_median := stats::runmed(tg, 11, na.action = "na.omit")]
        temp[, dif := bg_median - tg_median]
        x_correct <- median(temp$dif, na.rm = TRUE)
        
        targets[chromosome == "X" & is_target == TRUE, log2 := log2 + x_correct]
      },
      silent = TRUE
    )
  }
  
  return(targets)
}

#' Perform Segmentation
#'
#' @param targets Targets data.table
#' @param reference Reference object
#' @param alltargets Original targets before filtering
#' @return List with targets and segments
#' @keywords internal
perform_segmentation <- function(targets, reference, alltargets) {
  message("Segmenting...")
  
  alpha <- 0.02
  if (reference$target_bed_file == "wgs") alpha <- 1e-5
  
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
  
  # Merge back filtered bins
  common_cols <- intersect(names(alltargets), names(targets))
  targets <- merge(alltargets, targets, by = common_cols, all = TRUE)
  targets <- targets[order(bin)]
  
  targets[, smooth_log2 := stats::runmed(log2, k = 7), by = chromosome]
  
  return(list(targets = targets, segments = segments))
}

#' Compute GIS and Add MAF to Targets
#'
#' @param targets Targets data.table
#' @param snp_table SNP data.table
#' @param reference Reference object
#' @return List with gis_table and updated targets
#' @keywords internal
compute_gis_and_maf <- function(targets, snp_table, reference) {
  if (is.null(snp_table)) {
    targets[, maf := as.numeric(NA)]
    return(list(gis_table = NULL, targets = targets))
  }
  
  message("Computing GIS...")
  genome_version <- if (!is.null(reference$genome)) reference$genome else "hg19"
  gis_table <- compute_gis_table(targets, snp_table, genome = genome_version)
  
  # Map MAF to targets
  snps_sub <- snp_table[allele_ratio > .02 & allele_ratio < .98]
  snps_sub[, maf := 0.5 + abs(allele_ratio - 0.5)]
  
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
  
  hits <- GenomicRanges::findOverlaps(snpranges, targetranges)
  
  targets[, maf := as.numeric(NA)]
  targets[S4Vectors::subjectHits(hits), maf := snps_sub[S4Vectors::queryHits(hits)]$maf]
  
  return(list(gis_table = gis_table, targets = targets))
}

#' Generate and Save Plots
#'
#' @param targets Targets data.table
#' @param segments Segments data.table
#' @param reference Reference object
#' @param output_dir Output directory
#' @param sample_name Sample name
#' @param snp_table SNP table
#' @param gis_table GIS table
#' @param somatic Somatic variants table
#' @keywords internal
generate_plots <- function(targets, segments, reference, output_dir, 
                           sample_name, snp_table, gis_table, somatic, qc_metrics = NULL) {
  message("Plotting...")
  
  annotate_targets(targets)
  
  plot_title <- sample_name
  if (!is.null(qc_metrics)) {
    plot_title <- sprintf("%s | GC-Bias: %.2f | Noise: %.2f | Waviness: %.2f", 
                          sample_name, qc_metrics$gc_bias, qc_metrics$noise, qc_metrics$waviness)
  }
  
  plot_file <- file.path(output_dir, paste0(sample_name, ".png"))
  plot_results(targets, segments,
               reference = reference,
               output_file = plot_file, title = plot_title,
               snp_table = snp_table,
               gis_table = gis_table,
               somatic_table = somatic
  )
  
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
  
  # MSI VAF plot
  if (!is.null(somatic) && "MSI" %in% names(somatic)) {
    msi_plot_file <- file.path(output_dir, paste0(sample_name, ".msi.png"))
    tryCatch(
      {
        plot_msi_vaf(somatic, msi_plot_file, title = sample_name)
      },
      error = function(e) {
        warning("Failed to create MSI plot: ", conditionMessage(e))
      }
    )
  }
}

#' Save Analysis Results
#'
#' @param targets Targets data.table
#' @param segments Segments data.table
#' @param snp_table SNP table
#' @param gis_table GIS table
#' @param qc_metrics QC metrics table
#' @param output_dir Output directory
#' @param sample_name Sample name
#' @param bam_file BAM file path
#' @param reference_file Reference file path
#' @param snp_vcf SNP VCF path
#' @keywords internal
save_analysis_results <- function(targets, segments, snp_table, gis_table, qc_metrics,
                                  output_dir, sample_name, bam_file,
                                  reference_file, snp_vcf) {
  message("Saving results...")
  
  # Main CSV file
  data.table::fwrite(targets, file.path(output_dir, paste0(
    sample_name, ".jumble.csv"
  )), sep = ",")
  
  # Save SNP table
  if (!is.null(snp_table)) {
    snp_output_file <- file.path(output_dir, paste0(
      sample_name, ".jumble_snps.csv"
    ))
    data.table::fwrite(snp_table, snp_output_file, sep = ",")
  }
  
  # Save QC metrics
  qc_output_file <- file.path(output_dir, paste0(sample_name, ".qc.csv"))
  write_qc_metrics(qc_metrics, qc_output_file)
  
  # Save GIS table
  if (!is.null(gis_table)) {
    gis_output_file <- file.path(output_dir, paste0(sample_name, ".jumble_gis.csv"))
    data.table::fwrite(gis_table, gis_output_file, sep = ",")
  }
}

#' Export Analysis Files (CNR, CNS, SEG)
#'
#' @param targets Targets data.table
#' @param segments Segments data.table
#' @param output_dir Output directory
#' @param sample_name Sample name
#' @keywords internal
export_analysis_files <- function(targets, segments, output_dir, sample_name) {
  # Export CNR
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
  
  # Export CNS
  if (nrow(segments) > 0) {
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
  
  # Export SEG
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
    seg[, chrom := stringr::str_replace(chrom, "Y", "24")]
    seg[, chrom := stringr::str_replace(chrom, "X", "23")]
    seg[, chrom := as.numeric(chrom)]
    data.table::fwrite(seg, file.path(output_dir, paste0(
      sample_name, "_dnacopy.seg"
    )), sep = "\t")
  }
}

#' Run Jumble Analysis
#'
#' Runs the complete Jumble analysis pipeline on a query sample.
#'
#' @param bam_file Path to the input BAM file.
#' @param reference_file Path to the reference RDS file.
#' @param output_dir Directory to save results.
#' @param snp_vcf Optional path to a VCF file with SNPs.
#' @param somatic_vcf Optional path to a VCF file with somatic variants.
#' @param cores Number of cores (currently unused).
#' @param genome Genome version.
#' @param correction String indicating the method: "optim" (L1+TV, default) or "rlm" (Robust LM).
#' @param ... Additional arguments.
#' @return A list containing results (targets, segments, etc.).
#' @importFrom data.table fread fwrite
#' @importFrom ggplot2 ggsave
#' @export
run_jumble <- function(bam_file, reference_file, output_dir = ".",
                       snp_vcf = NULL, somatic_vcf = NULL, cores = 1, 
                       genome = NULL, correction = "optim", ...) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 1. Load Reference
  reference <- load_reference_data(reference_file)
  
  # 2. Reference PCA
  message("Computing reference PCA...")
  ref_pca <- compute_reference_pca(reference)
  
  # 3. Load or Generate Counts
  message("Generating counts for query sample...")
  counts <- load_or_generate_counts(bam_file, reference, output_dir)
  
  # 4. Prepare Targets
  targets <- prepare_targets(reference, counts, bam_file)
  
  # 5. Add Annotations
  targets <- add_gene_annotations(targets, reference)
  targets <- add_exon_annotations(targets, reference)
  targets <- define_backbone(targets)
  
  # Save all targets before processing
  alltargets <- data.table::copy(targets)
  
  # 6. Process SNPs
  snp_table <- NULL
  sample_name <- basename(bam_file)
  if (!is.null(snp_vcf) && file.exists(snp_vcf)) {
    message("Processing SNPs from VCF...")
    snp_table <- process_snps(snp_vcf, targets, sample_name = sample_name)
  }
  
  # 7. Normalize
  message("Normalizing...")
  targets <- normalize_sample(targets, ref_pca, correction)
  
  # Check if WGS
  is_wgs <- FALSE
  if (!is.null(reference$target_bed_file)) {
    if (grepl("wgs", reference$target_bed_file, ignore.case = TRUE)) {
      is_wgs <- TRUE
    }
  }
  
  # Apply X chromosome correction
  targets <- apply_x_chromosome_correction(targets, is_wgs)
  
  # 8. Segment
  seg_result <- perform_segmentation(targets, reference, alltargets)
  targets <- seg_result$targets
  segments <- seg_result$segments
  
  # Determine sample name
  sample_name <- basename(bam_file)
  sample_name <- sub("\\.counts\\.RDS$", "", sample_name, ignore.case = TRUE)
  sample_name <- sub("\\.bam$", "", sample_name, ignore.case = TRUE)
  
  # 9. Compute GIS
  gis_result <- compute_gis_and_maf(targets, snp_table, reference)
  gis_table <- gis_result$gis_table
  targets <- gis_result$targets
  
  # 9b. Process Somatic VCF
  somatic <- NULL
  if (!is.null(somatic_vcf)) {
    use_genome <- if (!is.null(genome)) genome else if (!is.null(reference$genome)) reference$genome else "hg19"
    somatic <- process_somatic_vcf(somatic_vcf, reference, genome = use_genome)
  }
  
  # 10. Compute QC Metrics
  contamination <- if (!is.null(snp_table)) {
    estimate_contamination(snp_table)
  } else {
    NA_real_
  }
  
  qc_metrics <- compute_qc_metrics(
    targets = targets,
    bam_file = bam_file,
    reference_file = reference_file,
    snp_vcf = snp_vcf,
    sample_name = sample_name,
    contamination = contamination,
    somatic = somatic
  )
  
  # 11. Plot
  generate_plots(targets, segments, reference, output_dir, sample_name,
                 snp_table, gis_table, somatic, qc_metrics)
  
  # 12. Save Results
  save_analysis_results(targets, segments, snp_table, gis_table, qc_metrics,
                        output_dir, sample_name, bam_file,
                        reference_file, snp_vcf)
  
  export_analysis_files(targets, segments, output_dir, sample_name)
  
  message("Done.")
  invisible(list(targets = targets, segments = segments))
}