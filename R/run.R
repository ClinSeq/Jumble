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
    counts <- sanitize_legacy_counts(counts)
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
  counts$count_medium <- bamsignals::bamCount(bam_file, reference$ranges,
                                              paired.end = "midpoint",
                                              mapq = 20, filteredFlag = flag, tlenFilter = c(0, 300),
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
#' @param exclude_long_fragments If TRUE, use count_medium (≤300bp) instead of
#'   count (all fragments) as the main depth signal.
#' @return Targets data.table
#' @keywords internal
prepare_targets <- function(reference, counts, bam_file, exclude_long_fragments = FALSE) {
  targets <- data.table::copy(reference$target_template)
  targets[, chromosome := clean_chrom_names(chromosome)]
  
  sample_name <- basename(bam_file)
  targets[, sample := sample_name]
  
  if (exclude_long_fragments && !is.null(counts$count_medium)) {
    targets[, count := counts$count_medium]
  } else {
    targets[, count := counts$count]
  }
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

#' Add Cytoband Annotations to Targets
#'
#' @param targets Targets data.table
#' @param reference Reference object
#' @return Targets with cytoband annotations
#' @keywords internal
add_cytoband_annotations <- function(targets, reference) {
  targets[, band := ""]
  if (!is.null(reference$cytobands) && nrow(reference$cytobands) > 0) {
    binranges <- GenomicRanges::makeGRangesFromDataFrame(targets,
                                                         seqnames.field = "chromosome",
                                                         start.field = "start", end.field = "end"
    )
    cyto_gr <- GenomicRanges::makeGRangesFromDataFrame(reference$cytobands,
                                                       seqnames.field = "chromosome",
                                                       start.field = "start", end.field = "end",
                                                       keep.extra.columns = TRUE
    )
    
    ov <- GenomicRanges::findOverlaps(binranges, cyto_gr)
    
    mapping <- data.table::data.table(
      bin_id = S4Vectors::queryHits(ov),
      band = reference$cytobands$band[S4Vectors::subjectHits(ov)]
    )
    mapping <- mapping[!is.na(band) & band != ""]
    
    if (nrow(mapping) > 0) {
      agg <- mapping[, .(band_coll = {
        bands <- unique(band)
        if (length(bands) == 1) bands else paste0(bands[1], "-", bands[length(bands)])
      }), by = bin_id]
      targets[agg$bin_id, band := paste0(targets$chromosome[agg$bin_id], agg$band_coll)]
    }
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
  
  alpha <- 0.01
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
                          cancerexons = cancerexons,
                          allgenes = reference$allgenes,
                          cytobands = reference$cytobands
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
#' @param hrd_model Optional HRD model
#' @return List with gis_table and updated targets
#' @keywords internal
compute_gis_and_maf <- function(targets, snp_table, reference, hrd_model = NULL) {
  if (is.null(snp_table)) {
    targets[, maf := as.numeric(NA)]
    targets[, long_median := as.numeric(NA)]
    targets[, local_snp_bg := as.numeric(NA)]
    return(list(gis_table = NULL, targets = targets))
  }
  
  message("Computing GIS...")
  genome_version <- if (!is.null(reference$genome)) reference$genome else "hg19"
  gis_table <- compute_gis_table(targets, snp_table, genome = genome_version, hrd_model = hrd_model)
  
  # Map long_median from 5MB background segments to targets
  bins_final <- attr(gis_table, "bins_final")
  targets[, long_median := as.numeric(NA)]
  targets[, local_snp_bg := as.numeric(NA)]
  if (!is.null(bins_final) && nrow(bins_final) > 0) {
    br <- GenomicRanges::makeGRangesFromDataFrame(bins_final, seqnames.field = "chromosome", start.field = "start", end.field = "end")
    tr <- GenomicRanges::makeGRangesFromDataFrame(targets, seqnames.field = "chromosome", start.field = "start", end.field = "end")
    near_idx <- GenomicRanges::nearest(tr, br, select = "arbitrary", ignore.strand = TRUE)
    
    valid <- !is.na(near_idx)
    if (any(valid)) {
      targets[valid, long_median := bins_final$long_median[near_idx[valid]]]
      targets[valid, local_snp_bg := bins_final$maf[near_idx[valid]]]
    }
  }
  
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
  
  # Compute local_snp_bg for somatic rare SNP filtering
  # 1. Smooth within segment
  targets[, local_snp_bg := as.numeric(NA)]
  
  safe_runmed <- function(x) {
    valid <- !is.na(x)
    n <- sum(valid)
    if (n == 0) return(x)
    if (n <= 2) {
      x[valid] <- median(x[valid])
      return(x)
    }
    k <- min(9, n)
    if (k %% 2 == 0) k <- k - 1
    x[valid] <- stats::runmed(x[valid], k)
    return(x)
  }
  
  targets[, local_snp_bg := safe_runmed(maf), by = segment]
  
  targets[, local_snp_bg := {
    if (sum(!is.na(local_snp_bg)) >= 2) {
      stats::approx(1:.N, local_snp_bg, xout = 1:.N, rule = 2)$y
    } else if (sum(!is.na(local_snp_bg)) == 1) {
      rep(mean(local_snp_bg, na.rm = TRUE), .N)
    } else {
      rep(as.numeric(NA), .N)
    }
  }, by = segment]
  
  # 2. Extrapolate missing segments strictly from structurally adjacent valid logR boundaries
  for (chr in unique(targets$chromosome)) {
    chr_idx <- which(targets$chromosome == chr)
    seg_ids <- unique(targets$segment[chr_idx])
    if (length(seg_ids) == 0) next
    
    seg_summary <- lapply(seg_ids, function(s) {
      idx <- which(targets$segment == s & targets$chromosome == chr)
      list(segment = s, log2 = median(targets$smooth_log2[idx], na.rm = TRUE), bg = median(targets$local_snp_bg[idx], na.rm = TRUE))
    })
    seg_dt <- data.table::rbindlist(seg_summary)
    
    valid_indices <- which(!is.na(seg_dt$bg))
    if (length(valid_indices) == 0) {
      targets[chromosome == chr, local_snp_bg := 0.5]
    } else {
      for (i in which(is.na(seg_dt$bg))) {
        left_idx <- valid_indices[valid_indices < i]
        left_idx <- if (length(left_idx) > 0) max(left_idx) else NA
        
        right_idx <- valid_indices[valid_indices > i]
        right_idx <- if (length(right_idx) > 0) min(right_idx) else NA
        
        if (is.na(left_idx)) {
          best_bg <- seg_dt$bg[right_idx]
        } else if (is.na(right_idx)) {
          best_bg <- seg_dt$bg[left_idx]
        } else {
          dist_left <- abs(seg_dt$log2[i] - seg_dt$log2[left_idx])
          dist_right <- abs(seg_dt$log2[i] - seg_dt$log2[right_idx])
          if (is.na(dist_left) && is.na(dist_right)) {
            best_bg <- seg_dt$bg[left_idx]
          } else if (is.na(dist_left)) {
            best_bg <- seg_dt$bg[right_idx]
          } else if (is.na(dist_right)) {
            best_bg <- seg_dt$bg[left_idx]
          } else {
            best_bg <- if (dist_left <= dist_right) seg_dt$bg[left_idx] else seg_dt$bg[right_idx]
          }
        }
        targets[chromosome == chr & segment == seg_dt$segment[i], local_snp_bg := best_bg]
      }
    }
  }
  
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
                          sample_name,
                          qc_metrics$gc_bias[[1]],
                          qc_metrics$noise[[1]],
                          qc_metrics$waviness[[1]])
    # Add MSI score (mono + di) if available
    msi_mono <- qc_metrics$MSI_mono[[1]]
    msi_di   <- qc_metrics$MSI_di[[1]]
    if (!is.null(msi_mono) && !is.na(msi_mono) && !is.null(msi_di) && !is.na(msi_di)) {
      msi_score <- msi_mono + msi_di
      plot_title <- sprintf("%s | Repeat Tract Indels: %d", plot_title, msi_score)
    }
    # Add TMB score if available and valid
    tmb_score <- qc_metrics$TMB_score[[1]]
    if (!is.null(tmb_score) && !is.na(tmb_score) && nchar(tmb_score) > 0 && tmb_score != "NA") {
      plot_title <- sprintf("%s | TMB: %s", plot_title, tmb_score)
    }
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
  
  # TMB VAF plot
  if (!is.null(somatic)) {
    tmb_plot_file <- file.path(output_dir, paste0(sample_name, ".tmb.png"))
    tryCatch(
      {
        plot_tmb_vaf(somatic, targets, tmb_plot_file, title = sample_name)
      },
      error = function(e) {
        warning("Failed to create TMB plot: ", conditionMessage(e))
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

#' Export Analysis Files (CNR, CNS, SEG, GENES)
#'
#' @param targets Targets data.table
#' @param segments Segments data.table
#' @param output_dir Output directory
#' @param sample_name Sample name
#' @param reference Reference object (needed for annotation)
#' @param snp_table SNP table (needed for annotation)
#' @keywords internal
export_analysis_files <- function(targets, segments, output_dir, sample_name, reference, snp_table) {
  # Export CNR
  cnr <- targets[!is.na(log2), .(
    chromosome = as.character(chromosome),
    start, end, gene,
    bed_name = if ("bed_name" %in% names(targets)) bed_name else "",
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
      band = if ("band" %in% names(segments)) ifelse(is.na(band) | band == "", "-", band) else "-",
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
  
  # Export GENES
  tgt <- targets[is_target == TRUE]
  if (!is.null(reference$allgenes) && nrow(tgt) > 0) {
    # Map segmented CBS mean to bins via coordinate overlap (nearest for edge cases)
    tgt[, seg_mean := as.numeric(NA)]
    seg_gr <- GenomicRanges::makeGRangesFromDataFrame(segments,
                seqnames.field = "chromosome", start.field = "start_pos", end.field = "end_pos")
    tgt_gr <- GenomicRanges::makeGRangesFromDataFrame(tgt,
                seqnames.field = "chromosome", start.field = "start", end.field = "end")
    seg_ov <- GenomicRanges::nearest(tgt_gr, seg_gr, select = "arbitrary", ignore.strand = TRUE)
    valid_seg <- !is.na(seg_ov)
    if (any(valid_seg)) {
      tgt[valid_seg, seg_mean := segments$mean[seg_ov[valid_seg]]]
    }
                 
    generanges <- GenomicRanges::makeGRangesFromDataFrame(reference$allgenes,
                                                          seqnames.field = "Chromosome/scaffold name",
                                                          start.field = "Gene start (bp)", end.field = "Gene end (bp)")
    tranges <- GenomicRanges::makeGRangesFromDataFrame(tgt,
                                                       seqnames.field = "chromosome",
                                                       start.field = "start", end.field = "end")
    ov <- data.table::as.data.table(GenomicRanges::findOverlaps(generanges, tranges))
    
    if (nrow(ov) > 0) {
      mapped <- tgt[ov$subjectHits]
      mapped$gene_symbol <- reference$allgenes$`Gene name`[ov$queryHits]
      mapped$gene_start <- reference$allgenes$`Gene start (bp)`[ov$queryHits]
      mapped$gene_end <- reference$allgenes$`Gene end (bp)`[ov$queryHits]
      
      safe_median <- function(x) if (length(x) == 0 || all(is.na(x))) as.numeric(NA) else median(x, na.rm = TRUE)
      safe_min <- function(x) if (length(x) == 0 || all(is.na(x))) as.numeric(NA) else min(x, na.rm = TRUE)
      safe_max <- function(x) if (length(x) == 0 || all(is.na(x))) as.numeric(NA) else max(x, na.rm = TRUE)
      
      gene_metrics <- mapped[, .(
        chromosome = chromosome[1],
        start = as.numeric(gene_start[1]),
        end = as.numeric(gene_end[1]),
        band = paste(unique(band[band != "" & band != "-"]), collapse = "-"),
        n_bins = as.integer(.N),
        exonic_bins = as.integer(sum(type == "exonic", na.rm = TRUE)),
        depth = as.integer(round(safe_median(count))),
        segments = as.integer(data.table::uniqueN(segment)),
        segmented_median = as.numeric(round(safe_median(2^seg_mean[if (any(type == "exonic", na.rm=TRUE)) type == "exonic" else TRUE]), 3)),
        segmented_min = as.numeric(round(safe_min(2^seg_mean[if (any(type == "exonic", na.rm=TRUE)) type == "exonic" else TRUE]), 3)),
        segmented_max = as.numeric(round(safe_max(2^seg_mean[if (any(type == "exonic", na.rm=TRUE)) type == "exonic" else TRUE]), 3)),
        segmented_armlevel = as.numeric(round(safe_median(2^long_median), 3)),
        smooth_maf = as.numeric(round(safe_median(local_snp_bg), 3))
      ), by = gene_symbol]
      data.table::setnames(gene_metrics, "gene_symbol", "gene")
      
      gene_metrics <- gene_metrics[gene != "" & !is.na(gene)]
      
      # Natively compute SNP metrics for gene table
      if (!is.null(snp_table) && nrow(snp_table) > 0) {
        gene_gr <- GenomicRanges::makeGRangesFromDataFrame(gene_metrics, seqnames.field = "chromosome", start.field = "start", end.field = "end")
        snps_sub <- snp_table[allele_ratio > .02 & allele_ratio < .98]
        snps_sub[, maf := 0.5 + abs(allele_ratio - 0.5)]
        snp_gr <- GenomicRanges::makeGRangesFromDataFrame(snps_sub, seqnames.field = "chromosome", start.field = "start", end.field = "end", ignore.strand=TRUE)
        
        snp_ov <- data.table::as.data.table(GenomicRanges::findOverlaps(gene_gr, snp_gr))
        if (nrow(snp_ov) > 0) {
          snp_merge <- data.table::data.table(g_idx = snp_ov$queryHits, maf = snps_sub$maf[snp_ov$subjectHits])
          snp_agg <- snp_merge[, .(snps = .N, maf = round(median(maf, na.rm = TRUE), 3)), by = g_idx]
          gene_metrics[, snps := 0]
          gene_metrics[, maf := as.numeric(NA)]
          gene_metrics[snp_agg$g_idx, snps := snp_agg$snps]
          gene_metrics[snp_agg$g_idx, maf := snp_agg$maf]
        } else {
          gene_metrics[, snps := 0]
          gene_metrics[, maf := as.numeric(NA)]
        }
      } else {
        gene_metrics[, snps := 0]
        gene_metrics[, maf := as.numeric(NA)]
      }
      
      # Order strictly by genomic location (factorizing chrom to avoid stringsort 1, 10, 2)
      chrom_levels <- c(as.character(1:22), "X", "Y")
      gene_metrics[, chr_factor := factor(gsub("chr", "", chromosome), levels = chrom_levels)]
      gene_metrics <- gene_metrics[order(chr_factor, start)]
      gene_metrics[, chr_factor := NULL]
      
      data.table::setcolorder(gene_metrics, c("gene", "chromosome", "start", "end", "band", "n_bins", "exonic_bins", "depth", "segments", "segmented_median", "segmented_min", "segmented_max", "segmented_armlevel", "snps", "maf", "smooth_maf"))
      
      data.table::fwrite(gene_metrics, file.path(output_dir, paste0(sample_name, ".genes.csv")), sep = ",")
    }
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
#' @param hrd_model_file Optional path to a custom HRD model object (e.g. randomForest, glm, or function) saved as an RDS file. When supplied, a supplementary Custom HRD score is added to GIS output.
#' @param exclude_long_fragments If TRUE, use count_medium (fragments ≤300bp)
#'   instead of count (all fragments) as the main depth signal. Set to TRUE when
#'   using clipoverlap BAMs where the midpoint calculation may be affected by
#'   TLEN inflation from overlap clipping.
#' @param ... Additional arguments.
#' @return A list containing results (targets, segments, etc.).
#' @importFrom data.table fread fwrite
#' @importFrom ggplot2 ggsave
#' @export
run_jumble <- function(bam_file, reference_file, output_dir = ".",
                       snp_vcf = NULL, somatic_vcf = NULL, cores = 1,
                       genome = NULL, correction = "optim", hrd_model_file = NULL,
                       exclude_long_fragments = FALSE, ...) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  hrd_model <- NULL
  if (!is.null(hrd_model_file)) {
    if (file.exists(hrd_model_file)) {
      hrd_model <- readRDS(hrd_model_file)
    } else {
      warning("Provided hrd_model_file does not exist. Ignoring custom HRD model.")
    }
  }
  
  # 1. Load Reference
  reference <- load_reference_data(reference_file)
  
  # 2. Load or Generate Counts
  message("Generating counts for query sample...")
  counts <- load_or_generate_counts(bam_file, reference, output_dir)
  
  # 3. Leave-Me-Out: Check if test sample exists in reference by counting data
  if (!is.null(reference$allcounts)) {
    match_idx <- integer(0)
    for (i in seq_along(reference$allcounts)) {
      if (isTRUE(all.equal(counts$count, reference$allcounts[[i]]$count)) &&
          isTRUE(all.equal(counts$count_short, reference$allcounts[[i]]$count_short))) {
        match_idx <- c(match_idx, i)
      }
    }
    if (length(match_idx) > 0) {
      remaining <- length(reference$allcounts) - length(match_idx)
      if (remaining >= 1) {
        message("Test sample matched in reference. Excluding it before PCA...")
        reference$allcounts <- reference$allcounts[-match_idx]
        reference$samples <- reference$samples[-match_idx]
      } else {
        warning("Test sample matched all reference samples. Skipping leave-me-out exclusion to preserve reference for PCA.")
      }
    }
  }

  if (exclude_long_fragments) {
    message("Using count_medium (≤300bp fragments) as main depth signal.")
  }
  
  # 4. Reference PCA
  message("Computing reference PCA...")
  ref_pca <- compute_reference_pca(reference, exclude_long_fragments = exclude_long_fragments)

  # 5. Prepare Targets
  targets <- prepare_targets(reference, counts, bam_file, exclude_long_fragments = exclude_long_fragments)
  
  # 5. Add Annotations
  targets <- add_gene_annotations(targets, reference)
  targets <- add_cytoband_annotations(targets, reference)
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
  gis_result <- compute_gis_and_maf(targets, snp_table, reference, hrd_model = hrd_model)
  gis_table <- gis_result$gis_table
  targets <- gis_result$targets
  
  # 9b. Process Somatic VCF
  somatic <- NULL
  if (!is.null(somatic_vcf)) {
    use_genome <- if (!is.null(genome)) genome else if (!is.null(reference$genome)) reference$genome else "hg19"
    somatic <- process_somatic_vcf(somatic_vcf, reference, genome = use_genome)
    
    # Rare germline SNP rejection flagging (LOH integration)
    if (!is.null(somatic) && nrow(somatic) > 0 && "local_snp_bg" %in% names(targets)) {
      somatic[, somatic_maf := 0.5 + abs(AF - 0.5)]
      somatic[, background_maf := targets$local_snp_bg[bin]]
      somatic[, is_indel := nchar(REF) != nchar(ALT)]
      somatic[, maf_diff := abs(somatic_maf - background_maf)]
      somatic[, is_rare_snp := !is.na(maf_diff) & maf_diff <= 0.05]
    }
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
    somatic_vcf = somatic_vcf,
    sample_name = sample_name,
    contamination = contamination,
    snp_table = snp_table,
    somatic = somatic
  )
  
  # 11. Plot
  generate_plots(targets, segments, reference, output_dir, sample_name,
                 snp_table, gis_table, somatic, qc_metrics)
  
  # 12. Save Results
  save_analysis_results(targets, segments, snp_table, gis_table, qc_metrics,
                        output_dir, sample_name, bam_file,
                        reference_file, snp_vcf)
  
  export_analysis_files(targets, segments, output_dir, sample_name, reference, snp_table)
  
  message("Done.")
  invisible(list(targets = targets, segments = segments))
}