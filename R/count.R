#' Generate Read Counts from BAM File
#'
#' Counts reads in target bins (or WGS bins) from a BAM file.
#'
#' @param bam_file Path to the input BAM file.
#' @param target_bed Path to a BED file defining target regions. If NULL, WGS
#'   mode is used.
#' @param wgs_bin_size Bin size for WGS mode (default: 10000).
#' @param format Output format. Currently returns a list object.
#'
#' @return A list containing counts, ranges, and metadata.
#' @importFrom bamsignals bamCount
#' @importFrom GenomicRanges makeGRangesFromDataFrame seqnames seqinfo width
#'   start<- end<- reduce gaps width<- sort
#' @importFrom GenomeInfoDb seqlengths seqlevelsStyle<-
#' @importFrom Rsamtools BamFile
#' @importFrom data.table fread data.table rbindlist as.data.table :=
#' @importFrom stringr str_detect str_remove
#' @importFrom methods is
#' @export
generate_counts <- function(bam_file, target_bed = NULL,
              wgs_bin_size = 10000) {
  # 1. Validate Inputs -------------------------------------------------------
  # Check if BAM file exists
  if (!file.exists(bam_file)) stop("BAM file not found: ", bam_file)

  bam_obj <- BamFile(bam_file)

  # Check chromosome style (UCSC vs Ensembl)
  seq_info <- seqinfo(bam_obj)
  bam_is_ucsc <- any(str_detect(seqnames(seq_info), "^chr"))

  if (bam_is_ucsc) {
    chromlength <- seqlengths(seq_info)[paste0("chr", c(1:22, "X", "Y"))]
    names(chromlength) <- str_remove(names(chromlength), "^chr")
  } else {
    chromlength <- seqlengths(seq_info)[c(1:22, "X", "Y")]
  }

  # 2. Define Ranges ---------------------------------------------------------
  # Determine mode
  wgs <- is.null(target_bed)

  bed <- NULL
  ranges <- NULL

  if (!wgs) {
    # Targeted Mode
    if (!file.exists(target_bed)) {
      stop("Target BED file not found: ", target_bed)
    }

    # Detect if BED has header
    first_line <- readLines(target_bed, n = 1)
    has_header <- grepl("^[a-zA-Z]", first_line) && !grepl("^chr[0-9XYM]", first_line, ignore.case = TRUE)
    
    bed <- fread(target_bed, header = has_header)
    
    # Ensure we have at least 3 columns
    if (ncol(bed) < 3) stop("BED file must have at least 3 columns (chr, start, end)")
    
    # Standardize column names for first 3 columns
    setnames(bed, 1:3, c("chromosome", "start", "end"))

    # Format BED data
    # Add target ID as first column
    targets_original <- cbind(
        data.table(target = 1:nrow(bed)),
        bed[, .(chromosome, start, end)]
    )
    
    # Ensure Ensembl style for processing: Remove 'chr' prefix and force to character
    targets_original[, chromosome := str_remove(as.character(chromosome), "^chr")]

    # --- Target Bins ---
    binsize <- 200
    ranges_target <- makeGRangesFromDataFrame(targets_original)

    # Extend targets
    start(ranges_target) <- start(ranges_target) - 50
    end(ranges_target) <- end(ranges_target) + 50

    # Extend again if too short
    widths <- width(ranges_target)
    if (any(widths < 200)) {
      start(ranges_target) <- start(ranges_target) - 50
      end(ranges_target) <- end(ranges_target) + 50
    }

    # Merge overlaps
    ranges_target <- reduce(ranges_target)

    # Shrink to multiple of binsize
    w <- width(ranges_target)
    mid <- start(ranges_target) + round(w / 2)
    new_width <- w - w %% binsize
    new_start <- mid - new_width / 2
    new_end <- mid + new_width / 2

    # Split longer targets
    # We need to do this per chromosome/range to handle vectorization correctly
    # or use a loop
    # The original code looped over chromosomes of the ranges.
    # Let's replicate the logic but safely.

    # Create a data.table to expand
    dt_ranges <- data.table(
      chromosome = as.character(seqnames(ranges_target)),
      start_orig = new_start,
      end_orig = new_end
    )

    targets_list <- lapply(1:nrow(dt_ranges), function(i) {
      s <- dt_ranges$start_orig[i]
      e <- dt_ranges$end_orig[i]
      if (e < s) {
        return(NULL)
      } # Should not happen if width >= binsize

      starts <- seq(s, e - binsize, by = binsize) + 1
      ends <- seq(s + binsize, e, by = binsize)

      data.table(
        chromosome = dt_ranges$chromosome[i],
        start = starts,
        end = ends,
        mid = ceiling((starts + ends) / 2)
      )
    })
    targets <- rbindlist(targets_list)
    targets[, length := end - start]

    target_ranges <- makeGRangesFromDataFrame(targets)

    # --- Background Bins ---
    # Dummy targets to ensure coverage of all chromosomes
    dummy_targets <- rbind(
      data.table(
        target = 0, chromosome = c(1:22, "X", "Y"), start = 1,
        end = 10
      ),
      data.table(
        target = 0, chromosome = c(1:22, "X", "Y"),
        start = chromlength - 10, end = chromlength
      )
    )

    bg_binsize <- 1e6
    bg_minsize <- 3e5

    # Combine original targets and dummy targets for background calculation
    combined_targets <- rbind(
      targets_original[, .(target,
        chromosome = as.character(chromosome),
        start, end
      )],
      dummy_targets
    )
    ranges_bg_base <- makeGRangesFromDataFrame(combined_targets)

    # Extend by 1k
    start(ranges_bg_base) <- start(ranges_bg_base) - 1000
    end(ranges_bg_base) <- end(ranges_bg_base) + 1000

    ranges_bg_base <- reduce(ranges_bg_base)

    # Invert to get background
    # We need to define seqlengths for gaps() to work correctly on the whole
    # genome
    # But gaps() on GRanges works per sequence if seqlengths are set.
    # However, we are working with Ensembl style names here (1..22, X, Y)
    # chromlength is already named 1..22, X, Y

    # We need to construct a GRanges with proper seqlengths for gaps
    sl <- chromlength
    # Filter to only standard chromosomes
    ranges_bg_base <- ranges_bg_base[seqnames(ranges_bg_base) %in% names(sl)]
    seqlevels(ranges_bg_base) <- names(sl)
    seqlengths(ranges_bg_base) <- sl

    ranges_bg <- gaps(ranges_bg_base)
    # gaps returns the whole genome minus the ranges.
    # It might include strand specific gaps, but we ignore strand (*).

    # Drop too short
    ranges_bg <- ranges_bg[width(ranges_bg) > bg_minsize]

    # Shrink to multiple of binsize
    w_bg <- width(ranges_bg)
    mid_bg <- start(ranges_bg) + round(w_bg / 2)
    new_width_bg <- w_bg
    # Only shrink if larger than binsize? Original code:
    # new_width[width > binsize] <- width[width > binsize] -
    #   width[width > binsize] %% binsize
    mask <- w_bg > bg_binsize
    new_width_bg[mask] <- w_bg[mask] - w_bg[mask] %% bg_binsize

    new_start_bg <- mid_bg - round(new_width_bg / 2)
    new_end_bg <- mid_bg + round(new_width_bg / 2)

    # Split longer antitargets
    dt_bg <- data.table(
      chromosome = as.character(seqnames(ranges_bg)),
      start_orig = new_start_bg,
      end_orig = new_end_bg,
      n_bins = ceiling(new_width_bg / bg_binsize)
    )

    antitargets_list <- lapply(1:nrow(dt_bg), function(i) {
      s <- dt_bg$start_orig[i]
      e <- dt_bg$end_orig[i]
      n <- dt_bg$n_bins[i]
      chrom <- dt_bg$chromosome[i]

      if (n > 1) {
        starts <- seq(s, e - bg_binsize, by = bg_binsize) + 1
        ends <- seq(s + bg_binsize, e, by = bg_binsize)
        data.table(chromosome = chrom, start = starts, end = ends)
      } else {
        data.table(chromosome = chrom, start = s + 1, end = e)
      }
    })
    background <- rbindlist(antitargets_list)
    background[, mid := (start + end) / 2]
    background[, length := end - start]

    background_ranges <- makeGRangesFromDataFrame(background)

    # Combine
    ranges <- sort(c(target_ranges, background_ranges))
  } else {
    # WGS Mode
    bins_list <- lapply(names(chromlength), function(chrom) {
      len <- chromlength[chrom]
      starts <- seq(0, len - wgs_bin_size * 2, by = wgs_bin_size) + 1
      data.table(
        chromosome = chrom,
        start = starts,
        end = starts + wgs_bin_size - 1
      )
    })
    bins <- rbindlist(bins_list)
    ranges <- makeGRangesFromDataFrame(bins)
  }

  # Set style to match BAM
  if (bam_is_ucsc) {
    seqlevelsStyle(ranges) <- "UCSC"
  } else {
    seqlevelsStyle(ranges) <- "Ensembl"
    # Note: 'Ensembl' style is 1, 2, ... which matches our construction
  }

  # 3. Count Reads -----------------------------------------------------------
  # QC fail, optical duplicate, supplementary alignment, secondary alignment,
  # unmapped
  flag <- 2816

  # bamsignals::bamCount
  counts_all <- bamCount(bam_file, ranges,
    paired.end = "midpoint", mapq = 20,
    filteredFlag = flag, verbose = FALSE
  )
  counts_short <- bamCount(bam_file, ranges,
    paired.end = "midpoint", mapq = 20,
    filteredFlag = flag, tlenFilter = c(0, 150), verbose = FALSE
  )
  counts_medium <- bamCount(bam_file, ranges,
    paired.end = "midpoint",
    mapq = 20, filteredFlag = flag, tlenFilter = c(0, 300), verbose = FALSE
  )

  # 4. Construct Output ------------------------------------------------------
  res <- list(
    input_bam_file = bam_file,
    target_bed_file = if (wgs) NULL else target_bed,
    date_count = date(),
    bed = bed,
    ranges = ranges,
    chromlength = chromlength,
    genome = detect_genome(bam_file),
    flag = flag,
    count = counts_all,
    count_short = counts_short,
    count_medium = counts_medium
  )

  return(res)
}
