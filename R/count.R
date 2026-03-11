#' Generate Read Counts from BAM File
#'
#' Counts reads in target bins (or WGS bins) from a BAM file.
#'
#' @param bam_file Path to the input BAM file.
#' @param target_bed Path to a BED file defining target regions. If NULL, WGS
#'   mode is used.
#' @param wgs_bin_size Bin size for WGS mode (default: 10000).
#'
#' @return A list containing counts, ranges, and metadata.
#' @importFrom Rsamtools BamFile
#' @importFrom GenomeInfoDb seqinfo seqnames seqlevelsStyle<-
#' @importFrom GenomicRanges sort
#' @importFrom GenomeInfoDb seqlengths
#' @export
generate_counts <- function(bam_file, target_bed = NULL, wgs_bin_size = 10000) {
  
  # 1. Validation and Setup
  if (!file.exists(bam_file)) stop("BAM file not found: ", bam_file)
  
  # Get chromosome lengths and style from BAM
  bam_info <- get_bam_info(bam_file)
  chromlength <- bam_info$chromlength
  bam_is_ucsc <- bam_info$is_ucsc
  
  # 2. Define Ranges
  if (!is.null(target_bed)) {
    # --- Targeted Mode ---
    if (!file.exists(target_bed)) stop("Target BED file not found: ", target_bed)
    
    # Load and standardize BED
    raw_targets <- load_target_bed(target_bed)
    
    # Generate Bins
    target_ranges <- create_target_bins(raw_targets)
    bg_ranges <- create_background_bins(raw_targets, chromlength)
    
    # Combine and Sort
    ranges <- sort(c(target_ranges, bg_ranges))
    
    # Pass input BED back for metadata
    bed_out <- raw_targets
    
  } else {
    # --- WGS Mode ---
    ranges <- create_wgs_bins(chromlength, wgs_bin_size)
    bed_out <- NULL
  }
  
  # 3. Match BAM Style
  # We processed in Ensembl (1, 2...) style; convert to UCSC if BAM is UCSC
  if (bam_is_ucsc) {
    seqlevelsStyle(ranges) <- "UCSC"
  } else {
    seqlevelsStyle(ranges) <- "Ensembl"
  }
  
  # 4. Count Reads
  counts <- perform_bin_counts(bam_file, ranges)
  
  # 5. Construct Output
  list(
    input_bam_file = bam_file,
    target_bed_file = target_bed,
    date_count = date(),
    bed = bed_out,
    ranges = ranges,
    chromlength = chromlength,
    genome = detect_genome(bam_file), # Assuming this helper exists in your package
    flag = counts$flag,
    count = counts$all,
    count_short = counts$short,
    count_medium = counts$medium
  )
}

# ==============================================================================
# Internal Helper Functions
# ==============================================================================

#' Extract Chromosome Lengths and Style from BAM
#' @keywords internal
#' @importFrom Rsamtools BamFile
#' @importFrom GenomicRanges seqinfo seqnames
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom stringr str_detect str_remove
get_bam_info <- function(bam_file) {
  bam_obj <- BamFile(bam_file)
  seq_info <- seqinfo(bam_obj)
  
  # Check style
  is_ucsc <- any(str_detect(seqnames(seq_info), "^chr"))
  
  # Extract lengths for standard chromosomes (1-22, X, Y)
  # We standardize names to Ensembl style (no "chr") for internal processing
  if (is_ucsc) {
    # Extract chr1...chrY
    desired_names <- paste0("chr", c(1:22, "X", "Y"))
    # Handle cases where BAM might miss a chromosome gracefully, or subset strictly
    existing_names <- intersect(desired_names, seqnames(seq_info))
    len_vec <- seqlengths(seq_info)[existing_names]
    names(len_vec) <- str_remove(names(len_vec), "^chr")
  } else {
    desired_names <- c(1:22, "X", "Y")
    existing_names <- intersect(desired_names, seqnames(seq_info))
    len_vec <- seqlengths(seq_info)[existing_names]
  }
  
  return(list(chromlength = len_vec, is_ucsc = is_ucsc))
}

#' Load and Standardize BED File
#' @keywords internal
#' @importFrom data.table fread setnames data.table :=
#' @importFrom stringr str_remove
load_target_bed <- function(target_bed) {
  # Detect header
  first_line <- readLines(target_bed, n = 1)
  has_header <- grepl("^[a-zA-Z]", first_line) && 
    !grepl("^chr[0-9XYM]", first_line, ignore.case = TRUE)
  
  bed <- fread(target_bed, header = has_header)
  if (ncol(bed) < 3) stop("BED file must have at least 3 columns")
  
  setnames(bed, 1:3, c("chromosome", "start", "end"))
  
  # Add IDs and normalize chromosome names
  bed_out <- cbind(
    data.table(target = 1:nrow(bed)),
    bed[, .(chromosome, start, end)]
  )
  bed_out[, chromosome := str_remove(as.character(chromosome), "^chr")]
  
  return(bed_out)
}

#' Create Target Bins
#' @keywords internal
#' @importFrom GenomicRanges makeGRangesFromDataFrame start<- end<- reduce width
create_target_bins <- function(raw_targets) {
  binsize <- 200
  
  # Create GRanges
  gr <- makeGRangesFromDataFrame(raw_targets)
  
  # Extend
  start(gr) <- start(gr) - 50
  end(gr) <- end(gr) + 50
  
  # Ensure minimum width
  widths <- width(gr)
  if (any(widths < 200)) {
    start(gr) <- start(gr) - 50
    end(gr) <- end(gr) + 50
  }
  
  # Merge overlaps
  gr <- reduce(gr)
  
  # Snap to bin size multiples (center aligned)
  w <- width(gr)
  mid <- start(gr) + round(w / 2)
  new_width <- w - w %% binsize
  
  # Split large ranges into smaller bins
  # We construct a data.table to pass to the generic splitter
  dt_ranges <- data.table(
    chromosome = as.character(seqnames(gr)),
    start = mid - new_width / 2,
    end = mid + new_width / 2
  )
  
  final_gr <- split_ranges_into_bins(dt_ranges, binsize)
  return(final_gr)
}

#' Create Background (Antitarget) Bins
#' @keywords internal
#' @importFrom GenomicRanges makeGRangesFromDataFrame reduce gaps start<- end<- width trim
#' @importFrom IRanges trim
#' @importFrom GenomeInfoDb seqlevels<- seqlengths<-
#' @importFrom data.table data.table rbindlist
create_background_bins <- function(raw_targets, chromlength) {
  bg_binsize <- 1e6
  bg_minsize <- 3e5
  
  # 1. Create blocked regions (Targets + Telomeres/Centromeres dummy)
  dummy_targets <- rbind(
    data.table(target = 0, chromosome = names(chromlength), start = 1, end = 10),
    data.table(target = 0, chromosome = names(chromlength), start = chromlength - 10, end = chromlength)
  )
  
  combined <- rbind(
    raw_targets[, .(target, chromosome = as.character(chromosome), start, end)],
    dummy_targets
  )
  
  # 2. Extend and Reduce blocked regions
  gr_base <- makeGRangesFromDataFrame(combined)
  start(gr_base) <- start(gr_base) - 1000
  end(gr_base) <- end(gr_base) + 1000
  gr_base <- reduce(gr_base)
  
  # 3. Invert to get gaps (The background)
  # Must set seqinfo for gaps() to work on the whole genome
  gr_base <- gr_base[as.character(seqnames(gr_base)) %in% names(chromlength)]
  seqlevels(gr_base) <- names(chromlength)
  suppressWarnings(seqlengths(gr_base) <- chromlength)
  gr_base <- IRanges::trim(gr_base)
  
  gr_bg <- gaps(gr_base)
  
  # 4. Filter and Resize
  gr_bg <- gr_bg[width(gr_bg) > bg_minsize]
  
  w <- width(gr_bg)
  mid <- start(gr_bg) + round(w / 2)
  
  # Logic: Only shrink if larger than bin size
  new_width <- w
  mask <- w > bg_binsize
  new_width[mask] <- w[mask] - w[mask] %% bg_binsize
  
  # 5. Split
  dt_ranges <- data.table(
    chromosome = as.character(seqnames(gr_bg)),
    start = mid - round(new_width / 2),
    end = mid + round(new_width / 2)
  )
  
  final_gr <- split_ranges_into_bins(dt_ranges, bg_binsize)
  return(final_gr)
}

#' Generic Range Splitter
#' Splits ranges in a data.table into smaller fixed-size bins
#' @keywords internal
#' @param dt_ranges data.table with chromosome, start, end
#' @param bin_size integer size of bins
split_ranges_into_bins <- function(dt_ranges, bin_size) {
  
  # Calculate how many bins fit in each range
  dt_ranges[, n_bins := ceiling((end - start) / bin_size)]
  
  # Vectorized expansion using lapply
  # (Much faster than row-wise operations if N is large, though similar to original logic)
  bins_list <- lapply(1:nrow(dt_ranges), function(i) {
    s <- dt_ranges$start[i]
    e <- dt_ranges$end[i]
    n <- dt_ranges$n_bins[i]
    chrom <- dt_ranges$chromosome[i]
    
    if (n <= 0) return(NULL)
    
    if (n > 1) {
      starts <- seq(s, e - bin_size, by = bin_size) + 1
      ends <- seq(s + bin_size, e, by = bin_size)
      data.table(chromosome = chrom, start = starts, end = ends)
    } else {
      # Keep as is if it fits or is smaller (though standard logic usually implies fixed size)
      # Original code handled "remainder" or single bins this way
      data.table(chromosome = chrom, start = s + 1, end = e)
    }
  })
  
  bins_dt <- rbindlist(bins_list)
  if (nrow(bins_dt) > 0) {
    bins_dt[, mid := (start + end) / 2]
    bins_dt[, length := end - start]
    return(makeGRangesFromDataFrame(bins_dt))
  } else {
    return(GRanges())
  }
}

#' Create WGS Bins
#' @keywords internal
#' @importFrom data.table data.table rbindlist
#' @importFrom GenomicRanges makeGRangesFromDataFrame
create_wgs_bins <- function(chromlength, bin_size) {
  bins_list <- lapply(names(chromlength), function(chrom) {
    len <- chromlength[chrom]
    starts <- seq(0, len - bin_size * 2, by = bin_size) + 1
    data.table(
      chromosome = chrom,
      start = starts,
      end = starts + bin_size - 1
    )
  })
  bins <- rbindlist(bins_list)
  makeGRangesFromDataFrame(bins)
}

#' Perform BAM Counting
#' @keywords internal
#' @importFrom bamsignals bamCount
perform_bin_counts <- function(bam_file, ranges) {
  flag <- 2816 # QC fail, optical dup, supp align, sec align, unmapped
  
  c_all <- bamCount(bam_file, ranges,
                    paired.end = "midpoint", mapq = 20,
                    filteredFlag = flag, verbose = FALSE
  )
  
  c_short <- bamCount(bam_file, ranges,
                      paired.end = "midpoint", mapq = 20,
                      filteredFlag = flag, tlenFilter = c(0, 150), verbose = FALSE
  )
  
  c_med <- bamCount(bam_file, ranges,
                    paired.end = "midpoint",
                    mapq = 20, filteredFlag = flag, tlenFilter = c(0, 300), verbose = FALSE
  )
  
  list(flag = flag, all = c_all, short = c_short, medium = c_med)
}