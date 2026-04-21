#' Detect Genome Build from BAM Header
#'
#' Inspects the BAM header to guess the genome build (hg19 or hg38).
#'
#' @param bam_file Path to the BAM file.
#' @return A string "hg19", "hg38", or NULL if undetermined.
#' @importFrom Rsamtools scanBamHeader
#' @keywords internal
detect_genome <- function(bam_file) {
  # 1. Header Check ----------------------------------------------------------
  header <- scanBamHeader(bam_file)
  sq <- header[[1]]$targets

  if (length(sq) == 0) {
    return(NULL)
  }

  # Check chr1 length
  # hg19 chr1: 249250621
  # hg38 chr1: 248956422

  len1 <- sq["chr1"]
  if (is.na(len1)) len1 <- sq["1"]

  if (!is.na(len1)) {
    if (len1 == 249250621) {
      return("hg19")
    }
    if (len1 == 248956422) {
      return("hg38")
    }
  }

  return(NULL)
}

#' Compute Backbone Weights
#'
#' Computes continuous weights for bins to be used in regression training (backbone).
#' Prevents densely targeted genes from dominating the training set.
#' - Bins with no gene: weight = 1.0
#' - Bins in genes with <= 25 bins: weight = 1.0
#' - Bins in genes with > 25 bins: weight = 25 / n
#' - Non-autosomal bins (X, Y): weight = 0.0
#'
#' @param targets Data.table containing bin information.
#' @return A numeric vector of weights.
#' @keywords internal
compute_backbone_weights <- function(targets) {
  weights <- rep(1.0, nrow(targets))
  
  # Non-autosomal bins get 0 weight for backbone training
  auto_chroms <- c(as.character(1:22), paste0("chr", 1:22))
  weights[!targets$chromosome %in% auto_chroms] <- 0.0
  
  if ("gene" %in% names(targets)) {
    # Count bins per gene
    gene_counts <- table(targets$gene)
    # Exclude empty gene names
    gene_counts <- gene_counts[names(gene_counts) != ""]
    
    # Identify large genes
    large_genes <- names(gene_counts)[gene_counts > 25]
    
    for (g in large_genes) {
      n_bins <- gene_counts[[g]]
      idx <- which(targets$gene == g)
      weights[idx] <- 25.0 / n_bins
    }
  }
  
  return(weights)
}

#' Define Backbone Bins
#'
#' Identifies bins to be used as backbone for normalization regression training.
#' Replaces the old binary logic with a continuous weight column.
#'
#' @param targets Data.table containing bin information with 'chromosome',
#'   'bin', and optionally 'gene' columns.
#' @return Targets with 'backbone_weight' column.
#' @importFrom data.table copy
#' @keywords internal
define_backbone <- function(targets) {
  # Add continuous backbone weights (replaces binary is_backbone)
  targets[, backbone_weight := compute_backbone_weights(targets)]
  
  return(targets)
}

#' Clean Chromosome Names
#'
#' Standardizes chromosome names by removing 'chr' prefixes.
#'
#' @param x Vector of chromosome names.
#' @return Cleaned chromosome names (1, 2, ..., X, Y).
#' @importFrom stringr str_remove
#' @keywords internal
clean_chrom_names <- function(x) {
  x <- as.character(x)
  stringr::str_remove(x, "^chr")
}

#' Sanitize Legacy Counts Object
#'
#' Removes phantom background bins produced by an early package bug where
#' \code{gaps()} was called without filtering for \code{strand == "*"}. This
#' caused background bins to be generated for the \code{+} and \code{-} strands
#' as well, producing spurious 1 Mb bins that span target regions.
#'
#' Two artefact types are handled:
#' \enumerate{
#'   \item Exact duplicate ranges (same start/end) -- kept once, duplicate removed.
#'   \item Phantom grid -- large bins whose start offset (mod 1 Mb) matches a
#'     confirmed phantom on the same chromosome. A confirmed phantom is any
#'     background bin that overlaps a target bin; valid background bins cannot
#'     overlap target bins by construction (targets are excluded with a 1 kb
#'     margin before gap inversion in \code{create_background_bins}).
#' }
#'
#' @param counts Counts list object as returned by \code{count_reads}.
#' @return Sanitized counts list, or the original list unchanged if no
#'   artefacts are detected.
#' @importFrom GenomicRanges width findOverlaps seqnames start
#' @importFrom S4Vectors queryHits subjectHits
#' @keywords internal
sanitize_legacy_counts <- function(counts) {

  if (is.null(counts$ranges) || length(counts$ranges) == 0) return(counts)

  gr <- counts$ranges

  # ---- Step 1: remove exact duplicate ranges --------------------------------
  dup_mask <- !duplicated(gr)
  if (!all(dup_mask)) {
    message(sprintf(
      "sanitize_legacy_counts: removing %d exact-duplicate bins.",
      sum(!dup_mask)
    ))
    gr <- gr[dup_mask]
    if (!is.null(counts$count))        counts$count        <- counts$count[dup_mask]
    if (!is.null(counts$count_short))  counts$count_short  <- counts$count_short[dup_mask]
    if (!is.null(counts$count_medium)) counts$count_medium <- counts$count_medium[dup_mask]
  }

  # ---- Step 2: detect and remove phantom grid -------------------------------
  bg_mask <- GenomicRanges::width(gr) >= 100000
  bg_gr   <- gr[bg_mask]
  tgt_gr  <- gr[!bg_mask]

  # Early exit: no background-background overlaps → already clean.
  if (length(bg_gr) == 0 || length(tgt_gr) == 0) {
    counts$ranges <- gr
    return(counts)
  }
  bg_bg_olaps <- GenomicRanges::findOverlaps(bg_gr, bg_gr)
  bg_bg_olaps <- bg_bg_olaps[S4Vectors::queryHits(bg_bg_olaps) != S4Vectors::subjectHits(bg_bg_olaps)]
  if (length(bg_bg_olaps) == 0) {
    counts$ranges <- gr
    return(counts)
  }

  # Confirmed phantoms: background bins overlapping any target bin.
  # Valid background bins cannot overlap targets by construction.
  olaps <- GenomicRanges::findOverlaps(bg_gr, tgt_gr)
  confirmed_bg_local <- unique(S4Vectors::queryHits(olaps))

  if (length(confirmed_bg_local) == 0) {
    counts$ranges <- gr
    return(counts)
  }

  # Per chromosome: phantom grid fingerprint = start %% 1e6 of confirmed phantoms.
  # Remove all background bins on that chromosome sharing that fingerprint.
  bg_chroms <- as.character(GenomicRanges::seqnames(bg_gr))
  bg_starts <- GenomicRanges::start(bg_gr)

  phantom_bg_local <- integer(0)
  for (ch in unique(bg_chroms[confirmed_bg_local])) {
    ch_confirmed  <- confirmed_bg_local[bg_chroms[confirmed_bg_local] == ch]
    offsets       <- unique(bg_starts[ch_confirmed] %% 1e6)
    ch_all        <- which(bg_chroms == ch)
    ch_phantoms   <- ch_all[bg_starts[ch_all] %% 1e6 %in% offsets]
    phantom_bg_local <- c(phantom_bg_local, ch_phantoms)
  }
  phantom_bg_local <- unique(phantom_bg_local)

  if (length(phantom_bg_local) > 0) {
    phantom_global <- which(bg_mask)[phantom_bg_local]
    message(sprintf(
      "sanitize_legacy_counts: removing %d phantom background bins (strand-artefact grid).",
      length(phantom_global)
    ))
    keep <- seq_along(gr)[-phantom_global]
    gr   <- gr[keep]
    if (!is.null(counts$count))        counts$count        <- counts$count[keep]
    if (!is.null(counts$count_short))  counts$count_short  <- counts$count_short[keep]
    if (!is.null(counts$count_medium)) counts$count_medium <- counts$count_medium[keep]
  }

  counts$ranges <- gr
  return(counts)
}
