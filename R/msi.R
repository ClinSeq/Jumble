#' Extract Flanking Sequences for Indel Variants
#'
#' Retrieves left and right flanking reference sequences from a BSgenome
#' object for each indel variant. Processes by chromosome for memory
#' efficiency with large variant sets (200k-1M+).
#'
#' This function extracts flanks only — the variant sequence itself
#' (deleted or inserted bases) is derived from VCF REF/ALT at the
#' calling layer, not from the reference.
#'
#' @param variants data.table with columns:
#'   \describe{
#'     \item{chrom}{Chromosome name (e.g., "chr1" or "1")}
#'     \item{pos}{1-based genomic position. For deletions: start of the
#'       deleted region (i.e., one past the VCF anchor base). For
#'       insertions: the position after the anchor (insertion point).}
#'     \item{length}{Deletion length in bases. Set to 0 for insertions
#'       (flanks are extracted around the insertion point).}
#'   }
#' @param genome A BSgenome object (e.g., BSgenome.Hsapiens.UCSC.hg19).
#'   The function automatically harmonises chromosome name styles.
#' @param flank_size Integer, number of flanking bases to extract on each
#'   side (default: 25). Values >= 10 are recommended.
#' @param max_del_length Integer, maximum deletion length to extract flanks
#'   for (default: 50). Variants longer than this get NA flanks.
#' @return The input data.table with two new columns: left_flank, right_flank.
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom Biostrings getSeq
#' @importFrom GenomeInfoDb seqnames seqlengths
#' @importFrom data.table set
#' @keywords internal
extract_flanks <- function(variants, genome, flank_size = 25L, max_del_length = 50L) {

  # Get chromosome info from the BSgenome
  genome_chrnames <- GenomeInfoDb::seqnames(genome)
  genome_lengths  <- GenomeInfoDb::seqlengths(genome)
  genome_has_chr  <- any(grepl("^chr", genome_chrnames))

  # Initialise output columns
  data.table::set(variants, j = "left_flank",  value = NA_character_)
  data.table::set(variants, j = "right_flank", value = NA_character_)

  # Eligible: length <= max_del_length (includes length=0 insertions)
  eligible <- which(variants$length <= max_del_length)
  if (length(eligible) == 0) return(variants)

  sub <- variants[eligible]

  # Harmonise chromosome names to match the genome
  sub_chroms <- sub$chrom
  if (genome_has_chr && !any(grepl("^chr", sub_chroms))) {
    sub_chroms <- paste0("chr", sub_chroms)
  } else if (!genome_has_chr && any(grepl("^chr", sub_chroms))) {
    sub_chroms <- sub("^chr", "", sub_chroms)
  }

  # Process by chromosome
  unique_chroms <- unique(sub_chroms)
  for (chr in unique_chroms) {
    if (!(chr %in% genome_chrnames)) {
      warning("Chromosome '", chr, "' not in reference genome, skipping.")
      next
    }
    chr_len <- genome_lengths[chr]
    idx <- which(sub_chroms == chr)

    pos_vec <- sub$pos[idx]
    len_vec <- sub$length[idx]

    # Left flank: ends just before the variant position
    lf_start <- pmax(1L, pos_vec - flank_size)
    lf_end   <- pmax(1L, pos_vec - 1L)

    # Right flank: starts just after the variant region
    # For insertions (length=0), right flank starts at pos itself
    rf_start <- pos_vec + len_vec
    rf_end   <- pmin(chr_len, pos_vec + len_vec + flank_size - 1L)

    # Build GRanges: left flanks + right flanks in one call
    n <- length(idx)
    all_starts <- c(lf_start, rf_start)
    all_ends   <- c(lf_end,   rf_end)

    # Clamp to valid ranges
    all_starts <- pmax(1L, all_starts)
    all_ends   <- pmin(chr_len, all_ends)
    valid <- all_ends >= all_starts
    all_starts[!valid] <- 1L
    all_ends[!valid]   <- 1L

    gr <- GenomicRanges::GRanges(
      seqnames = chr,
      ranges   = IRanges::IRanges(start = all_starts, end = all_ends)
    )

    seqs <- as.character(Biostrings::getSeq(genome, gr))

    lf_seqs <- seqs[1:n]
    rf_seqs <- seqs[(n + 1):(2 * n)]

    # NA for invalid ranges
    lf_seqs[!valid[1:n]]          <- NA_character_
    rf_seqs[!valid[(n+1):(2*n)]]  <- NA_character_

    # Write back
    global_idx <- eligible[idx]
    data.table::set(variants, i = global_idx, j = "left_flank",  value = lf_seqs)
    data.table::set(variants, i = global_idx, j = "right_flank", value = rf_seqs)
  }

  return(variants)
}


#' Classify Indels as MSI-Like (Mono/Di/Trinucleotide Repeats)
#'
#' Classifies each indel as MSI-like based on whether it occurs within a
#' repeat tract, reporting the repeat period. This is the core MSI calling
#' engine: the MSI score for a sample is sum(classify_msi(...) > 0).
#'
#' @param nonanchored_alt Character vector of variant sequences with the
#'   VCF anchor base STRIPPED. For deletions: substring(REF, 2). For
#'   insertions: substring(ALT, 2). The name emphasises that the anchor
#'   base must not be included.
#' @param left_flank  Character vector of upstream flanking sequences.
#' @param right_flank Character vector of downstream flanking sequences.
#' @param max_del_length Integer, maximum variant sequence length to
#'   consider (default: 10). Longer variants return 0.
#'
#' @return Integer vector:
#'   \describe{
#'     \item{0}{Not MSI-like (no qualifying repeat tract found)}
#'     \item{1}{Mononucleotide repeat (e.g., poly-A). Indicates
#'       MLH1/MSH2/PMS2 loss.}
#'     \item{2}{Dinucleotide repeat (e.g., CA-repeat). Indicates
#'       MSH3 loss.}
#'     \item{3}{Trinucleotide repeat (e.g., CAG-repeat).}
#'   }
#'   Shortest qualifying period wins (period 1 checked first).
#'
#' @details
#'   **Algorithm (per variant):**
#'   For each repeat period k in {1, 2, 3}:
#'   1. Extract the candidate motif: the first k characters of the
#'      variant sequence.
#'   2. Check if the variant sequence is > 80% composed of repeats
#'      of that motif.
#'   3. If yes, count how many times the motif repeats consecutively
#'      extending backward into left_flank and forward into right_flank.
#'   4. Compute tract_len = left_repeats + var_repeats + right_repeats
#'      (in motif units, not bases).
#'   5. If tract_len >= threshold for period k, return k.
#'   The first (shortest) qualifying period wins.
#'
#'   **Thresholds (in motif repeat units):**
#'   - Period 1: >= 6 repeats (6bp). Validated: AUC=0.930.
#'   - Period 2: >= 4 repeats (8bp). Provisional — validate on MSH3 data.
#'   - Period 3: >= 3 repeats (9bp). Provisional — validate.
#'
#'   **No external dependencies** — uses base R string operations only.
#'
#' @examples
#' \dontrun{
#'   # After extract_flanks():
#'   # Strip anchor from VCF alleles:
#'   variants[, nonanchored_alt := ifelse(
#'     nchar(ref) > nchar(alt),      # deletion
#'     substring(ref, 2),
#'     substring(alt, 2)             # insertion
#'   )]
#'
#'   result <- classify_msi(variants$nonanchored_alt,
#'                          variants$left_flank,
#'                          variants$right_flank)
#'
#'   msi_count  <- sum(result > 0)   # total MSI-like indels
#'   mono_count <- sum(result == 1)   # mononucleotide (classic MSI)
#'   di_count   <- sum(result == 2)   # dinucleotide (MSH3-like)
#' }
#' @keywords internal
classify_msi <- function(nonanchored_alt, left_flank, right_flank, max_del_length = 10L) {

  # Thresholds: minimum repeat units for each period
  # Period 1: validated (Youden-optimal, AUC=0.930)
  # Period 2-3: provisional — validate on appropriate cohorts
  THRESHOLDS <- c(6L, 4L, 3L)  # index = period

  n <- length(nonanchored_alt)
  stopifnot(length(left_flank) == n, length(right_flank) == n)

  result <- integer(n)  # defaults to 0

  for (i in seq_len(n)) {
    seq_i <- nonanchored_alt[i]
    lf    <- left_flank[i]
    rf    <- right_flank[i]

    # Skip NA or overlength
    if (is.na(seq_i) || is.na(lf) || is.na(rf)) next
    seq_len_i <- nchar(seq_i)
    if (seq_len_i > max_del_length || seq_len_i == 0) next

    # Try each period (shortest first -> first match wins)
    for (period in 1:3) {
      if (seq_len_i < period) next  # sequence shorter than motif

      # Extract candidate motif (first k bases)
      motif <- substr(seq_i, 1, period)

      # Check if sequence is composed of this motif (> 80%)
      # Count complete motif occurrences
      n_motifs_possible <- seq_len_i %/% period
      if (n_motifs_possible == 0) next

      motif_matches <- 0L
      for (m in seq_len(n_motifs_possible)) {
        start_pos <- (m - 1L) * period + 1L
        chunk <- substr(seq_i, start_pos, start_pos + period - 1L)
        if (chunk == motif) motif_matches <- motif_matches + 1L
      }

      # For dominance, count matched bases vs total
      matched_bases <- motif_matches * period
      if (matched_bases / seq_len_i <= 0.8) next  # strict > 0.8

      # Count motif repeats in the variant sequence
      var_repeats <- motif_matches

      # Count consecutive motif repeats extending backward in left_flank
      lf_len <- nchar(lf)
      left_repeats <- 0L
      pos <- lf_len - period + 1L
      while (pos >= 1L) {
        chunk <- substr(lf, pos, pos + period - 1L)
        if (chunk == motif) {
          left_repeats <- left_repeats + 1L
          pos <- pos - period
        } else {
          break
        }
      }

      # Count consecutive motif repeats extending forward in right_flank
      rf_len <- nchar(rf)
      right_repeats <- 0L
      pos <- 1L
      while (pos + period - 1L <= rf_len) {
        chunk <- substr(rf, pos, pos + period - 1L)
        if (chunk == motif) {
          right_repeats <- right_repeats + 1L
          pos <- pos + period
        } else {
          break
        }
      }

      tract_len <- left_repeats + var_repeats + right_repeats

      if (tract_len >= THRESHOLDS[period]) {
        result[i] <- period
        break  # shortest period wins
      }
    }
  }

  return(result)
}


#' Load BSgenome Object by Genome Version
#'
#' Loads the appropriate BSgenome object for the specified genome version.
#'
#' @param genome Character string: "hg19" or "hg38".
#' @return A BSgenome object, or NULL if the package is not available.
#' @keywords internal
load_bsgenome <- function(genome = "hg19") {
  if (genome == "hg19") {
    if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
      return(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
    }
  } else if (genome == "hg38") {
    if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
      return(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
    }
  }
  warning("BSgenome for ", genome, " not available. MSI classification skipped.")
  return(NULL)
}
