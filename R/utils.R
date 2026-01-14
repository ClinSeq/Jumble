#' Detect Genome Build from BAM Header
#'
#' Inspects the BAM header to guess the genome build (hg19 or hg38).
#'
#' @param bam_file Path to the BAM file.
#' @return A string "hg19", "hg38", or NULL if undetermined.
#' @importFrom Rsamtools scanBamHeader
#' @export
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

#' Define Backbone Bins
#'
#' Identifies bins to be used as backbone for normalization.
#' Uses the variance across the reference set to select stable bins.
#'
#' @param targets Data.table containing bin information.
#' @param reference Reference object containing 'mat' and 'bins'.
#' @return Targets with 'is_backbone' column.
#' @importFrom data.table copy
#' @importFrom stats sd quantile
#' @export
define_backbone <- function(targets, reference = NULL) {
  # 1. Default Backbone ------------------------------------------------------
  # Default: autosomes are backbone candidates
  targets[, is_backbone := chromosome %in% c(as.character(1:22), 1:22)]

  # 2. Refine with Variance --------------------------------------------------
  # If reference is provided, use variance to refine backbone
  if (!is.null(reference) && !is.null(reference$mat) &&
    !is.null(reference$bins)) {
    # Calculate SD for each bin in reference
    mat <- reference$mat

    # Ensure mat rows match reference$bins
    if (nrow(mat) == length(reference$bins)) {
      sds <- apply(mat, 1, sd, na.rm = TRUE)

      # Select bins with lowest 50% variance
      # Note: This assumes 'mat' rows correspond to 'reference$bins' in order.
      # And 'reference$bins' are the bin IDs.

      cutoff <- quantile(sds, 0.5, na.rm = TRUE)
      backbone_bins <- reference$bins[sds < cutoff]

      # Update is_backbone
      # Only consider bins that are ALREADY autosomes (from above) AND in low
      # variance set
      targets[, is_backbone := is_backbone & (bin %in% backbone_bins)]
    }
  }

  return(targets)
}

#' Clean Chromosome Names
#'
#' Standardizes chromosome names by removing 'chr' prefixes.
#'
#' @param x Vector of chromosome names.
#' @return Cleaned chromosome names (1, 2, ..., X, Y).
#' @importFrom stringr str_remove
#' @export
clean_chrom_names <- function(x) {
  x <- as.character(x)
  stringr::str_remove(x, "^chr")
}
