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
