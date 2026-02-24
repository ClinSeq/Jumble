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
#' Identifies bins to be used as backbone for normalization regression training.
#' All autosomal bins are included by default. Densely-targeted genes (BRCA2,
#' PTEN) are capped to a maximum number of backbone bins to prevent their
#' gene bodies from dominating the regression fit when they harbour CNAs.
#'
#' @param targets Data.table containing bin information with 'chromosome',
#'   'bin', and optionally 'gene' columns.
#' @return Targets with 'is_backbone' column.
#' @importFrom data.table copy
#' @export
define_backbone <- function(targets) {
  # 1. Default Backbone ------------------------------------------------------
  # All autosomes are backbone candidates
  targets[, is_backbone := chromosome %in% c(as.character(1:22), 1:22)]

  # 2. Cap Dense Genes -------------------------------------------------------
  # Restrict the contribution of densely-targeted genes to prevent them from
  # dominating the regression training set. Hard-coded for BRCA2 and PTEN as
  # per the published Jumble method.
  max_in_b <- 10
  genes_to_cap <- c("BRCA2", "PTEN")

  if ("gene" %in% names(targets) && any(targets$gene != "")) {
    set.seed(25)
    for (g in genes_to_cap) {
      genebins <- targets[gene == g]$bin
      n <- length(genebins)
      if (n > max_in_b) {
        targets[bin %in% genebins, is_backbone := FALSE]
        genebins <- sample(genebins, max_in_b, replace = FALSE)
        targets[bin %in% genebins, is_backbone := TRUE]
      }
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
