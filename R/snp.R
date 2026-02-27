#' Match and Subset VCF by Sample Name
#'
#' Determines the correct sample column in a multi-sample VCF based on an exact,
#' fuzzy, or reverse-fuzzy match of the provided sample name.
#'
#' @param vcf A VCF object.
#' @param sample_name Optional string representing the sample name to find.
#' @return A VCF object subsetted to a single sample column.
#' @keywords internal
match_vcf_sample <- function(vcf, sample_name = NULL) {
  vcf_names <- colnames(vcf)
  
  if (length(vcf_names) <= 1) {
    return(vcf)
  }
  
  ix <- 1 # Default to the first column
  
  if (!is.null(sample_name)) {
    exact_match <- which(vcf_names == sample_name)
    fuzzy_match <- which(grepl(sample_name, vcf_names, fixed = TRUE))
    fuzzy_match_rev <- which(sapply(vcf_names, function(vn) grepl(vn, sample_name, fixed = TRUE)))
    
    if (length(exact_match) == 1) {
      ix <- exact_match
      message("Matched VCF sample column '", vcf_names[ix], "' to sample '", sample_name, "'")
    } else if (length(fuzzy_match) >= 1) {
      ix <- fuzzy_match[1]
      message("Fuzzy matched VCF sample column '", vcf_names[ix], "' to sample '", sample_name, "'")
    } else if (length(fuzzy_match_rev) >= 1) {
      ix <- fuzzy_match_rev[1]
      message("Reverse fuzzy matched VCF sample column '", vcf_names[ix], "' to sample '", sample_name, "'")
    } else {
      warning("Could not match sample name '", sample_name, "' to any VCF column. Using first column: ", vcf_names[1])
    }
  } else {
    warning("Multiple samples in VCF but no sample_name provided. Using first column: ", vcf_names[1])
  }
  
  return(vcf[, ix])
}

#' Extract and Filter SNP Table from VCF
#'
#' Filters a single-sample VCF for biallelic SNPs and constructs a data.table 
#' containing allele depths, calculating allele ratios.
#'
#' @param vcf A single-sample VCF object.
#' @return A processed data.table containing SNP information.
#' @importFrom VariantAnnotation alt ref geno
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom S4Vectors elementNROWS
#' @importFrom GenomeInfoDb seqnames
#' @importFrom stringr str_remove str_detect str_length
#' @importFrom data.table data.table as.data.table :=
#' @keywords internal
extract_snp_table <- function(vcf) {
  # 1. Filter for biallelic sites with exactly 2 AD values (REF + ALT)
  n_alt <- elementNROWS(alt(vcf))
  ad_len <- elementNROWS(geno(vcf)$AD)
  vcf <- vcf[n_alt == 1 & ad_len == 2]
  
  sample_colname <- colnames(vcf)
  g <- geno(vcf)
  rr <- rowRanges(vcf)
  
  # 2. Extract genomic coordinates
  chroms <- as.character(seqnames(rr))
  chroms <- str_remove(chroms, "^chr")
  
  snp_table <- data.table(
    sample = sample_colname,
    id = names(rr),
    chromosome = chroms,
    start = start(rr),
    end = end(rr)
  )
  
  # 3. Extract alleles
  snp_table$ref_allele <- as.character(ref(vcf))
  snp_table$alt_allele <- as.character(unlist(alt(vcf)))
  snp_table$n_alt_alleles <- 1 
  
  # 4. Standardize SNP types
  snp_table[!str_detect(ref_allele, "^[ACGT]$"), ref_allele := "other"]
  snp_table[!str_detect(alt_allele, "^[ACGT]$"), alt_allele := "other"]
  snp_table[, type := paste0(ref_allele, ">", alt_allele)]
  snp_table[str_length(ref_allele) != 1 | str_length(alt_allele) != 1, type := "other"]
  
  # 5. Get read counts and calculate ratios
  ad_list <- g$AD[, 1] # Extract vector of lists for the single sample
  snp_table$RD <- sapply(ad_list, "[[", 1)
  snp_table$AD <- sapply(ad_list, "[[", 2)
  
  snp_table[, DP := AD + RD]
  snp_table[, logDP := log2(DP)]
  
  raw_allele_ratio <- round(snp_table$AD / snp_table$DP, 4)
  raw_allele_ratio[is.nan(raw_allele_ratio)] <- 0
  snp_table$allele_ratio <- raw_allele_ratio
  
  snp_table[, snp := paste(id, ref_allele, alt_allele)]
  
  return(snp_table)
}

#' Map SNPs to Target Bins
#'
#' Finds overlaps between extracted SNPs and valid genomic target bins, 
#' dropping off-target SNPs.
#'
#' @param snp_table A data.table of processed SNPs.
#' @param targets A data.table of target regions containing 'chromosome', 'start', 'end', 'bin', and 'is_target'.
#' @return A data.table of SNPs filtered and assigned to target bins.
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom data.table :=
#' @keywords internal
map_snps_to_bins <- function(snp_table, targets) {
  snp_gr <- makeGRangesFromDataFrame(
    snp_table,
    seqnames.field = "chromosome",
    start.field = "start",
    end.field = "end"
  )
  
  targets_gr <- makeGRangesFromDataFrame(
    targets,
    seqnames.field = "chromosome",
    start.field = "start",
    end.field = "end"
  )
  
  overlap <- findOverlaps(snp_gr, targets_gr)
  matched_indices <- subjectHits(overlap)
  snp_indices <- queryHits(overlap)
  
  # Map mapped bin IDs to SNPs
  snp_table[, bin := as.integer(NA)]
  snp_table[snp_indices, bin := targets$bin[matched_indices]]
  
  # Filter for SNPs falling within valid target bins
  valid_bins <- targets[is_target == TRUE]$bin
  snp_table <- snp_table[bin %in% valid_bins]
  
  return(snp_table)
}

#' Process SNPs from VCF
#'
#' Main wrapper: Loads and processes SNPs from a VCF file, calculates allele 
#' ratios, and maps the results to genomic target bins.
#'
#' @param vcf_file Path to the VCF file.
#' @param targets A data.table of target regions.
#' @param sample_name Optional sample name to match in VCF columns.
#' @return A data.table containing processed SNP information, or NULL if processing fails.
#' @importFrom VariantAnnotation readVcf
#' @importFrom stringr str_detect
#' @keywords internal
process_snps <- function(vcf_file, targets, sample_name = NULL) {
  if (is.null(vcf_file) || !file.exists(vcf_file)) {
    return(NULL)
  }
  
  # 1. Validation
  if (!str_detect(vcf_file, "\\.[vV][cC][fF]$") && !str_detect(vcf_file, "\\.[vV][cC][fF]\\.[gG][zZ]$")) {
    warning("SNP vcf file appears incorrect (extension check)")
    return(NULL)
  }
  
  message("Processing SNPs from VCF: ", vcf_file)
  
  # 2. Execution block with error handling
  tryCatch({
    # Read VCF
    vcf <- readVcf(vcf_file)
    
    # Subset to correct sample
    vcf <- match_vcf_sample(vcf, sample_name)
    
    # Extract, filter, and format SNPs
    snp_table <- extract_snp_table(vcf)
    
    # Map to genomic target bins
    final_snp_table <- map_snps_to_bins(snp_table, targets)
    
    return(final_snp_table)
  }, error = function(e) {
    warning("Failed to process SNP VCF: ", e$message)
    return(NULL)
  })
}