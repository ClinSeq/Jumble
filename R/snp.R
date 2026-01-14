#' Process SNPs from VCF
#'
#' Loads and processes SNPs from a VCF file, calculating allele ratios and mapping to bins.
#'
#' @param vcf_file Path to the VCF file.
#' @param sample_name Optional sample name to match in VCF columns.
#' @return A data.table containing processed SNP information, or NULL if processing fails.
#' @importFrom VariantAnnotation readVcf ref alt geno
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom data.table as.data.table := data.table
#' @importFrom stringr str_detect str_remove str_length
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom S4Vectors queryHits subjectHits elementNROWS
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges IRanges
#' @importFrom methods as
#' @export
process_snps <- function(vcf_file, targets, sample_name = NULL) {
  if (is.null(vcf_file) || !file.exists(vcf_file)) {
    return(NULL)
  }

  message("Processing SNPs from VCF: ", vcf_file)

  # 1. Validation ------------------------------------------------------------
  # Check file extension
  if (!str_detect(vcf_file, ".[vV][cC][fF]$") &
    !str_detect(vcf_file, ".[vV][cC][fF].[gG][zZ]$")) {
    warning("SNP vcf file appears incorrect (extension check)")
    return(NULL)
  }

  # 2. Read VCF --------------------------------------------------------------
  tryCatch(
    {
      vcf <- readVcf(vcf_file)

      # If there is more than one sample in the VCF, which one to use?
      ix <- 1 # base assumption: first is this sample
      vcf_names <- colnames(vcf)
      
      if (length(vcf_names) > 1) {
        if (!is.null(sample_name)) {
            # Try to find a match for the provided sample_name
            # Exact match first
            exact_match <- which(vcf_names == sample_name)
            if (length(exact_match) == 1) {
                ix <- exact_match
                message("Matched VCF sample column '", vcf_names[ix], "' to sample '", sample_name, "'")
            } else {
                # Fuzzy match (substring)
                # Ensure we don't match multiple (pick first if so)
                fuzzy_match <- which(grepl(sample_name, vcf_names, fixed = TRUE))
                if (length(fuzzy_match) >= 1) {
                    ix <- fuzzy_match[1]
                    message("Fuzzy matched VCF sample column '", vcf_names[ix], "' to sample '", sample_name, "'")
                } else {
                    # Try reverse fuzzy (vcf name inside sample name - rare but possible)
                    fuzzy_match_rev <- which(sapply(vcf_names, function(vn) grepl(vn, sample_name, fixed = TRUE)))
                    if (length(fuzzy_match_rev) >= 1) {
                        ix <- fuzzy_match_rev[1]
                        message("Reverse fuzzy matched VCF sample column '", vcf_names[ix], "' to sample '", sample_name, "'")
                    } else {
                        warning("Could not match sample name '", sample_name, "' to any VCF column (", paste(vcf_names, collapse=","), "). Using first column: ", vcf_names[1])
                    }
                }
            }
        } else {
             warning("Multiple samples in VCF but no sample_name provided. Using first column: ", vcf_names[1])
        }
        vcf <- vcf[, ix]
      }

      # 3. Filter Variants -------------------------------------------------------
      # Filter SNPs to ensure clean allele ratio calculation
      # We require:
      # 1. Exactly one ALT allele (n_alt == 1), meaning the site is biallelic
      #    (REF + 1 ALT).
      # 2. The AD (Allele Depth) field must have exactly 2 values (depth for
      #    REF and depth for ALT).
      # This simplifies downstream processing to a straightforward REF vs ALT
      # comparison.

      # Check number of alt alleles per site
      n_alt <- elementNROWS(alt(vcf))

      # Check AD field length
      ad_len <- elementNROWS(geno(vcf)$AD)

      # Apply filter
      vcf <- vcf[n_alt == 1 & ad_len == 2]

      name <- colnames(vcf)
      g <- geno(vcf)

      # 4. Create SNP Table ------------------------------------------------------
      rr <- rowRanges(vcf)

      # Extract chromosome and handle chr prefix
      # Note: seqnames returns factor-Rle, convert to character
      chroms <- as.character(seqnames(rr))
      chroms <- str_remove(chroms, "^chr")

      snp_table <- data.table(
        chromosome = chroms,
        start = start(rr),
        end = end(rr)
      )

      snp_table <- cbind(
        data.table(
          sample = name,
          id = names(rr) # row names of VCF are IDs
        ),
        snp_table
      )

      # Get alleles
      snp_table$ref_allele <- as.character(ref(vcf))

      # Get alt allele (first one, since we filtered for length 1)
      # alt(vcf) is DNAStringSetList. unlist or sapply.
      # Original: alt_allele <- as.data.table(alt(vcf))[, .(values =
      # list(value)), by = group][,values]
      # Simpler:
      alts <- alt(vcf)
      snp_table$alt_allele <- as.character(unlist(alts))

      snp_table$n_alt_alleles <- 1 # We filtered for this

      # Get SNP type
      # Original logic:
      # from <- rep('C/G',length(vcf)); from[snp_table$ref_allele %in%
      # c('A','T')] <- 'A/T'
      # ...
      # snp_table[!str_detect(ref_allele,'^[ACGT]$'),ref_allele:='other']
      # ...

      # We can replicate exactly:
      snp_table[!str_detect(ref_allele, "^[ACGT]$"), ref_allele := "other"]
      snp_table[!str_detect(alt_allele, "^[ACGT]$"), alt_allele := "other"]
      snp_table[, type := paste0(ref_allele, ">", alt_allele)]
      snp_table[
        str_length(ref_allele) != 1 | str_length(alt_allele) != 1,
        type := "other"
      ]

      # Get read counts
      # g$AD is a list of integer vectors (per site)
      # We need to extract 1st (Ref) and 2nd (Alt) values
      # Since we have 1 sample, g$AD is a list (length = n_sites)
      # But wait, geno(vcf)$AD returns a matrix if multiple samples, or list?
      # It returns a matrix of lists. Since we subsetted to 1 sample:
      ad_list <- g$AD[, 1] # Vector of lists

      snp_table$RD <- sapply(ad_list, "[[", 1)
      snp_table$AD <- sapply(ad_list, "[[", 2)

      snp_table[, DP := AD + RD]
      snp_table[, logDP := log2(DP)]

      # Compute raw allele ratio
      raw_allele_ratio <- round(snp_table$AD / snp_table$DP, 4)
      raw_allele_ratio[is.nan(raw_allele_ratio)] <- 0
      snp_table$allele_ratio <- raw_allele_ratio

      # SNP id with details
      snp_table[, snp := paste(id, ref_allele, alt_allele)]

      # 5. Map to Bins -----------------------------------------------------------
      # Overlap with bins, drop SNPs not on targets
      # targets needs to be GRanges for overlap
      # But targets is data.table.
      # We assume targets has chromosome, start, end.

      # Ensure targets chromosome format matches snp_table (no chr prefix)
      # targets usually has 1, 2, ... X, Y. snp_table has same (we removed chr).

      snp_gr <- makeGRangesFromDataFrame(snp_table,
        seqnames.field = "chromosome",
        start.field = "start",
        end.field = "end"
      )

      targets_gr <- makeGRangesFromDataFrame(targets,
        seqnames.field = "chromosome",
        start.field = "start",
        end.field = "end"
      )

      overlap <- findOverlaps(snp_gr, targets_gr)

      # Map SNPs to bins
      # We use findOverlaps to identify which target bin each SNP falls into.
      # 'subjectHits(overlap)' gives the index of the target bin in the
      # 'targets' table.
      # We use this index to retrieve the actual 'bin' ID from targets$bin.

      matched_indices <- subjectHits(overlap)
      snp_indices <- queryHits(overlap)

      # Initialize bin column
      snp_table[, bin := as.integer(NA)]
      snp_table[snp_indices, bin := targets$bin[matched_indices]]

      # Filter for SNPs that map to valid target bins
      # This removes SNPs that fall into off-target regions or unmapped bins.
      valid_bins <- targets[is_target == TRUE]$bin
      snp_table <- snp_table[bin %in% valid_bins]

      return(snp_table)
    },
    error = function(e) {
      warning("Failed to process SNP VCF: ", e$message)
      return(NULL)
    }
  )
}
