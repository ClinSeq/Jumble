#' Extract Base Somatic Variants and Allele Metrics
#'
#' Extracts chromosome, position, reference, and alternate alleles from a VCF, 
#' and calculates Allele Frequency (AF), Allele Observation (AO), and Depth (DP).
#'
#' @param vcf A VCF object.
#' @param tumor_idx Integer indicating the column index of the tumor sample.
#' @return A data.table containing unannotated somatic variants.
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomicRanges seqnames start end
#' @importFrom VariantAnnotation geno
#' @importFrom data.table data.table set :=
#' @keywords internal
extract_base_variants <- function(vcf, tumor_idx) {
  r <- SummarizedExperiment::rowRanges(vcf)
  g <- VariantAnnotation::geno(vcf)
  samples <- colnames(vcf)
  
  # Initialize base data.table
  somatic <- data.table::data.table(
    chromosome = gsub("^chr", "", as.character(GenomicRanges::seqnames(r))),
    start = GenomicRanges::start(r),
    end = GenomicRanges::end(r),
    REF = as.character(r$REF),
    ALT = as.character(unlist(lapply(r$ALT, function(x) as.character(x)[1])))
  )
  
  # Extract AF
  if (!is.null(g$VAF)) {
    af_dat <- g$VAF[, tumor_idx]
    somatic$AF <- if (is.list(af_dat)) sapply(af_dat, function(x) if (length(x) > 0) as.numeric(x[1]) else 0) else as.numeric(af_dat)
  } else if (!is.null(g$AF)) {
    af_dat <- g$AF[, tumor_idx]
    somatic$AF <- if (is.list(af_dat)) sapply(af_dat, function(x) if (length(x) > 0) x[1] else 0) else as.numeric(af_dat)
  } else if (!is.null(g$AD)) {
    ad_dat <- g$AD[, tumor_idx]
    somatic$AF <- sapply(ad_dat, function(x) if (length(x) >= 2 && (x[1] + x[2]) > 0) x[2] / (x[1] + x[2]) else 0)
  } else {
    warning("No AF, VAF, or AD field found in VCF. Setting AF to 0.")
    somatic$AF <- 0
  }
  
  # Extract AO and DP
  if (!is.null(g$AD)) {
    ad_dat <- g$AD[, tumor_idx]
    somatic$AO <- sapply(ad_dat, function(x) if (length(x) >= 2 && !is.na(x[2])) x[2] else 0)
    somatic$DP <- sapply(ad_dat, function(x) if (length(x) > 0) sum(x, na.rm = TRUE) else 0)
  } else if (!is.null(g$DP)) {
    dp_val <- if (length(dim(g$DP)) == 3 && tumor_idx <= dim(g$DP)[2]) {
      dp_slice <- g$DP[, tumor_idx, ]
      if (is.matrix(dp_slice)) rowSums(dp_slice, na.rm = TRUE) else sum(dp_slice, na.rm = TRUE)
    } else {
      dp_raw <- as.numeric(g$DP)
      if (length(dp_raw) == nrow(somatic) * length(samples)) {
        matrix(dp_raw, nrow = nrow(somatic), ncol = length(samples))[, tumor_idx]
      } else if (length(dp_raw) == nrow(somatic)) {
        dp_raw
      } else {
        rep(NA, nrow(somatic))
      }
    }
    somatic$DP <- dp_val
    somatic$AO <- as.integer(round(somatic$AF * somatic$DP))
  } else {
    somatic$AO <- 5
    somatic$DP <- 10
  }
  
  # Clean up NAs
  somatic[is.na(AO), AO := 0]
  somatic[is.na(DP), DP := 0]
  
  return(somatic)
}

#' Parse VEP CSQ Annotations
#'
#' Parses the CSQ string from the VCF INFO column for a subset of variants.
#'
#' @param vcf A VCF object.
#' @param keep_mask Logical vector indicating which rows to process.
#' @return A data.table containing the parsed CSQ columns.
#' @importFrom VariantAnnotation info header
#' @importFrom data.table data.table set
#' @keywords internal
parse_vep_csq <- function(vcf, keep_mask) {
  vcf_header_info <- VariantAnnotation::info(VariantAnnotation::header(vcf))
  csq_desc <- if ("CSQ" %in% rownames(vcf_header_info)) vcf_header_info$Description[rownames(vcf_header_info) == "CSQ"] else character(0)
  
  n_rows <- sum(keep_mask)
  cols <- c("SYMBOL", "Consequence", "Protein_position", "CANONICAL", "IMPACT", "CLIN_SIG", "Amino_acids", "MAX_AF")
  res <- data.table::data.table(matrix(NA_character_, nrow = n_rows, ncol = length(cols)))
  names(res) <- cols
  
  if (length(csq_desc) == 0 || n_rows == 0) {
    # Fallback: try INFO/GENE and INFO/EFFECT when VEP CSQ is absent
    if (n_rows > 0) {
      info_data <- VariantAnnotation::info(vcf)
      if ("GENE" %in% names(info_data)) {
        res$SYMBOL <- as.character(info_data$GENE[keep_mask])
      }
      if ("EFFECT" %in% names(info_data)) {
        effects <- as.character(info_data$EFFECT[keep_mask])
        # Parse "Missense_p.G245S" -> Consequence + Amino_acids + Protein_position
        parts <- strsplit(effects, "_p\\.")
        res$Consequence <- vapply(parts, function(x) {
          ct <- x[1]
          ct <- gsub("Missense", "missense_variant", ct)
          ct <- gsub("Frameshift", "frameshift_variant", ct)
          ct <- gsub("Nonsense", "stop_gained", ct)
          ct <- gsub("Splice", "splice_region_variant", ct)
          ct <- gsub("Synonymous", "synonymous_variant", ct)
          ct
        }, character(1))
        res$IMPACT <- ifelse(grepl("frameshift|stop_gained", res$Consequence), "HIGH",
                      ifelse(grepl("missense", res$Consequence), "MODERATE",
                      ifelse(grepl("splice", res$Consequence), "LOW", "MODIFIER")))
        res$CANONICAL <- "YES"
        # Extract amino acids and protein position
        pchanges <- vapply(parts, function(x) if (length(x) >= 2) x[2] else "", character(1))
        aa_matches <- regmatches(pchanges, regexec("^([A-Z*])(\\d+)(.*)", pchanges))
        res$Protein_position <- vapply(aa_matches, function(m) if (length(m) == 4) m[3] else "", character(1))
        res$Amino_acids <- vapply(aa_matches, function(m) if (length(m) == 4) paste0(m[2], "/", m[4]) else "", character(1))
      }
      if ("CLINVAR" %in% names(info_data)) {
        res$CLIN_SIG <- as.character(info_data$CLINVAR[keep_mask])
      }
    }
    return(res)
  }
  
  csq_format <- sub(".*Format: ", "", csq_desc)
  csq_fields <- strsplit(csq_format, "\\|")[[1]]
  
  csq_data <- VariantAnnotation::info(vcf)$CSQ
  if (length(csq_data) != length(keep_mask)) {
    warning("Length mismatch between VCF info and Rows. CSQ parsing may be misaligned.")
    return(res)
  }
  
  csq_data <- csq_data[keep_mask]
  
  idx <- sapply(cols, function(x) which(csq_fields == x)[1]) # Named vector of indices
  
  parsed_csq <- lapply(csq_data, function(csq_strings) {
    if (length(csq_strings) == 0) return(as.list(setNames(rep("", length(cols)), cols)))
    
    parsed_list <- strsplit(csq_strings, "\\|")
    is_canonical <- sapply(parsed_list, function(x) if (!is.na(idx["CANONICAL"]) && length(x) >= idx["CANONICAL"]) x[idx["CANONICAL"]] == "YES" else FALSE)
    
    sel_idx <- if (any(is_canonical)) which(is_canonical)[1] else 1
    sel <- parsed_list[[sel_idx]]
    
    lapply(idx, function(i) if (!is.na(i) && length(sel) >= i) sel[i] else "")
  })
  
  for (col in cols) {
    data.table::set(res, j = col, value = sapply(parsed_csq, "[[", col))
  }
  
  return(res)
}

#' Annotate Hotspots and Variant Effects
#'
#' Maps somatic mutations to known hotspot tables and classifies their impact.
#'
#' @param somatic A filtered, CSQ-annotated somatic data.table.
#' @param genome Genome build version ("hg19" or "hg38").
#' @return The mutated data.table with `hotkey`, `is_hotspot`, and `effect` columns.
#' @importFrom stringr str_match str_remove
#' @importFrom data.table :=
#' @keywords internal
annotate_hotspots <- function(somatic, genome) {
  if (nrow(somatic) == 0) return(somatic)
  
  hotspots_loaded <- FALSE
  tryCatch({
    data("hotspots_snvs", package = "Jumble", envir = environment())
    data("hotspots_inframes", package = "Jumble", envir = environment())
    
    if (genome == "hg38") {
      data("hotspots_splice_hg38", package = "Jumble", envir = environment())
      if (exists("hotspots_splice_hg38")) hotspots_splice <- hotspots_splice_hg38
    } else {
      data("hotspots_splice_hg19", package = "Jumble", envir = environment())
      if (exists("hotspots_splice_hg19")) hotspots_splice <- hotspots_splice_hg19
    }
    if (exists("hotspots_snvs") && exists("hotspots_inframes")) hotspots_loaded <- TRUE
  }, error = function(e) {
    warning("Hotspot data not found or failed to load. Skipping hotspot annotation.")
  })
  
  if (hotspots_loaded) {
    somatic[, hotkey := paste(SYMBOL, stringr::str_match(Protein_position, "^([0-9]*[-]*[0-9]*)/")[, 2])]
    somatic[SYMBOL == "" | is.na(SYMBOL), hotkey := NA_character_]
    somatic[grepl(" NA$", hotkey), hotkey := NA_character_]
    
    somatic[grepl("frameshift", Consequence), hotkey := paste(hotkey, "fs")]
    somatic[grepl("inframe", Consequence), hotkey_start := stringr::str_remove(hotkey, "-[0-9]+")]
    somatic[grepl("inframe", Consequence), hotkey_end := stringr::str_remove(hotkey, "[0-9]+-")]
    
    somatic[, is_hotspot := FALSE]
    if (exists("hotspots_snvs")) somatic[hotkey %in% hotspots_snvs, is_hotspot := TRUE]
    if (exists("hotspots_inframes")) {
      somatic[hotkey_start %in% hotspots_inframes | hotkey_end %in% hotspots_inframes, is_hotspot := TRUE]
    }
    if (exists("hotspots_splice")) {
      somatic[paste(chromosome, start, sep = ":") %in% hotspots_splice, is_hotspot := TRUE]
    }
    
    # TERT promoter
    tert_coords <- if (genome == "hg19") c("5:1295228", "5:1295250") else c("5:1295113", "5:1295135")
    somatic[paste(chromosome, start, sep = ":") %in% tert_coords, is_hotspot := TRUE]
    
    somatic[, effect := NA_character_]
    ix <- (somatic$CANONICAL == "YES" | !is.na(somatic$SYMBOL))
    
    if (!all(is.na(somatic$IMPACT))) {
      somatic[ix & IMPACT == "MODERATE", effect := "uncertain"]
      somatic[ix & IMPACT == "HIGH", effect := "high-impact"]
    }
    if (!all(is.na(somatic$CLIN_SIG))) {
      somatic[ix & grepl("pathogenic", CLIN_SIG), effect := "high-impact"]
    }
    somatic[ix & is_hotspot == TRUE, effect := "hotspot"]
  }
  return(somatic)
}

#' Classify Somatic Indels for MSI
#'
#' Identifies indels in the somatic table and classifies them using the
#' MSI repeat-tract classifier. Adds an MSI column (0/1/2/3) for indels,
#' NA for SNVs.
#'
#' @param somatic A processed somatic data.table (must have REF, ALT columns).
#' @param genome Genome build version ("hg19" or "hg38").
#' @return The somatic data.table with an appended `MSI` column.
#' @importFrom data.table := set
#' @keywords internal
classify_somatic_msi <- function(somatic, genome) {
  somatic[, MSI := NA_integer_]

  if (nrow(somatic) == 0) return(somatic)

  # Identify indels (REF and ALT differ in length)
  is_indel <- nchar(somatic$REF) != nchar(somatic$ALT)
  if (sum(is_indel) == 0) return(somatic)

  # Load BSgenome
  bsg <- load_bsgenome(genome)
  if (is.null(bsg)) return(somatic)

  # Prepare inputs for extract_flanks
  indel_idx <- which(is_indel)
  indel_sub <- somatic[indel_idx]

  is_del <- nchar(indel_sub$REF) > nchar(indel_sub$ALT)

  flank_input <- data.table::data.table(
    chrom  = indel_sub$chromosome,
    pos    = indel_sub$start + 1L,          # skip VCF anchor base
    length = ifelse(is_del, nchar(indel_sub$REF) - 1L, 0L)
  )

  # Extract flanking sequences
  flank_input <- extract_flanks(flank_input, bsg)

  # Derive nonanchored_alt: strip VCF anchor base
  nonanchored_alt <- ifelse(
    is_del,
    substring(indel_sub$REF, 2),
    substring(indel_sub$ALT, 2)
  )

  # Classify
  msi_result <- classify_msi(
    nonanchored_alt,
    flank_input$left_flank,
    flank_input$right_flank
  )

  data.table::set(somatic, i = indel_idx, j = "MSI", value = msi_result)

  n_msi <- sum(msi_result > 0)
  message("MSI classification: ", sum(is_indel), " indels, ",
          n_msi, " MSI-like (", sum(msi_result == 1), " mono, ",
          sum(msi_result == 2), " di, ", sum(msi_result == 3), " tri)")

  return(somatic)
}


#' Map Somatic Variants to Reference Bins
#'
#' @param somatic A processed data.table of somatic variants.
#' @param reference The Jumble reference object containing `$ranges`.
#' @return The data.table with an appended `bin` column.
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors subjectHits queryHits
#' @importFrom data.table :=
#' @keywords internal
map_variants_to_bins <- function(somatic, reference) {
  if (!is.null(reference$ranges) && nrow(somatic) > 0) {
    som_gr <- GenomicRanges::GRanges(
      seqnames = somatic$chromosome,
      ranges = IRanges::IRanges(start = somatic$start, end = somatic$end)
    )

    # Harmonise chromosome naming to match somatic (bare names, no "chr" prefix)
    ref_gr <- reference$ranges
    ref_chroms <- as.character(GenomeInfoDb::seqnames(ref_gr))
    if (any(grepl("^chr", ref_chroms))) {
      GenomeInfoDb::seqlevelsStyle(ref_gr) <- "NCBI"
    }

    ol <- suppressWarnings(
      GenomicRanges::findOverlaps(ref_gr, som_gr)
    )

    somatic[, bin := NA_integer_]
    somatic[S4Vectors::subjectHits(ol), bin := S4Vectors::queryHits(ol)]
  } else {
    somatic[, bin := NA_integer_]
  }
  return(somatic)
}

#' Process Somatic VCF
#'
#' Parses a somatic VCF file, extracts allele depths, annotates with VEP CSQ fields 
#' and internal hotspots, and maps variants to reference bins.
#'
#' @param vcf_file Path to the somatic VCF file.
#' @param reference The Jumble reference object (containing $ranges).
#' @param genome Genome build version ("hg19" or "hg38").
#' @return A data.table with somatic mutations (chromosome, start, AF, bin), or NULL on failure.
#' @importFrom VariantAnnotation readVcf
#' @keywords internal
process_somatic_vcf <- function(vcf_file, reference, genome = "hg19") {
  if (is.null(vcf_file) || !file.exists(vcf_file)) {
    warning("Somatic VCF file not found or NULL.")
    return(NULL)
  }
  message("Processing somatic VCF: ", vcf_file, " (Genome: ", genome, ")")

  vcf <- tryCatch({
    suppressWarnings(VariantAnnotation::readVcf(vcf_file))
  }, error = function(e) {
    warning("Failed to read VCF: ", e$message)
    return(NULL)
  })

  if (is.null(vcf) || length(vcf) == 0) return(NULL)

  # Determine tumor index
  samples <- colnames(vcf)
  tumor_idx <- length(samples)
  if ("TUMOR" %in% toupper(samples)) {
    tumor_idx <- which(toupper(samples) == "TUMOR")[1]
  } else if (length(samples) >= 2) {
    tumor_idx <- 2
  }

  # 1. Extract Base Data
  somatic <- extract_base_variants(vcf, tumor_idx)

  # 2. Filter (AO >= 5)
  keep_mask <- somatic$AO >= 5
  somatic <- somatic[keep_mask]

  # 3. Parse VEP Annotations for Kept Rows
  csq_dt <- parse_vep_csq(vcf, keep_mask)
  somatic <- cbind(somatic, csq_dt)

  # 3.5 Aggressive Population Frequency Extirpation
  if ("MAX_AF" %in% names(somatic)) {
    pop_af <- suppressWarnings(as.numeric(somatic$MAX_AF))
    pop_af[is.na(pop_af)] <- 0
    somatic <- somatic[pop_af == 0]
  }

  # 4. Annotate Hotspots
  somatic <- annotate_hotspots(somatic, genome)

  # 5. MSI Classification (indels only)
  somatic <- classify_somatic_msi(somatic, genome)

  # 6. Map to Bins
  somatic <- map_variants_to_bins(somatic, reference)

  return(somatic)
}