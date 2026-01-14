#' Process Somatic VCF
#'
#' Parses a somatic VCF file and maps variants to reference bins.
#'
#' @param vcf_file Path to the somatic VCF file.
#' @param reference The Jumble reference object (containing $ranges).
#' @param genome Genome build version ("hg19" or "hg38").
#' @return A data.table with somatic mutations (chromosome, start, AF, bin).
#' @importFrom VariantAnnotation readVcf geno info header
#' @importFrom GenomicRanges GRanges findOverlaps subjectHits seqnames start end
#' @importFrom IRanges IRanges
#' @importFrom data.table data.table := set
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom S4Vectors subjectHits queryHits
#' @importFrom stringr str_match str_remove
#' @export
process_somatic_vcf <- function(vcf_file, reference, genome = "hg19") {
  if (is.null(vcf_file) || !file.exists(vcf_file)) {
  warning("Somatic VCF file not found or NULL.")
  return(NULL)
  }

  message("Processing somatic VCF: ", vcf_file, " (Genome: ", genome, ")")

  # Load VCF
  vcf <- tryCatch({
  suppressWarnings(VariantAnnotation::readVcf(vcf_file))
  }, error = function(e) {
  warning("Failed to read VCF: ", e$message)
  return(NULL)
  })

  if (is.null(vcf) || length(vcf) == 0) return(NULL)

  # Extract Ranges
  r <- SummarizedExperiment::rowRanges(vcf)
  chr <- as.character(GenomicRanges::seqnames(r))
  pos <- GenomicRanges::start(r)

  # Extract AF (Tumor)
  samples <- colnames(vcf)
  tumor_idx <- length(samples)
  if ("TUMOR" %in% toupper(samples)) {
  tumor_idx <- which(toupper(samples) == "TUMOR")[1]
  } else if (length(samples) >= 2) {
  tumor_idx <- 2
  }

  g <- VariantAnnotation::geno(vcf)
  af_values <- NULL

  if (!is.null(g$VAF)) {
  af_dat <- g$VAF[, tumor_idx]
  if (is.list(af_dat)) {
  af_values <- sapply(af_dat, function(x) if (length(x) > 0) as.numeric(x[1]) else 0)
  } else {
  af_values <- as.numeric(af_dat)
  }
  } else if (!is.null(g$AF)) {
  af_dat <- g$AF[, tumor_idx]
  if (is.list(af_dat)) {
  af_values <- sapply(af_dat, function(x) if (length(x) > 0) x[1] else 0)
  } else {
  af_values <- as.numeric(af_dat)
  }
  } else if (!is.null(g$AD)) {
  ad_dat <- g$AD[, tumor_idx]
  af_values <- sapply(ad_dat, function(x) {
  if (length(x) >= 2) {
  if ((x[1] + x[2]) > 0) return(x[2] / (x[1] + x[2])) else return(0)
  } else {
  return(0)
  }
  })
  } else {
  warning("No AF or AD field found in VCF. Cannot extract Allele Frequency.")
  return(NULL)
  }

  somatic <- data.table::data.table(
  chromosome = chr,
  start = pos,
  end = GenomicRanges::end(r),
  REF = as.character(r$REF),
  ALT = as.character(unlist(lapply(r$ALT, function(x) as.character(x)[1]))),
  AF = af_values
  )

  

  
  # Fix chromosome names
  set(somatic, i = NULL, j = "chromosome", value = gsub("^chr", "", somatic$chromosome))

  # AO/DP
  # AO/DP
  if (!is.null(g$AD)) {
  ad_dat <- g$AD[, tumor_idx]
  if(inherits(ad_dat, "list") || inherits(ad_dat, "IntegerList")) {
     # Pass
  }
  
  set(somatic, i = NULL, j = "AO", value = sapply(ad_dat, function(x) {
    if (length(x) >= 2 && !is.na(x[2])) x[2] else 0
  }))
  set(somatic, i = NULL, j = "DP", value = sapply(ad_dat, function(x) {
     if (length(x) > 0) sum(x, na.rm = TRUE) else 0
  }))
  } else if (!is.null(g$DP)) {
    # Check for 3D array (e.g. VarScan/Mutect sometimes output DP as 4 values: RefFwd, RefRev, AltFwd, AltRev)
    if (length(dim(g$DP)) == 3 && tumor_idx <= dim(g$DP)[2]) {
      # Slice the tumor sample: resulting in Variants x Values matrix
      dp_slice <- g$DP[, tumor_idx, ]
      if (is.matrix(dp_slice)) {
        dp_val <- rowSums(dp_slice, na.rm = TRUE)
      } else {
        # Fallback if slice is not a matrix (e.g. single variant case -> vector)
        dp_val <- sum(dp_slice, na.rm = TRUE)
      }
    } else {
      # Standard matrix or flattened logic
      dp_val <- as.numeric(g$DP)
      n_vars <- nrow(somatic)
      n_samples <- length(samples)
      
      if (length(dp_val) == n_vars) {
        # Already correct length
      } else if (length(dp_val) == n_vars * n_samples) {
        dp_mat <- matrix(dp_val, nrow = n_vars, ncol = n_samples)
        if (tumor_idx <= ncol(dp_mat)) {
          dp_val <- dp_mat[, tumor_idx]
        } else {
          dp_val <- rep(NA, n_vars)
        }
      } else {
        # Dimension mismatch
        dp_val <- rep(NA, n_vars)
      }
    }
    
    set(somatic, i = NULL, j = "DP", value = dp_val)
    set(somatic, i = NULL, j = "AO", value = as.integer(round(somatic$AF * somatic$DP)))
  } else {
    set(somatic, i = NULL, j = "AO", value = 5) # Dummy
    set(somatic, i = NULL, j = "DP", value = 10)
  }

  # Force 0 for NAs
  somatic[is.na(AO), AO := 0]
  somatic[is.na(DP), DP := 0]
  
  # Filter low support variants
  keep_mask <- somatic$AO >= 5
  somatic <- somatic[keep_mask] # Filtering rows creates a new table, careful with sync

  # VEP Parsing
  vcf_header_info <- VariantAnnotation::info(VariantAnnotation::header(vcf))
  csq_desc <- if ("CSQ" %in% rownames(vcf_header_info)) vcf_header_info$Description[rownames(vcf_header_info) == "CSQ"] else character(0)

  if (length(csq_desc) > 0) {
  csq_format <- sub(".*Format: ", "", csq_desc)
  csq_fields <- strsplit(csq_format, "\\|")[[1]]

  idx_symbol <- which(csq_fields == "SYMBOL")
  idx_cons <- which(csq_fields == "Consequence")
  idx_pos <- which(csq_fields == "Protein_position")
  idx_canonical <- which(csq_fields == "CANONICAL")
  idx_impact <- which(csq_fields == "IMPACT")
  idx_clinsig <- which(csq_fields == "CLIN_SIG")
  idx_aa <- which(csq_fields == "Amino_acids")

  # Get CSQ data and subset to match filtered somatic rows
  csq_data <- VariantAnnotation::info(vcf)$CSQ
  if (length(csq_data) == length(keep_mask)) {
  csq_data <- csq_data[keep_mask]
  } else {
  # Fallback if length mismatch (should rarely happen unless vcf header vs data mismatch)
  warning("Length mismatch between VCF info and Rows. CSQ parsing may be misaligned.")
  }

  # Parse CSQ for each variant
  parsed_csq <- lapply(csq_data, function(csq_strings) {
  if (length(csq_strings) == 0) {
  return(list(
    SYMBOL = "", Consequence = "", Protein_position = "",
    CANONICAL = "", IMPACT = "", CLIN_SIG = "", Amino_acids = ""
  ))
  }

  parsed_list <- strsplit(csq_strings, "\\|")

  # Prioritize CANONICAL='YES'
  is_canonical <- sapply(parsed_list, function(x) {
  if (length(x) >= idx_canonical) x[idx_canonical] == "YES" else FALSE
  })
  # Select best entry: Canonical, else first
  sel_idx <- if (any(is_canonical)) which(is_canonical)[1] else 1
  sel <- parsed_list[[sel_idx]]

  get_val <- function(idx) if (length(sel) >= idx) sel[idx] else ""

  list(
  SYMBOL = get_val(idx_symbol),
  Consequence = get_val(idx_cons),
  Protein_position = get_val(idx_pos),
  CANONICAL = get_val(idx_canonical),
  IMPACT = get_val(idx_impact),
  CLIN_SIG = get_val(idx_clinsig),
  Amino_acids = get_val(idx_aa)
  )
  })

  # Bind annotation columns individually to avoid parse issues
  if (length(parsed_csq) > 0) {
  set(somatic, i = NULL, j = "SYMBOL", value = sapply(parsed_csq, "[[", "SYMBOL"))
  set(somatic, i = NULL, j = "Consequence", value = sapply(parsed_csq, "[[", "Consequence"))
  set(somatic, i = NULL, j = "Protein_position", value = sapply(parsed_csq, "[[", "Protein_position"))
  set(somatic, i = NULL, j = "CANONICAL", value = sapply(parsed_csq, "[[", "CANONICAL"))
  set(somatic, i = NULL, j = "IMPACT", value = sapply(parsed_csq, "[[", "IMPACT"))
  set(somatic, i = NULL, j = "CLIN_SIG", value = sapply(parsed_csq, "[[", "CLIN_SIG"))
  set(somatic, i = NULL, j = "Amino_acids", value = sapply(parsed_csq, "[[", "Amino_acids"))
  } else {
   # Fallback if parsed_csq is empty but code reached here (e.g. 0 rows)
   # If 0 rows, we still need to add the columns
   set(somatic, i = NULL, j = "SYMBOL", value = rep(NA_character_, nrow(somatic)))
   set(somatic, i = NULL, j = "Consequence", value = rep(NA_character_, nrow(somatic)))
   set(somatic, i = NULL, j = "Protein_position", value = rep(NA_character_, nrow(somatic)))
   set(somatic, i = NULL, j = "CANONICAL", value = rep(NA_character_, nrow(somatic)))
   set(somatic, i = NULL, j = "IMPACT", value = rep(NA_character_, nrow(somatic)))
   set(somatic, i = NULL, j = "CLIN_SIG", value = rep(NA_character_, nrow(somatic)))
   set(somatic, i = NULL, j = "Amino_acids", value = rep(NA_character_, nrow(somatic)))
  }
  } else {
  # Fallback columns
  set(somatic, i = NULL, j = "SYMBOL", value = NA_character_)
  set(somatic, i = NULL, j = "Consequence", value = NA_character_)
  set(somatic, i = NULL, j = "Protein_position", value = NA_character_)
  set(somatic, i = NULL, j = "CANONICAL", value = NA_character_)
  set(somatic, i = NULL, j = "IMPACT", value = NA_character_)
  set(somatic, i = NULL, j = "CLIN_SIG", value = NA_character_)
  set(somatic, i = NULL, j = "Amino_acids", value = NA_character_)
  }

  # Hotspot Annotation ---------------------------------------------------------
  # Load internal data
  hotspots_loaded <- FALSE
  tryCatch({
  data("hotspots_snvs", package = "Jumble", envir = environment())
  data("hotspots_inframes", package = "Jumble", envir = environment())

  # Load Splice hotspots based on genome
  if (genome == "hg38") {
  data("hotspots_splice_hg38", package = "Jumble", envir = environment())
  if (exists("hotspots_splice_hg38")) {
  hotspots_splice <- hotspots_splice_hg38
  hotspots_loaded <- TRUE
  }
  } else {
  data("hotspots_splice_hg19", package = "Jumble", envir = environment())
  if (exists("hotspots_splice_hg19")) {
  hotspots_splice <- hotspots_splice_hg19
  hotspots_loaded <- TRUE
  }
  }
  
  # Check if others exist
  if (exists("hotspots_snvs") && exists("hotspots_inframes")) hotspots_loaded <- TRUE
  
  }, error = function(e) {
  warning("Hotspot data not found or failed to load. Skipping hotspot annotation.")
  })

  if (hotspots_loaded && nrow(somatic) > 0) {
  somatic[, hotkey := paste(SYMBOL, stringr::str_match(Protein_position, "^([0-9]*[-]*[0-9]*)/")[, 2])]
  somatic[SYMBOL == "" | is.na(SYMBOL), hotkey := NA_character_]
  # Clean up hotkey format "NA$"
  somatic[grepl(" NA$", hotkey), hotkey := NA_character_]

  somatic[grepl("frameshift", Consequence), hotkey := paste(hotkey, "fs")]
  somatic[grepl("inframe", Consequence), hotkey_start := stringr::str_remove(hotkey, "-[0-9]+")]
  somatic[grepl("inframe", Consequence), hotkey_end := stringr::str_remove(hotkey, "[0-9]+-")]

  somatic[, is_hotspot := FALSE]
  if (exists("hotspots_snvs")) somatic[hotkey %in% hotspots_snvs, is_hotspot := TRUE]
  if (exists("hotspots_inframes")) {
  somatic[hotkey_start %in% hotspots_inframes, is_hotspot := TRUE]
  somatic[hotkey_end %in% hotspots_inframes, is_hotspot := TRUE]
  }

  # Splice logic
  if (exists("hotspots_splice")) {
  somatic[paste(chromosome, start, sep = ":") %in% hotspots_splice, is_hotspot := TRUE]
  }

  # TERT promoter (hardcoded)
  if (genome == "hg19") {
  somatic[paste(chromosome, start, sep = ":") %in% c("5:1295228", "5:1295250"), is_hotspot := TRUE]
  } else if (genome == "hg38") {
  somatic[paste(chromosome, start, sep = ":") %in% c("5:1295113", "5:1295135"), is_hotspot := TRUE]
  }

  # Effect labeling
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

  # Map to Reference Bins

  if (!is.null(reference$ranges) && nrow(somatic) > 0) {
  som_gr <- GenomicRanges::GRanges(
  seqnames = somatic$chromosome,
  ranges = IRanges::IRanges(start = somatic$start, end = somatic$end)
  )

  ol <- GenomicRanges::findOverlaps(reference$ranges, som_gr)
  
  somatic[, bin := NA_integer_]
  somatic[S4Vectors::subjectHits(ol), bin := S4Vectors::queryHits(ol)]
  } else {
  somatic[, bin := NA_integer_]
  }

  return(somatic)
}
