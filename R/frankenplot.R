# ============================================================================
# Frankenplot — Standalone genome report for Jumble output files
# ============================================================================

# ---------------------------------------------------------------------------
# Cancer gene shortlist for Frankenplot labeling
# ---------------------------------------------------------------------------
#' Cancer gene label list for Frankenplot
#'
#' Returns a curated vector of ~100 cancer-relevant gene symbols used for
#' labeling variants in Frankenplot genome reports. This is a broader list
#' than the 20-gene label set used in standard Jumble plots.
#'
#' @return Character vector of gene symbols.
#' @keywords internal
frankenplot_cancer_genes <- function() {
  c("ABL1", "AKT1", "ALK", "APC", "AR", "ARAF", "ATM", "BAP1",
    "BIRC7", "BRAF", "BRCA1", "BRCA2", "CBFB", "CCND1", "CCND2",
    "CCND3", "CHD1", "CCNE1", "CD274", "CDK12", "CDK4", "CDK6", "CDKN2A",
    "CEBPA", "DCC", "DNMT3A", "EGFR", "ERBB2", "ERBB4", "ERCC2",
    "ERG", "ESR1", "EWSR1", "FGF3", "FGFR1", "FGFR2", "FGFR3", "FLT3",
    "GATA2", "IDH1", "IDH2", "JAK2", "KIT", "KRAS", "LRP1B", "MAP2K1",
    "MET", "MLH1", "MTOR", "MYC", "NOTCH1", "NPM1", "NRAS", "NT5C2",
    "PAX8", "PDGFRA", "PIK3CA", "PML", "PRKACA", "PTEN", "RAC1",
    "RB1", "RET", "ROS1", "RUNX1", "SF3B1", "TMPRSS2", "TP53", "U2AF1",
    "VHL", "WT1", "NTRK3", "POLE", "MSH3", "ETV6", "BRIP1", "KMT2D",
    "CTNNA1", "MSH2", "KMT2C", "FAT1", "MITF", "KMT2A", "RAF1", "BARD1",
    "CHEK2", "NF1", "CDH1", "COL5A1", "PALB2", "FH", "SPEN", "USP9X",
    "SPTA1", "MGA", "MED12", "ZFHX3", "ATRX", "BMPR1A", "ATR", "CHD3",
    "POLD1", "NOTCH2", "CREBBP", "EP300", "NCOR1", "TAF1")
}


# ---------------------------------------------------------------------------
# Frankenplot chromosome naming helpers
# ---------------------------------------------------------------------------
#' Standardize chromosome names for Frankenplot
#'
#' Removes UCSC-style chr prefixes and keeps the package-wide chromosome naming
#' convention used by Jumble outputs: 1-22, X, Y, MT/M.
#'
#' @param x Vector of chromosome names.
#' @return Character vector with standardized chromosome names.
#' @keywords internal
fp_clean_chrom_names <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x[x %in% c("23", "x")] <- "X"
  x[x %in% c("24", "y")] <- "Y"
  x[x %in% c("m", "M", "MT", "mt")] <- "MT"
  x
}

#' Standardize a chromosome column in-place
#'
#' @param dt data.table-like object.
#' @param col Column name to standardize.
#' @return The input object with standardized chromosome names.
#' @keywords internal
fp_standardize_chrom_col <- function(dt, col = "chromosome") {
  if (!is.null(dt) && col %in% names(dt)) {
    dt[, (col) := fp_clean_chrom_names(get(col))]
  }
  dt
}


# ---------------------------------------------------------------------------
# SNP VCF smart sample selection
# ---------------------------------------------------------------------------
#' Select sample columns from SNP VCF for tumor and normal
#'
#' Uses heuristic pattern matching on VCF column names to identify tumor
#' vs normal samples. Patterns: -T-, -CFDNA-, TUMOR for tumor; -N-, NORMAL
#' for normal.
#'
#' @param vcf A VCF object (VariantAnnotation).
#' @param sample_hint Optional hint for sample name (exact/fuzzy match).
#' @param role Either "tumor" or "normal" — which sample to return.
#' @return Integer index of the selected sample column.
#' @keywords internal
fp_select_snp_sample <- function(vcf, sample_hint = NULL, role = "tumor") {
  vcf_names <- colnames(vcf)
  if (length(vcf_names) <= 1) return(1L)

  # Try exact match first
  if (!is.null(sample_hint)) {
    exact <- which(vcf_names == sample_hint)
    if (length(exact) == 1) return(exact)
    fuzzy <- which(grepl(sample_hint, vcf_names, fixed = TRUE))
    if (length(fuzzy) >= 1) return(fuzzy[1])
    fuzzy_rev <- which(sapply(vcf_names, function(vn) grepl(vn, sample_hint, fixed = TRUE)))
    if (length(fuzzy_rev) >= 1) return(fuzzy_rev[1])
  }

  # Heuristic: look for -T-, -N-, -CFDNA-, TUMOR, NORMAL patterns
  n_pat <- grepl("-N-|NORMAL", vcf_names, ignore.case = TRUE)
  t_pat <- grepl("-T-|TUMOR|-CFDNA-", vcf_names, ignore.case = TRUE)

  if (role == "tumor") {
    if (sum(t_pat) == 1) return(which(t_pat))
    if (sum(n_pat) == 1) return(which(!n_pat)[1])
    return(length(vcf_names))
  } else {
    if (sum(n_pat) == 1) return(which(n_pat))
    if (sum(t_pat) == 1) return(which(!t_pat)[1])
    return(1L)
  }
}


# ---------------------------------------------------------------------------
# Parse SNP VCF for allele ratios (following original Frankenplot logic)
# ---------------------------------------------------------------------------
#' Parse SNP VCF for allele ratios
#'
#' Reads a germline SNP VCF and extracts allele ratios for biallelic SNPs
#' with rs IDs and sufficient depth. Follows the original Frankenplot logic.
#'
#' @param vcf_file Path to SNP VCF file.
#' @param sample_idx Integer index of the sample column to use.
#' @param genome Genome build ("hg19" or "hg38").
#' @return data.table with chromosome, start, end, allele_ratio, maf, depth.
#' @keywords internal
fp_parse_snp_vcf <- function(vcf_file, sample_idx, genome = "hg19") {
  if (is.null(vcf_file) || !file.exists(vcf_file)) return(NULL)

  vcf <- tryCatch(
    suppressWarnings(VariantAnnotation::readVcf(vcf_file)),
    error = function(e) { warning("Failed to read SNP VCF: ", e$message); NULL }
  )
  if (is.null(vcf) || length(vcf) == 0) return(NULL)

  # Subset to selected sample
  if (sample_idx > ncol(vcf)) sample_idx <- ncol(vcf)
  vcf_sub <- vcf[, sample_idx]

  g <- VariantAnnotation::geno(vcf_sub)
  r <- SummarizedExperiment::rowRanges(vcf_sub)
  chr <- as.character(GenomicRanges::seqnames(r))
  pos_df <- as.data.frame(GenomicRanges::ranges(r))

  if (nrow(pos_df) == 0) return(NULL)

  alf <- data.table::data.table(
    chromosome = chr,
    start = as.double(pos_df$start),
    end = as.double(pos_df$end)
  )
  alf[, pos := start]

  # Extract allele frequency
  ad <- g$AD
  dp <- g$DP
  if (is.null(ad) || is.null(dp)) return(NULL)

  alf$allele_ratio <- as.numeric(sapply(ad[, 1], "[", 2)) / as.numeric(dp[, 1])
  alf$allele_ratio[is.nan(alf$allele_ratio)] <- NA
  alf$depth <- as.numeric(dp[, 1])

  alf <- fp_standardize_chrom_col(alf)

  # Filter: biallelic SNPs with rs IDs, good depth
  is_snp <- GenomicRanges::width(r) == 1
  has_rs <- grepl("rs", names(r))
  valid_chr <- alf$chromosome %in% c(as.character(1:22), "X", "Y")

  alf <- alf[is_snp & has_rs & valid_chr]
  alf <- alf[depth > 30 & depth * allele_ratio > 5 & depth * (1 - allele_ratio) > 5]
  alf[, maf := abs(allele_ratio - 0.5) + 0.5]

  return(alf)
}


# ---------------------------------------------------------------------------
# Parse somatic VCF
# ---------------------------------------------------------------------------
#' Parse somatic VCF for Frankenplot
#'
#' Reads a somatic VCF and extracts variant information including allele
#' frequencies, VEP annotations, hotspot status, and effect classification.
#' Uses the package's bundled hotspot data via data() calls.
#'
#' @param vcf_file Path to somatic VCF file.
#' @param genome Genome build ("hg19" or "hg38").
#' @return data.table with somatic variant annotations.
#' @keywords internal
fp_parse_somatic_vcf <- function(vcf_file, genome = "hg19") {
  if (is.null(vcf_file) || !file.exists(vcf_file)) return(NULL)

  tryCatch({
    # Read VCF
    vcf <- suppressWarnings(VariantAnnotation::readVcf(vcf_file))
    if (is.null(vcf) || length(vcf) == 0) return(NULL)

    # Filter to PASS/LowQual
    filt <- SummarizedExperiment::rowRanges(vcf)$FILTER
    keep <- filt %in% c("PASS", "LowQual", ".")
    if (sum(keep) == 0) keep <- rep(TRUE, length(filt))
    vcf <- vcf[keep]
    if (length(vcf) > 50000) {
      keep2 <- SummarizedExperiment::rowRanges(vcf)$FILTER %in% c("PASS", ".")
      if (sum(keep2) > 0) vcf <- vcf[keep2]
    }

    # Determine tumor index (heuristic: normal first, tumor second)
    samples <- colnames(vcf)
    tumor_idx <- length(samples)
    if (length(samples) >= 2) {
      n_pat <- grepl("-N-|NORMAL", samples, ignore.case = TRUE)
      if (sum(n_pat) == 1 && which(n_pat) == 1) tumor_idx <- 2
      if (sum(n_pat) == 1 && which(n_pat) == 2) tumor_idx <- 1
    }

    g <- VariantAnnotation::geno(vcf)
    r <- SummarizedExperiment::rowRanges(vcf)
    chr <- as.character(GenomicRanges::seqnames(r))
    pos_df <- as.data.frame(GenomicRanges::ranges(r))

    if (nrow(pos_df) == 0) return(NULL)

    salf <- data.table::data.table(
      chromosome = chr,
      start = as.double(pos_df$start),
      end = as.double(pos_df$end)
    )
    salf <- fp_standardize_chrom_col(salf)
    salf$REF <- as.character(VariantAnnotation::ref(vcf))
    alt_list <- VariantAnnotation::alt(vcf)
    salf$ALT <- vapply(alt_list, function(x) as.character(x)[1], character(1))
    salf$FILTER <- SummarizedExperiment::rowRanges(vcf)$FILTER

    # Extract AF
    if (!is.null(g$VAF)) {
      salf[, `AF.T` := as.numeric(g$VAF[, tumor_idx])]
    } else if (!is.null(g$AF)) {
      salf[, `AF.T` := as.numeric(g$AF[, tumor_idx])]
    } else if (!is.null(g$AD) && !is.null(g$DP)) {
      ad_alt <- sapply(g$AD[, tumor_idx], "[", 2)
      salf[, `AF.T` := as.numeric(ad_alt) / as.numeric(g$DP[, tumor_idx])]
    } else {
      salf[, `AF.T` := NA_real_]
    }

    # Extract AO and DP
    if (!is.null(g$AD)) {
      salf[, `AO.T` := as.numeric(sapply(g$AD[, tumor_idx], "[", 2))]
      salf[, `DP.T` := as.numeric(g$DP[, tumor_idx])]
    } else if (!is.null(g$DP4)) {
      if (nrow(salf) > 1) {
        salf[, `AO.T` := as.numeric(apply(g$DP4[, tumor_idx, 3:4], 1, sum))]
        salf[, `DP.T` := as.numeric(apply(g$DP4[, tumor_idx, ], 1, sum))]
      } else {
        salf[, `AO.T` := NA_real_]
        salf[, `DP.T` := NA_real_]
      }
    } else {
      salf[, `AO.T` := NA_real_]
      salf[, `DP.T` := NA_real_]
    }

    # Variant type
    salf[, type := "other"]
    salf[nchar(REF) == 1 & nchar(ALT) == 1, type := "snv"]
    salf[nchar(REF) > nchar(ALT), type := "del"]
    salf[nchar(REF) < nchar(ALT), type := "ins"]
    salf[, `point mutation` := type]

    # Parse VEP CSQ annotations
    salf <- fp_parse_vep_csq(vcf, salf)

    # Filter: AO >= 5 when allele observation counts are available.
    # Tumor-only or WIP caller outputs may provide AF without AO/AD; keep those
    # rows so they can still be displayed in Frankenplot.
    if (any(!is.na(salf$`AO.T`))) {
      salf <- salf[!is.na(`AO.T`) & `AO.T` >= 5]
    }

    # Annotate hotspots using package data
    salf <- fp_annotate_hotspots(salf, genome)

    # Annotate effect
    salf <- fp_annotate_effect(salf)

    return(salf)
  }, error = function(e) {
    warning("Failed to parse somatic VCF: ", e$message)
    return(NULL)
  })
}


# ---------------------------------------------------------------------------
# Parse VEP CSQ from VCF info field
# ---------------------------------------------------------------------------
#' Parse VEP CSQ annotations from VCF
#'
#' Extracts Consequence annotations from the VCF INFO CSQ field (VEP format).
#' Picks the canonical transcript annotation when available.
#'
#' @param vcf A VCF object.
#' @param dt data.table to add annotation columns to.
#' @return The input data.table with added VEP annotation columns.
#' @keywords internal
fp_parse_vep_csq <- function(vcf, dt) {
  key_cols <- c("SYMBOL", "Consequence", "IMPACT", "CANONICAL",
                "Protein_position", "Amino_acids", "CLIN_SIG",
                "Existing_variation", "SIFT", "PolyPhen", "BIOTYPE")

  tryCatch({
    hdr <- VariantAnnotation::info(VariantAnnotation::header(vcf))
    desc <- hdr$Description
    ix <- grep("Consequence annotations from Ensembl", desc)
    if (length(ix) == 0) {
      for (col in key_cols) dt[, (col) := ""]
      return(dt)
    }

    header_str <- strsplit(desc[ix], "\\|")[[1]]
    header_str[1] <- "Allele"

    vep <- VariantAnnotation::info(vcf)$CSQ
    if (is.null(vep)) {
      for (col in key_cols) dt[, (col) := ""]
      return(dt)
    }

    # For each variant, pick the first canonical or first annotation
    n_vars <- length(vep)
    result_list <- vector("list", n_vars)

    for (i in seq_len(n_vars)) {
      entries <- vep[[i]]
      if (length(entries) == 0) {
        result_list[[i]] <- rep("", length(header_str))
        next
      }
      parsed <- lapply(entries, function(e) strsplit(e, "\\|")[[1]])
      canon_idx <- which(sapply(parsed, function(p) {
        if (length(p) >= which(header_str == "CANONICAL")) {
          p[which(header_str == "CANONICAL")] == "YES"
        } else FALSE
      }))
      if (length(canon_idx) >= 1) {
        chosen <- parsed[[canon_idx[1]]]
      } else {
        chosen <- parsed[[1]]
      }
      if (length(chosen) < length(header_str)) {
        chosen <- c(chosen, rep("", length(header_str) - length(chosen)))
      }
      result_list[[i]] <- chosen[seq_along(header_str)]
    }

    csq_mat <- do.call(rbind, result_list)
    colnames(csq_mat) <- header_str

    for (col in key_cols) {
      if (col %in% colnames(csq_mat)) {
        dt[, (col) := csq_mat[, col]]
      } else {
        dt[, (col) := ""]
      }
    }

    return(dt)
  }, error = function(e) {
    for (col in key_cols) dt[, (col) := ""]
    return(dt)
  })
}


# ---------------------------------------------------------------------------
# Annotate hotspots — uses package bundled hotspot data
# ---------------------------------------------------------------------------
#' Annotate hotspot status using package hotspot data
#'
#' Uses the package's bundled hotspot_snvs, hotspots_inframes, and
#' hotspots_splice data loaded via data() calls, matching the approach
#' in annotate_hotspots() from somatic.R.
#'
#' @param dt data.table with SYMBOL, Protein_position, Consequence columns.
#' @param genome Genome build ("hg19" or "hg38").
#' @return The input data.table with is_hotspot column added.
#' @keywords internal
fp_annotate_hotspots <- function(dt, genome = "hg19") {
  if (nrow(dt) == 0) return(dt)

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
    warning("Hotspot data not found. Skipping hotspot annotation.")
  })

  dt[, hotkey := paste(SYMBOL, stringr::str_match(Protein_position, "^([0-9]*[-]*[0-9]*)/")[, 2])]
  dt[SYMBOL == "" | is.na(SYMBOL), hotkey := NA_character_]
  dt[grepl(" NA$", hotkey), hotkey := NA_character_]

  dt[, is_hotspot := FALSE]

  if (hotspots_loaded) {
    if (exists("hotspots_snvs")) dt[hotkey %in% hotspots_snvs, is_hotspot := TRUE]
    if (exists("hotspots_inframes")) {
      dt[grepl("inframe", Consequence), hotkey_start := stringr::str_remove(hotkey, "-[0-9]+")]
      dt[grepl("inframe", Consequence), hotkey_end := stringr::str_remove(hotkey, "[0-9]+-")]
      dt[hotkey_start %in% hotspots_inframes | hotkey_end %in% hotspots_inframes, is_hotspot := TRUE]
    }
    if (exists("hotspots_splice")) {
      dt[paste(chromosome, start, sep = ":") %in% hotspots_splice, is_hotspot := TRUE]
    }
    # TERT promoter
    tert_coords <- if (genome == "hg19") c("5:1295228", "5:1295250") else c("5:1295113", "5:1295135")
    dt[paste(chromosome, start, sep = ":") %in% tert_coords, is_hotspot := TRUE]
  }

  return(dt)
}


# ---------------------------------------------------------------------------
# Annotate effect classification
# ---------------------------------------------------------------------------
#' Classify variant effect
#'
#' Classifies variants as uncertain, high-impact, or hotspot based on
#' VEP IMPACT, ClinVar significance, and hotspot status.
#'
#' @param dt data.table with IMPACT, CLIN_SIG, CANONICAL, SYMBOL, is_hotspot.
#' @return The input data.table with effect column added.
#' @keywords internal
fp_annotate_effect <- function(dt) {
  if (nrow(dt) == 0) return(dt)

  dt[, effect := NA_character_]
  ix <- dt$CANONICAL == "YES" | dt$SYMBOL != ""
  dt[ix & IMPACT == "MODERATE", effect := "uncertain"]
  dt[ix & IMPACT == "HIGH", effect := "high-impact"]
  dt[ix & grepl("pathogenic", CLIN_SIG), effect := "high-impact"]
  dt[ix & is_hotspot == TRUE, effect := "hotspot"]

  return(dt)
}


# ---------------------------------------------------------------------------
# Parse germline VCF
# ---------------------------------------------------------------------------
#' Parse germline VCF for Frankenplot
#'
#' Reads a germline VCF, pre-filters to HIGH impact or pathogenic variants,
#' extracts normal and tumor allele fractions, and filters to cancer genes.
#' Returns two data.tables: galf_n (normal AF) and galf_t (tumor AF for
#' blue overlay in tumor panels).
#'
#' @param vcf_file Path to germline VCF file.
#' @param tumor_sample_name Optional tumor sample name hint.
#' @param normal_sample_name Optional normal sample name hint.
#' @param genome Genome build ("hg19" or "hg38").
#' @return List with galf_n and galf_t data.tables.
#' @keywords internal
fp_parse_germline_vcf <- function(vcf_file, tumor_sample_name = NULL,
                                   normal_sample_name = NULL, genome = "hg19") {
  if (is.null(vcf_file) || !file.exists(vcf_file)) return(list(galf_n = NULL, galf_t = NULL))

  tryCatch({
    vcf <- suppressWarnings(VariantAnnotation::readVcf(vcf_file))
    if (is.null(vcf) || length(vcf) == 0) return(list(galf_n = NULL, galf_t = NULL))

    samples <- colnames(vcf)

    # Determine normal and tumor indices — original logic:
    # If -N- or NORMAL is in column 2, swap so normal is column 1
    normal_idx <- 1L
    tumor_idx <- NULL

    if (length(samples) >= 2) {
      n_pat <- grepl("-N-|NORMAL", samples, ignore.case = TRUE)
      t_pat <- grepl("-T-|TUMOR|-CFDNA-", samples, ignore.case = TRUE)

      if (!is.null(normal_sample_name)) {
        m <- which(grepl(normal_sample_name, samples, fixed = TRUE))
        if (length(m) >= 1) normal_idx <- m[1]
      } else if (sum(n_pat) == 1) {
        normal_idx <- which(n_pat)
      }

      if (!is.null(tumor_sample_name)) {
        m <- which(grepl(tumor_sample_name, samples, fixed = TRUE))
        if (length(m) >= 1) tumor_idx <- m[1]
      } else if (sum(t_pat) == 1) {
        tumor_idx <- which(t_pat)
      } else {
        tumor_idx <- setdiff(seq_along(samples), normal_idx)[1]
      }
    }

    # Pre-filter to HIGH impact or pathogenic (original: csq filtering)
    csq_info <- tryCatch(VariantAnnotation::info(vcf)$CSQ, error = function(e) NULL)
    if (!is.null(csq_info)) {
      keep <- sapply(csq_info, function(entries) {
        any(grepl("HIGH|pathogenic", entries))
      })
      if (sum(keep) > 0) vcf <- vcf[keep]
    }

    if (length(vcf) == 0) return(list(galf_n = NULL, galf_t = NULL))

    # Extract normal sample data
    g <- VariantAnnotation::geno(vcf)
    r <- SummarizedExperiment::rowRanges(vcf)
    chr <- as.character(GenomicRanges::seqnames(r))
    pos_df <- as.data.frame(GenomicRanges::ranges(r))

    if (nrow(pos_df) == 0) return(list(galf_n = NULL, galf_t = NULL))

    galf <- data.table::data.table(
      chromosome = chr,
      start = as.double(pos_df$start),
      end = as.double(pos_df$end)
    )
    galf <- fp_standardize_chrom_col(galf)
    galf$REF <- as.character(VariantAnnotation::ref(vcf))
    alt_list <- VariantAnnotation::alt(vcf)
    galf$ALT <- vapply(alt_list, function(x) as.character(x)[1], character(1))

    # Normal AF (from normal_idx column)
    if (!is.null(g$AD) && !is.null(g$DP)) {
      galf$AF <- as.numeric(sapply(g$AD[, normal_idx], "[", 2)) / as.numeric(g$DP[, normal_idx])
      galf$AO <- as.numeric(sapply(g$AD[, normal_idx], "[", 2))
      galf$DP <- as.numeric(g$DP[, normal_idx])
    } else {
      galf[, AF := NA_real_]
      galf[, AO := NA_real_]
      galf[, DP := NA_real_]
    }

    # Tumor AF (for blue overlay — from tumor_idx column)
    galf[, AF_t := NA_real_]
    if (!is.null(tumor_idx) && !is.null(g$AD) && !is.null(g$DP)) {
      galf$AF_t <- as.numeric(sapply(g$AD[, tumor_idx], "[", 2)) / as.numeric(g$DP[, tumor_idx])
    }

    # Variant type
    galf[, type := "other"]
    galf[nchar(REF) == 1 & nchar(ALT) == 1, type := "snv"]
    galf[nchar(REF) > nchar(ALT), type := "del"]
    galf[nchar(REF) < nchar(ALT), type := "ins"]
    galf[, `point mutation` := type]

    # Parse VEP
    galf <- fp_parse_vep_csq(vcf, galf)

    # Filter: AF > 0
    galf <- galf[!is.na(AF) & AF > 0]

    # Annotate effect (following original logic with SIFT/PolyPhen)
    galf[, effect := NA_character_]
    ix <- galf$CANONICAL == "YES" | galf$SYMBOL != ""
    galf[ix & IMPACT == "MODERATE", effect := "uncertain"]
    galf[ix & IMPACT == "HIGH", effect := "high-impact"]
    galf[ix & grepl("pathogenic", CLIN_SIG), effect := "high-impact"]
    if ("SIFT" %in% names(galf)) {
      galf[ix & grepl("deleterious", SIFT), effect := "high-impact"]
    }
    if ("PolyPhen" %in% names(galf)) {
      galf[ix & grepl("damaging", PolyPhen), effect := "high-impact"]
    }

    # gnomAD frequency filtering (original logic)
    if (!is.null(galf$gnomADe_AF)) galf[, gnomAD_AF := as.numeric(gnomADe_AF)]
    if (!is.null(galf$gnomADg_AF)) galf[, gnomAD_AF := as.numeric(gnomADg_AF)]
    if (!is.null(galf$gnomAD_AF) && !is.numeric(galf$gnomAD_AF)) {
      galf[, gnomAD_AF := as.numeric(gnomAD_AF)]
    }
    if (!"gnomAD_AF" %in% names(galf)) galf[, gnomAD_AF := 0]
    galf[is.na(gnomAD_AF), gnomAD_AF := 0]

    # Filter to cancer genes with low gnomAD frequency
    cancer_genes <- frankenplot_cancer_genes()
    galf_n <- galf[gnomAD_AF < 0.01 & SYMBOL %in% cancer_genes]
    galf_n[, allele_ratio := AF]

    galf_t <- data.table::copy(galf[gnomAD_AF < 0.01 & SYMBOL %in% cancer_genes])
    galf_t[, allele_ratio := AF_t]

    return(list(galf_n = galf_n, galf_t = galf_t))
  }, error = function(e) {
    warning("Failed to parse germline VCF: ", e$message)
    return(list(galf_n = NULL, galf_t = NULL))
  })
}


# ---------------------------------------------------------------------------
# Read Jumble CSV (targets)
# ---------------------------------------------------------------------------
#' Read Jumble CSV output for Frankenplot
#'
#' Reads a .jumble.csv file and prepares it for Frankenplot visualization.
#' Adds bin numbering, clamps outliers, computes smoothed log2, and
#' classifies target types.
#'
#' @param csv_file Path to .jumble.csv file.
#' @return data.table with bin, log2, smooth_log2, target type columns.
#' @keywords internal
fp_read_jumble_csv <- function(csv_file) {
  if (is.null(csv_file) || !file.exists(csv_file)) {
    stop("Jumble CSV file not found: ", csv_file)
  }

  bins <- data.table::fread(csv_file)
  bins[, bin := 1:.N]

  # Create 'depth' alias from 'count' if needed (original Frankenplot uses 'depth')
  if ("count" %in% names(bins) && !"depth" %in% names(bins)) {
    bins[, depth := count]
  }

  bins <- fp_standardize_chrom_col(bins)
  if ("seqnames" %in% names(bins)) bins <- fp_standardize_chrom_col(bins, "seqnames")

  if ("gene" %in% names(bins)) {
    bins[gene == "Antitarget", gene := "Background"]
  }

  # Fix outliers (original logic)
  bins[2^log2 < 0.05, log2 := log2(0.05)]
  bins[2^log2 > 100, log2 := log2(100)]

  # Smooth
  if (nrow(bins) >= 25) {
    bins[, smooth_log2 := stats::runmed(log2, k = 25)]
  } else {
    bins[, smooth_log2 := log2]
  }

  # Classify target type based on spacing
  d <- bins[, (end + start + 1) / 2]
  d <- abs(diff(d))
  d_dt <- data.table::data.table(left = c(1e6, d), right = c(d, 1e6))
  d_dt[, min_d := left][right < left, min_d := right]
  bins[, `target type` := "target"]
  bins[d_dt$min_d < 400, `target type` := "exon"]
  if ("gene" %in% names(bins)) {
    bins[gene == "Background", `target type` := "background"]
  }

  return(bins)
}


# ---------------------------------------------------------------------------
# Read CNS (segments)
# ---------------------------------------------------------------------------
#' Read CNS segments for Frankenplot
#'
#' Reads a .cns file and maps segments to bins via genomic overlap.
#' Updates segment start/end to bin indices and assigns segment IDs to bins.
#'
#' @param cns_file Path to .cns file.
#' @param bins data.table of bins (from fp_read_jumble_csv).
#' @return List with updated bins and segments data.tables.
#' @keywords internal
fp_read_cns <- function(cns_file, bins) {
  if (is.null(cns_file) || !file.exists(cns_file)) {
    stop("CNS file not found: ", cns_file)
  }

  segments <- data.table::fread(cns_file)
  segments <- fp_standardize_chrom_col(segments)
  segments[, start_pos := start][, end_pos := end]

  # Map segments to bins via overlap
  binranges <- GenomicRanges::makeGRangesFromDataFrame(
    bins[, .(chromosome, start, end)],
    seqnames.field = "chromosome"
  )
  segmentranges <- GenomicRanges::makeGRangesFromDataFrame(
    segments[, .(chromosome, start, end)],
    seqnames.field = "chromosome"
  )
  overlap <- data.table::as.data.table(GenomicRanges::findOverlaps(segmentranges, binranges))
  starts <- overlap[, min(subjectHits), by = queryHits]
  ends <- overlap[, max(subjectHits), by = queryHits]
  segments[, start := NA_integer_][starts$queryHits, start := starts$V1]
  segments[, end := NA_integer_][ends$queryHits, end := ends$V1]
  bins[overlap$subjectHits, segment := overlap$queryHits]

  return(list(bins = bins, segments = segments))
}


# ---------------------------------------------------------------------------
# Map SNPs to bins
# ---------------------------------------------------------------------------
#' Map SNP allele ratios to bins
#'
#' Maps SNP VCF allele ratios to Jumble bins via genomic overlap (with
#' 50bp padding). Adds allele_ratio and maf columns to bins.
#'
#' @param alf data.table of SNP allele ratios (from fp_parse_snp_vcf).
#' @param bins data.table of bins.
#' @return List with updated bins and alf data.tables.
#' @keywords internal
fp_map_snps_to_bins <- function(alf, bins) {
  bins <- fp_standardize_chrom_col(bins)
  if (is.null(alf) || nrow(alf) == 0) {
    bins[, allele_ratio := as.numeric(NA)]
    bins[, maf := as.numeric(NA)]
    return(list(bins = bins, alf = alf))
  }

  alf <- fp_standardize_chrom_col(alf)

  binranges <- GenomicRanges::makeGRangesFromDataFrame(
    bins[, .(chromosome, start = start - 50, end = end + 50)],
    seqnames.field = "chromosome"
  )
  snpranges <- GenomicRanges::makeGRangesFromDataFrame(
    alf[, .(chromosome, start, end)],
    seqnames.field = "chromosome"
  )
  overlap <- GenomicRanges::findOverlaps(snpranges, binranges)
  alf[S4Vectors::queryHits(overlap), bin := bins[S4Vectors::subjectHits(overlap)]$bin]

  bins[, allele_ratio := as.numeric(NA)]
  bins[, maf := as.numeric(NA)]
  alf_valid <- alf[!is.na(allele_ratio)]
  if (nrow(alf_valid) > 0) {
    bins[alf_valid$bin, allele_ratio := alf_valid$allele_ratio]
    bins[alf_valid$bin, maf := alf_valid$maf]
  }

  return(list(bins = bins, alf = alf))
}


# ---------------------------------------------------------------------------
# Map mutations to bins
# ---------------------------------------------------------------------------
#' Map somatic or germline mutations to bins
#'
#' Maps variant positions to Jumble bins via genomic overlap (50bp padding).
#' Adds a bin column to the variant data.table.
#'
#' @param variants data.table with chromosome, start, end columns.
#' @param bins data.table of bins.
#' @return The variant data.table with bin column added.
#' @keywords internal
fp_map_variants_to_bins <- function(variants, bins) {
  if (is.null(variants) || nrow(variants) == 0) return(variants)

  variants <- fp_standardize_chrom_col(variants)
  bins <- fp_standardize_chrom_col(bins)

  binranges <- GenomicRanges::makeGRangesFromDataFrame(
    bins[, .(chromosome, start = start - 50, end = end + 50)],
    seqnames.field = "chromosome"
  )
  mutranges <- GenomicRanges::makeGRangesFromDataFrame(
    variants[, .(chromosome, start, end)],
    seqnames.field = "chromosome"
  )
  overlap <- GenomicRanges::findOverlaps(mutranges, binranges)
  variants[S4Vectors::queryHits(overlap), bin := bins[S4Vectors::subjectHits(overlap)]$bin]

  return(variants)
}


# ---------------------------------------------------------------------------
# Parse DPYD files
# ---------------------------------------------------------------------------
#' Parse DPYD genotyping files
#'
#' Reads DPYD JSON result and/or CSV evidence files. JSON is read via
#' RJSONIO::readJSONStream if available, otherwise jsonlite::fromJSON.
#'
#' @param dpyd_json Path to DPYD JSON result file (optional).
#' @param dpyd_csv Path to DPYD CSV evidence file (optional).
#' @return List with dpyd_result (data.table) and dpyd_table (data.table).
#' @keywords internal
fp_parse_dpyd <- function(dpyd_json = NULL, dpyd_csv = NULL) {
  dpyd_result <- NULL
  dpyd_table <- NULL

  if (!is.null(dpyd_json) && file.exists(dpyd_json)) {
    tryCatch({
      if (requireNamespace("RJSONIO", quietly = TRUE)) {
        dpyd_result <- data.table::as.data.table(
          t(RJSONIO::readJSONStream(dpyd_json))
        )
      } else if (requireNamespace("jsonlite", quietly = TRUE)) {
        j <- jsonlite::fromJSON(dpyd_json)
        dpyd_result <- data.table::as.data.table(t(j))
      }
    }, error = function(e) {
      warning("Failed to read DPYD JSON: ", e$message)
    })
  }

  if (!is.null(dpyd_csv) && file.exists(dpyd_csv)) {
    tryCatch({
      dpyd_table <- data.table::fread(dpyd_csv)
      if ("file" %in% names(dpyd_table)) {
        dpyd_table[, file := basename(file)]
      }
    }, error = function(e) {
      warning("Failed to read DPYD CSV: ", e$message)
    })
  }

  return(list(dpyd_result = dpyd_result, dpyd_table = dpyd_table))
}


# ===========================================================================
# Main exported function
# ===========================================================================

#' Generate Frankenplot Genome Report
#'
#' Creates an interactive HTML genome report from Jumble output files.
#' The report includes genome-wide and per-chromosome copy number plots
#' with SNP allele ratios, somatic mutation overlays (red), germline
#' mutation overlays (blue), GIS curves, and optional DPYD genotyping.
#'
#' @param tumor_jumble_csv Path to tumor .jumble.csv file (required).
#' @param tumor_cns Path to tumor .cns segment file (required).
#' @param output_file Path for the output HTML report (required).
#' @param normal_jumble_csv Path to normal .jumble.csv file (optional).
#' @param normal_cns Path to normal .cns segment file (optional).
#' @param tumor_snp_vcf Path to tumor SNP VCF file (optional).
#' @param normal_snp_vcf Path to normal SNP VCF file (optional).
#' @param somatic_vcf Path to somatic mutation VCF file (optional).
#' @param germline_vcf Path to germline mutation VCF file (optional).
#' @param tumor_sample_name Tumor sample name hint for VCF column selection.
#' @param normal_sample_name Normal sample name hint for VCF column selection.
#' @param genome Genome build version ("hg19" or "hg38"). Auto-detected if NULL.
#' @param hrdtable Precomputed GIS/HRD table to display in the report (optional).
#'   Frankenplot does not compute GIS internally.
#' @param dpyd_json Path to DPYD JSON result file (optional).
#' @param dpyd_csv Path to DPYD CSV evidence file (optional).
#'
#' @return Invisibly returns the path to the generated HTML report.
#'
#' @details
#' The function reads Jumble output files (targets and segments), optionally
#' parses VCF files for SNP allele ratios, somatic mutations, and germline
#' mutations, and renders an interactive HTML report using an Rmd template.
#' GIS/HRD results are displayed only when a precomputed table is supplied via
#' \code{hrdtable}; Frankenplot does not compute GIS internally.
#'
#' The report structure follows the original Frankenplot design:
#' \itemize{
#'   \item \strong{Genome tab}: Overview plots (log2 ratio, allele ratio, depth)
#'   \item \strong{HRD tab}: GIS score curves across tumor fractions
#'   \item \strong{DPYD tab}: DPYD genotyping results (if provided)
#'   \item \strong{Chromosome tabs}: Per-chromosome detail plots
#'   \item \strong{Data tab}: Somatic variant table and session info
#' }
#'
#' @examples
#' \dontrun{
#' frankenplot(
#'   tumor_jumble_csv = "sample.jumble.csv",
#'   tumor_cns = "sample.cns",
#'   output_file = "sample_frankenplot.html",
#'   tumor_snp_vcf = "sample.snp.vcf.gz",
#'   somatic_vcf = "sample.somatic.vcf.gz"
#' )
#' }
#'
#' @export
frankenplot <- function(tumor_jumble_csv,
                        tumor_cns,
                        output_file,
                        normal_jumble_csv = NULL,
                        normal_cns = NULL,
                        tumor_snp_vcf = NULL,
                        normal_snp_vcf = NULL,
                        somatic_vcf = NULL,
                        germline_vcf = NULL,
                        tumor_sample_name = NULL,
                        normal_sample_name = NULL,
                        genome = NULL,
                        hrdtable = NULL,
                        qc_file = NULL,
                        qc_metrics = NULL,
                        dpyd_json = NULL,
                        dpyd_csv = NULL) {

  # ---- Input validation ----
  if (missing(tumor_jumble_csv) || !file.exists(tumor_jumble_csv)) {
    stop("tumor_jumble_csv is required and must exist.")
  }
  if (missing(tumor_cns) || !file.exists(tumor_cns)) {
    stop("tumor_cns is required and must exist.")
  }
  if (missing(output_file)) {
    stop("output_file is required.")
  }

  # Validate optional file arguments
  opt_files <- list(
    normal_jumble_csv = normal_jumble_csv,
    normal_cns = normal_cns,
    tumor_snp_vcf = tumor_snp_vcf,
    normal_snp_vcf = normal_snp_vcf,
    somatic_vcf = somatic_vcf,
    germline_vcf = germline_vcf,
    qc_file = qc_file,
    dpyd_json = dpyd_json,
    dpyd_csv = dpyd_csv
  )
  for (nm in names(opt_files)) {
    f <- opt_files[[nm]]
    if (!is.null(f) && !file.exists(f)) {
      stop(nm, " file not found: ", f)
    }
  }

  # Normal requires both CSV and CNS
  has_normal <- !is.null(normal_jumble_csv) && !is.null(normal_cns)

  # ---- Auto-detect genome ----
  if (is.null(genome)) {
    # Try to detect from explicit seqnames first, then fall back to chromosome
    # naming. Bare chromosome names remain ambiguous, so preserve historical
    # behavior and default to hg19 unless the caller specifies hg38.
    test_dt <- data.table::fread(tumor_jumble_csv, nrows = 100)
    chrom_source <- if ("seqnames" %in% names(test_dt)) test_dt$seqnames else test_dt$chromosome
    if (any(grepl("^chr", chrom_source))) {
      genome <- "hg38"
    } else {
      genome <- "hg19"
    }
    message("Auto-detected genome: ", genome)
  }

  message("=== Frankenplot Report Generation ===")
  message("Tumor CSV: ", tumor_jumble_csv)
  message("Tumor CNS: ", tumor_cns)
  if (has_normal) message("Normal CSV: ", normal_jumble_csv)
  message("Genome: ", genome)

  # ---- 1. Read tumor copy number data ----
  message("Reading tumor copy number data...")
  bins_t <- fp_read_jumble_csv(tumor_jumble_csv)
  res_t <- fp_read_cns(tumor_cns, bins_t)
  bins_t <- res_t$bins
  segments_t <- res_t$segments

  # ---- 2. Parse SNP VCF (tumor) ----
  alf <- NULL
  if (!is.null(tumor_snp_vcf)) {
    message("Parsing tumor SNP VCF...")
    snp_vcf <- tryCatch(
      suppressWarnings(VariantAnnotation::readVcf(tumor_snp_vcf)),
      error = function(e) NULL
    )
    if (!is.null(snp_vcf)) {
      t_idx <- fp_select_snp_sample(snp_vcf, tumor_sample_name, "tumor")
      alf <- fp_parse_snp_vcf(tumor_snp_vcf, t_idx, genome)
    }
  }

  # Map SNPs to tumor bins
  snp_res <- fp_map_snps_to_bins(alf, bins_t)
  bins_t <- snp_res$bins
  alf <- snp_res$alf

  # ---- 3. Parse somatic VCF ----
  salf <- NULL
  if (!is.null(somatic_vcf)) {
    message("Parsing somatic VCF...")
    salf <- fp_parse_somatic_vcf(somatic_vcf, genome)
  }

  # Map somatic mutations to tumor bins
  if (!is.null(salf) && nrow(salf) > 0) {
    salf <- fp_map_variants_to_bins(salf, bins_t)
  }

  # ---- 4. Parse germline VCF ----
  galf_n <- NULL
  galf_t <- NULL
  if (!is.null(germline_vcf)) {
    message("Parsing germline VCF...")
    germ <- fp_parse_germline_vcf(germline_vcf, tumor_sample_name,
                                   normal_sample_name, genome)
    galf_n <- germ$galf_n
    galf_t <- germ$galf_t
  }

  # Map germline mutations to tumor bins (for blue overlay)
  if (!is.null(galf_t) && nrow(galf_t) > 0) {
    galf_t <- fp_map_variants_to_bins(galf_t, bins_t)
  }

  # ---- 5. Read normal copy number data ----
  bins_n <- NULL
  segments_n <- NULL
  alf_n <- NULL
  if (has_normal) {
    message("Reading normal copy number data...")
    bins_n <- fp_read_jumble_csv(normal_jumble_csv)
    res_n <- fp_read_cns(normal_cns, bins_n)
    bins_n <- res_n$bins
    segments_n <- res_n$segments

    # Parse normal SNP VCF
    if (!is.null(normal_snp_vcf)) {
      message("Parsing normal SNP VCF...")
      snp_vcf_n <- tryCatch(
        suppressWarnings(VariantAnnotation::readVcf(normal_snp_vcf)),
        error = function(e) NULL
      )
      if (!is.null(snp_vcf_n)) {
        n_idx <- fp_select_snp_sample(snp_vcf_n, normal_sample_name, "normal")
        alf_n <- fp_parse_snp_vcf(normal_snp_vcf, n_idx, genome)
      }
    }

    # Map SNPs to normal bins
    snp_res_n <- fp_map_snps_to_bins(alf_n, bins_n)
    bins_n <- snp_res_n$bins

    # Map germline mutations to normal bins
    if (!is.null(galf_n) && nrow(galf_n) > 0) {
      galf_n <- fp_map_variants_to_bins(galf_n, bins_n)
    }
  }

  # ---- 6. Use externally supplied GIS/HRD results if available ----
  if (!is.null(hrdtable)) {
    hrdtable <- data.table::as.data.table(hrdtable)
    message("Using externally supplied GIS/HRD table.")
  } else {
    message("No GIS/HRD table supplied; HRD tab will report missing GIS data.")
  }

  # ---- 7. Read QC metrics ----
  if (is.null(qc_metrics)) {
    if (is.null(qc_file)) {
      qc_candidate <- sub("\\.jumble\\.csv$", ".qc.csv", tumor_jumble_csv)
      if (!identical(qc_candidate, tumor_jumble_csv) && file.exists(qc_candidate)) {
        qc_file <- qc_candidate
      }
    }
    if (!is.null(qc_file) && file.exists(qc_file)) {
      message("Reading QC metrics: ", qc_file)
      qc_metrics <- tryCatch(
        data.table::fread(qc_file),
        error = function(e) {
          warning("Failed to read QC metrics: ", e$message)
          NULL
        }
      )
    }
  } else {
    qc_metrics <- data.table::as.data.table(qc_metrics)
  }

  # ---- 8. Parse DPYD ----
  dpyd_data <- fp_parse_dpyd(dpyd_json, dpyd_csv)

  # ---- 9. Render report ----
  message("Rendering Frankenplot report...")
  fp_render_report(
    bins_t = bins_t,
    segments_t = segments_t,
    bins_n = bins_n,
    segments_n = segments_n,
    alf = alf,
    alf_n = alf_n,
    salf = salf,
    galf_t = galf_t,
    galf_n = galf_n,
    hrdtable = hrdtable,
    qc_metrics = qc_metrics,
    dpyd_result = dpyd_data$dpyd_result,
    dpyd_table = dpyd_data$dpyd_table,
    genome = genome,
    output_file = output_file,
    tumor_jumble_csv = tumor_jumble_csv
  )

  message("Frankenplot report generated: ", output_file)
  invisible(output_file)
}


# ---------------------------------------------------------------------------
# Render Frankenplot report
# ---------------------------------------------------------------------------
#' Render Frankenplot HTML report
#'
#' Internal function that renders the Frankenplot Rmd template with all
#' pre-processed data passed as parameters.
#'
#' @param bins_t Tumor bins data.table.
#' @param segments_t Tumor segments data.table.
#' @param bins_n Normal bins data.table (or NULL).
#' @param segments_n Normal segments data.table (or NULL).
#' @param alf Tumor SNP allele ratio data.table (or NULL).
#' @param alf_n Normal SNP allele ratio data.table (or NULL).
#' @param salf Somatic variant data.table (or NULL).
#' @param galf_t Germline variants with tumor AF (or NULL).
#' @param galf_n Germline variants with normal AF (or NULL).
#' @param hrdtable GIS table data.table (or NULL).
#' @param qc_metrics QC metrics data.table (or NULL).
#' @param dpyd_result DPYD result data.table (or NULL).
#' @param dpyd_table DPYD evidence data.table (or NULL).
#' @param genome Genome build string.
#' @param output_file Output HTML path.
#' @param tumor_jumble_csv Path to tumor CSV (for title).
#' @keywords internal
fp_render_report <- function(bins_t, segments_t, bins_n, segments_n,
                              alf, alf_n, salf, galf_t, galf_n,
                              hrdtable, qc_metrics, dpyd_result, dpyd_table,
                              genome, output_file, tumor_jumble_csv) {

  # Find the Rmd template
  template_path <- system.file(
    "rmarkdown", "templates", "frankenplot", "skeleton", "skeleton.Rmd",
    package = "Jumble"
  )
  if (template_path == "" || !file.exists(template_path)) {
    # Fallback for development
    template_path <- file.path("inst", "rmarkdown", "templates",
                                "frankenplot", "skeleton", "skeleton.Rmd")
  }
  if (!file.exists(template_path)) {
    stop("Frankenplot Rmd template not found. Is the Jumble package installed correctly?")
  }

  # Derive sample name from tumor CSV filename
  sample_name <- basename(tumor_jumble_csv)
  sample_name <- sub("\\.jumble\\.csv$", "", sample_name)
  sample_name <- sub("\\.csv$", "", sample_name)

  # Ensure output directory exists
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Render
  rmarkdown::render(
    input = template_path,
    output_file = basename(output_file),
    output_dir = output_dir,
    params = list(
      bins_t = bins_t,
      segments_t = segments_t,
      bins_n = bins_n,
      segments_n = segments_n,
      alf = alf,
      alf_n = alf_n,
      salf = salf,
      galf_t = galf_t,
      galf_n = galf_n,
      hrdtable = hrdtable,
      qc_metrics = qc_metrics,
      dpyd_result = dpyd_result,
      dpyd_table = dpyd_table,
      genome = genome,
      sample_name = sample_name
    ),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )
}
