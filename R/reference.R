#' Get Ensembl Mart Setup
#' @keywords internal
get_ensembl_mart <- function(genome, mirror = NULL) {
  if (requireNamespace("httr", quietly = TRUE)) {
    httr::set_config(httr::config(timeout = 3600))
  }
  if (genome == "hg19") {
    biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37, mirror = mirror)
  } else if (genome == "hg38") {
    biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = mirror)
  } else {
    stop("Unsupported genome: ", genome)
  }
}

#' Fetch Ensembl Data (Genes or Exons)
#' @keywords internal
fetch_ensembl_data <- function(mart, genome, attributes, type = "genes") {
  chromosomes <- c(as.character(1:22), "X", "Y")
  
  if (genome == "hg19") {
    tryCatch({
      df <- biomaRt::getBM(attributes = attributes, mart = mart)
      as.data.table(df)
    }, error = function(e) {
      warning("Bulk fetch failed for ", type, ": ", e$message)
      data.table::data.table()
    })
  } else {
    message("Iterating chromosomes for stability...")
    res_list <- list()
    for (chrom in chromosomes) {
      tryCatch({
        chr_data <- biomaRt::getBM(attributes = attributes, filters = "chromosome_name", values = chrom, mart = mart)
        if (nrow(chr_data) > 0) res_list[[chrom]] <- as.data.table(chr_data)
      }, error = function(e) warning("Failed to fetch ", type, " for chr ", chrom))
    }
    data.table::rbindlist(res_list, fill = TRUE)
  }
}

#' Process Local Cancer Genes
#' @keywords internal
process_cancer_genes <- function(allgenes) {
  cgc_path <- system.file("extdata", "cancer_genes.csv", package = "Jumble")
  if (cgc_path == "") cgc_path <- "inst/extdata/cancer_genes.csv"
  
  if (!file.exists(cgc_path)) {
    warning("cancer_genes.csv not found. Returning empty table.")
    return(data.table(hugo_symbol = character(), ensembl_gene_id_version = character(), ANNOT = character(), chromosome = character(), start = integer(), end = integer()))
  }
  
  cgenes <- fread(cgc_path)
  if (!all(c("hugo_symbol", "alteration") %in% names(cgenes))) return(data.table())
  
  cgenes[, ANNOT := "AMBI"]
  cgenes[alteration == "amp", ANNOT := "ONCO"]
  cgenes[alteration == "del", ANNOT := "TSG"]
  
  mapping <- unique(allgenes[, .(`Gene name`, `Gene stable ID`, `Chromosome/scaffold name`, `Gene start (bp)`, `Gene end (bp)`)])
  cgenes <- merge(cgenes, mapping, by.x = "hugo_symbol", by.y = "Gene name", all.x = FALSE)
  
  setnames(cgenes, 
           old = c("Gene stable ID", "Chromosome/scaffold name", "Gene start (bp)", "Gene end (bp)"),
           new = c("ensembl_gene_id_version", "chromosome", "start", "end"))
  
  cgenes[, .(hugo_symbol, ensembl_gene_id_version, ANNOT, alteration, chromosome, start, end)]
}

#' Generate Gene Annotation
#'
#' Fetches gene and exon information from Ensembl via biomaRt and integrates with
#' a bundled CSV file of clinically relevant cancer genes.
#'
#' @param genome Genome version ("hg19" or "hg38").
#' @return A list containing `cancergenes`, `allgenes`, and `allexons`.
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom data.table fread as.data.table setnames :=
#' @importFrom data.table fread as.data.table setnames :=
#' @importFrom utils data
#' @importFrom httr set_config config timeout
#' @export
generate_gene_annotation <- function(genome = "hg19", mirror = NULL) {
  mart <- get_ensembl_mart(genome, mirror)
  chromosomes <- c(as.character(1:22), "X", "Y")
  
  # Fetch and format genes
  attr_genes <- c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "gene_biotype")
  message("Fetching gene data from Ensembl...")
  allgenes <- fetch_ensembl_data(mart, genome, attr_genes, "genes")
  allgenes <- allgenes[chromosome_name %in% chromosomes & gene_biotype == "protein_coding"]
  setnames(allgenes, old = attr_genes[1:5], new = c("Gene stable ID", "Gene name", "Chromosome/scaffold name", "Gene start (bp)", "Gene end (bp)"))
  
  # Fetch and format exons
  attr_exons <- c("ensembl_gene_id", "external_gene_name", "chromosome_name", "exon_chrom_start", "exon_chrom_end", "rank")
  message("Fetching exon data from Ensembl...")
  allexons <- fetch_ensembl_data(mart, genome, attr_exons, "exons")
  allexons <- allexons[chromosome_name %in% chromosomes & ensembl_gene_id %in% allgenes$`Gene stable ID`]
  setnames(allexons, old = attr_exons, new = c("Gene stable ID", "Gene name", "Chromosome/scaffold name", "Exon region start (bp)", "Exon region end (bp)", "Exon rank in transcript"))
  
  # Process cancer genes
  cancergenes <- process_cancer_genes(allgenes)
  
  list(cancergenes_clinseq = cancergenes, allgenes = allgenes, allexons = allexons)
}


#' Load Count Files
#' @keywords internal
load_count_files <- function(count_files) {
  allcounts <- list()
  bed_files <- character(length(count_files))
  detected_genomes <- character(length(count_files))
  
  for (i in seq_along(count_files)) {
    counts <- readRDS(count_files[i])
    if (is.null(counts$input_bam_file)) counts$input_bam_file <- count_files[i]
    bed_files[i] <- if (!is.null(counts$target_bed_file)) basename(counts$target_bed_file) else "wgs"
    detected_genomes[i] <- if (!is.null(counts$genome)) counts$genome else NA_character_
    allcounts[[i]] <- counts
  }
  list(allcounts = allcounts, bed_files = bed_files, detected_genomes = detected_genomes)
}

#' Create Target Template & Order Data
#' @keywords internal
create_target_template <- function(counts1, is_wgs, allcounts) {
  targets <- as.data.table(counts1$ranges)
  targets[, `:=`(sample = "", is_target = TRUE, type = "target", is_tiled = FALSE, 
                 chromosome = as.character(seqnames), mid = round((end + start) / 2), gene = "")]
  
  # Robust Sorting
  chrom_levels <- c(as.character(1:22), "X", "Y")
  targets[, clean_chr := gsub("chr", "", chromosome)] # Example simple clean
  targets[, chr_factor := factor(clean_chr, levels = chrom_levels)]
  
  ordering <- order(targets$chr_factor, targets$start)
  targets <- targets[ordering]
  targets[, bin := 1:.N]
  targets[, c("clean_chr", "chr_factor") := NULL]
  
  # Define Background
  if (is_wgs) {
    min_w <- min(targets$width)
    targets[width != min_w, is_target := FALSE]
  } else {
    targets[width > 100000, is_target := FALSE]
  }
  targets[is_target == FALSE, type := "background"]
  
  # Apply ordering to all counts
  for (i in seq_along(allcounts)) {
    if (!is.null(allcounts[[i]]$count)) allcounts[[i]]$count <- allcounts[[i]]$count[ordering]
    if (!is.null(allcounts[[i]]$count_short)) allcounts[[i]]$count_short <- allcounts[[i]]$count_short[ordering]
    allcounts[[i]]$ranges <- allcounts[[i]]$ranges[ordering]
  }
  
  list(targets = targets, allcounts = allcounts)
}

#' Calculate GC Content
#' @keywords internal
calculate_gc_content <- function(targets, ucsc_ranges, genome) {
  GenomeInfoDb::seqlevelsStyle(ucsc_ranges) <- "UCSC"
  
  if (genome == "hg19") {
    bsgenome <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  } else if (genome == "hg38") {
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) stop("Install BSgenome.Hsapiens.UCSC.hg38")
    bsgenome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
  } else {
    stop("Unsupported genome for GC: ", genome)
  }
  
  seqs <- Biostrings::getSeq(bsgenome, ucsc_ranges)
  gc_values <- Biostrings::letterFrequency(seqs, letters = "GC", as.prob = TRUE)[,1]
  targets[is_target %in% c(TRUE, FALSE), gc := gc_values]
  return(targets)
}

#' Build Reference Panel
#'
#' Creates a reference object from a set of count files.
#'
#' @param count_files Vector of paths to count RDS files.
#' @param annotation_source Source of annotation ("biomart" or path to folder).
#' @param genome Genome version ("hg19" or "hg38"). If NULL, detected from
#'   count files.
#' @param output_file Path to save the output RDS file. If NULL, defaults to '<target>.reference.RDS' or 'jumble.WGS.reference.RDS'.
#' @param cores Number of cores for parallel processing.
#' @return The reference object (invisibly if saved).
#' @importFrom data.table fread data.table rbindlist as.data.table := setkey
#' @importFrom GenomicRanges makeGRangesFromDataFrame seqnames width start end
#' @importFrom GenomeInfoDb seqlevelsStyle<-
#' @importFrom Biostrings getSeq letterFrequency
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 Hsapiens
#' @importFrom stringr str_remove
#' @export
build_reference <- function(count_files, annotation_source = "biomart", genome = NULL, output_file = NULL, cores = 1, mirror = NULL) {
  if (length(count_files) == 0) stop("No count files provided.")
  
  # 1. Load Data
  message("Reading count files...")
  cf_data <- load_count_files(count_files)
  
  # 2. Determine Genome
  unique_genomes <- unique(na.omit(cf_data$detected_genomes))
  if (is.null(genome)) {
    genome <- if (length(unique_genomes) == 1) unique_genomes else "hg19"
  }
  
  # 3. Load Annotation
  if (annotation_source == "biomart") {
    annot <- generate_gene_annotation(genome = genome, mirror = mirror)
  } else if (file.exists(annotation_source) && grepl("\\.RDS$", annotation_source, ignore.case = TRUE)) {
    annot <- readRDS(annotation_source)
  } else {
    stop("Invalid annotation source.")
  }

  is_wgs <- (cf_data$bed_files[1] == "wgs")
  
  # 4. Create Template
  template_data <- create_target_template(cf_data$allcounts[[1]], is_wgs, cf_data$allcounts)
  targets <- template_data$targets
  allcounts <- template_data$allcounts
  
  # 5. Calculate GC
  message("Calculating GC content...")
  targets <- calculate_gc_content(targets, allcounts[[1]]$ranges, genome)
  
  # 6. Build and Save Object
  reference <- list(
    target_bed_file = if (is_wgs) "wgs" else allcounts[[1]]$target_bed_file,
    chromlength = allcounts[[1]]$chromlength,
    ranges = allcounts[[1]]$ranges,
    flag = allcounts[[1]]$flag,
    date = date(),
    samples = count_files,
    target_template = targets,
    allcounts = allcounts,
    allgenes = annot$allgenes,
    allexons = annot$allexons,
    cancergenes_clinseq = annot$cancergenes_clinseq,
    genome = genome
  )
  
  if (is.null(output_file)) {
    name <- stringr::str_remove(stringr::str_remove(reference$target_bed_file %||% "jumble.WGS", ".*/"), "\\.[^.]+$")
    output_file <- paste0(name, ".reference.RDS")
  }
  
  saveRDS(reference, output_file)
  message("Reference saved to ", output_file)
  invisible(reference)
}