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
  # Increase timeout to 10 minutes
  if (requireNamespace("httr", quietly = TRUE)) {
    httr::set_config(httr::config(timeout = 3600))
  }
  # 1. Select Dataset --------------------------------------------------------
  if (genome == "hg19") {
    mart <- useEnsembl(
      biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
      GRCh = 37, mirror = mirror
    )
  } else if (genome == "hg38") {
    mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = mirror)
  } else {
    stop("Unsupported genome: ", genome)
  }

  chromosomes <- c(as.character(1:22), "X", "Y")

  # 2. Fetch Genes -----------------------------------------------------------
  # Note: 'external_gene_name' is the current attribute for symbol, formerly
  # 'hgnc_symbol' or 'gene_name'
  # 'ensembl_gene_id' is stable ID
  attributes_genes <- c(
    "ensembl_gene_id", "external_gene_name",
    "chromosome_name", "start_position", "end_position", "gene_biotype"
  )

  message("Fetching gene data from Ensembl...")
  if (genome == "hg19") {
      # Bulk fetch for hg19 (historically works fine and is faster)
      tryCatch({
          allgenes_df <- getBM(attributes = attributes_genes, mart = mart)
          allgenes <- as.data.table(allgenes_df)
      }, error = function(e) {
          warning("Bulk fetch failed, falling back to iterative: ", e$message)
          # Fallback logic could be added here, but for now just fail/warn
          allgenes <- data.table() 
      })
  } else {
       # Iterative fetch for hg38 (prone to timeouts)
       message("Iterating chromosomes for stability...")
       allgenes_list <- list()
       chromosomes <- c(as.character(1:22), "X", "Y")
       for (chrom in chromosomes) {
         message("  Fetching genes for chr ", chrom, "...")
         tryCatch({
            genes_chr <- getBM(attributes = attributes_genes, 
                               filters = "chromosome_name", 
                               values = chrom, 
                               mart = mart)
            if (nrow(genes_chr) > 0) {
                allgenes_list[[chrom]] <- as.data.table(genes_chr)
            }
         }, error = function(e) warning("Failed to fetch genes for chr ", chrom))
       }
       allgenes <- rbindlist(allgenes_list, fill = TRUE)
  }
  
  # Filter for standard chromosomes and protein coding
  # (Already filtered by fetching specific chromosomes, but keeping safety check)
  allgenes <- allgenes[chromosome_name %in% chromosomes]
  allgenes <- allgenes[gene_biotype == "protein_coding"]

  # Rename columns to match Jumble's internal naming convention
  # Jumble expects specific column names for downstream processing.
  # We map Ensembl attributes to these expected names:
  # 'ensembl_gene_id' -> 'Gene stable ID'
  # 'external_gene_name' -> 'Gene name'
  # 'chromosome_name' -> 'Chromosome/scaffold name'
  # 'start_position' -> 'Gene start (bp)'
  # 'end_position' -> 'Gene end (bp)'
  setnames(allgenes,
    old = c(
      "ensembl_gene_id", "external_gene_name", "chromosome_name",
      "start_position", "end_position"
    ),
    new = c(
      "Gene stable ID", "Gene name", "Chromosome/scaffold name",
      "Gene start (bp)", "Gene end (bp)"
    )
  )

  # 3. Fetch Exons -----------------------------------------------------------
  attributes_exons <- c(
    "ensembl_gene_id", "external_gene_name",
    "chromosome_name", "exon_chrom_start", "exon_chrom_end", "rank"
  )

  message("Fetching exon data from Ensembl...")
  if (genome == "hg19") {
       tryCatch({
          allexons_df <- getBM(attributes = attributes_exons, mart = mart)
          allexons <- as.data.table(allexons_df)
       }, error = function(e) {
          warning("Bulk fetch failed: ", e$message)
          allexons <- data.table()
       })
  } else {
      # Iterative fetch for hg38
      message("Iterating chromosomes for stability...")
      allexons_list <- list()
      chromosomes <- c(as.character(1:22), "X", "Y") # Re-define just in case
      for (chrom in chromosomes) {
        message("  Fetching exons for chr ", chrom, "...")
        tryCatch({
           exons_chr <- getBM(attributes = attributes_exons,
                              filters = "chromosome_name",
                              values = chrom,
                              mart = mart)
           if (nrow(exons_chr) > 0) {
               allexons_list[[chrom]] <- as.data.table(exons_chr)
           }
        }, error = function(e) warning("Failed to fetch exons for chr ", chrom))
      }
      allexons <- rbindlist(allexons_list, fill = TRUE)
  }

  # Filter
  allexons <- allexons[chromosome_name %in% chromosomes]
  allexons <- allexons[ensembl_gene_id %in% allgenes$`Gene stable ID`]

  # Rename
  setnames(allexons,
    old = c(
      "ensembl_gene_id", "external_gene_name", "chromosome_name",
      "exon_chrom_start", "exon_chrom_end", "rank"
    ),
    new = c(
      "Gene stable ID", "Gene name", "Chromosome/scaffold name",
      "Exon region start (bp)", "Exon region end (bp)",
      "Exon rank in transcript"
    )
  )


  # 4. Process Cancer Genes --------------------------------------------------
  cgc_path <- system.file("extdata", "cancer_genes.csv",
    package = "Jumble"
  )
  if (cgc_path == "") {
    # Fallback for development/testing when package not installed
    cgc_path <- "inst/extdata/cancer_genes.csv"
  }

  if (file.exists(cgc_path)) {
    cancergenes <- fread(cgc_path)
    
    # Ensure expected columns
    if (!all(c("hugo_symbol", "alteration") %in% names(cancergenes))) {
         warning("Cancer genes CSV missing required columns. Using empty list.")
         cancergenes <- data.table(
             hugo_symbol = character(),
             alteration = character(),
             ANNOT = character(),
             chromosome = character(),
             start = integer(),
             end = integer()
         )
    } else {
        # Map Alteration to ANNOT and keep Alteration
        cancergenes[, ANNOT := "AMBI"]
        cancergenes[alteration == "amp", ANNOT := "ONCO"]
        cancergenes[alteration == "del", ANNOT := "TSG"]
    
        # Map Ensembl ID and Coordinates
        # Join with allgenes
        mapping <- unique(allgenes[, .(
          `Gene name`, `Gene stable ID`,
          `Chromosome/scaffold name`, `Gene start (bp)`, `Gene end (bp)`
        )])
        
        cancergenes <- merge(cancergenes, mapping,
          by.x = "hugo_symbol",
          by.y = "Gene name", all.x = FALSE
        )
    
        setnames(cancergenes,
          old = c(
            "Gene stable ID", "Chromosome/scaffold name", "Gene start (bp)",
            "Gene end (bp)"
          ),
          new = c("ensembl_gene_id_version", "chromosome", "start", "end")
        )
    
        # Select columns - keep coordinates and alteration!
        cancergenes <- cancergenes[, .(
          hugo_symbol, ensembl_gene_id_version, ANNOT, alteration,
          chromosome, start, end
        )]
    }
  } else {
    warning("Cancer genes file (cancer_genes.csv) not found. Returning empty cancer genes table.")
    cancergenes <- data.table(
      hugo_symbol = character(), ensembl_gene_id_version = character(),
      ANNOT = character(), chromosome = character(), start = integer(),
      end = integer()
    )
  }

  return(list(
    cancergenes_clinseq = cancergenes,
    allgenes = allgenes,
    allexons = allexons
  ))
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
build_reference <- function(count_files, annotation_source = "biomart",
                            genome = NULL, output_file = NULL, cores = 1, mirror = NULL) {
  if (length(count_files) == 0) stop("No count files provided.")
  
  # 1. Read Count Files ------------------------------------------------------
  message("Reading count files...")
  allcounts <- list()
  ntargets <- numeric(length(count_files))
  bed_files <- character(length(count_files))
  detected_genomes <- character(length(count_files))
  
  for (i in seq_along(count_files)) {
    counts <- readRDS(count_files[i])
    if (is.null(counts$input_bam_file)) counts$input_bam_file <- count_files[i]
    
    ntargets[i] <- length(counts$count)
    
    bed <- counts$target_bed_file
    if (!is.null(bed)) {
      bed_files[i] <- basename(bed)
    } else {
      bed_files[i] <- "wgs"
    }
    
    # Collect genome if present
    if (!is.null(counts$genome)) {
      detected_genomes[i] <- counts$genome
    } else {
      detected_genomes[i] <- NA_character_
    }
    
    allcounts[[i]] <- counts
  }
  
  # 2. Determine Genome ------------------------------------------------------
  # 1. Argument takes precedence
  # 2. Else check counts files
  # 3. Default to hg19
  
  if (is.null(genome)) {
    unique_genomes <- unique(na.omit(detected_genomes))
    
    if (length(unique_genomes) == 1) {
      genome <- unique_genomes
      message("Detected genome from count files: ", genome)
    } else if (length(unique_genomes) > 1) {
      stop("Count files contain different reference genomes: ", 
           paste(unique_genomes, collapse = ", "))
    } else {
      genome <- "hg19"
      message("Genome not detected in count files, defaulting to: ", genome)
    }
  } else {
    # Check if provided genome matches detected ones (warn only?)
    unique_genomes <- unique(na.omit(detected_genomes))
    if (length(unique_genomes) > 0 && !all(unique_genomes == genome)) {
      warning("Provided genome ('", genome, 
              "') differs from genome detected in count files ('", 
              paste(unique_genomes, collapse = ", "), "')")
    }
  }

  # 3. Load Annotation -------------------------------------------------------
  if (annotation_source == "biomart") {
    annot <- generate_gene_annotation(genome = genome, mirror = mirror)
  } else if (file.exists(annotation_source) && grepl("\\.RDS$", annotation_source, ignore.case = TRUE)) {
    message("Loading annotation from RDS: ", annotation_source)
    annot <- readRDS(annotation_source)
  } else if (dir.exists(annotation_source)) {
    # Legacy mode: read from folder
    annot <- list(
      cancergenes_clinseq = fread(file.path(
        annotation_source,
        "cancergenes.txt"
      )),
      allgenes = fread(file.path(annotation_source, "allgenes.txt")),
      allexons = fread(file.path(annotation_source, "allexons.txt"))
    )
  } else {
    stop("Invalid annotation source.")
  }

  if (length(unique(ntargets)) > 1) {
    stop("Number of bins differs between samples.")
  }
  if (length(unique(bed_files)) > 1) stop("BED file differs between samples.")

  wgs <- (bed_files[1] == "wgs")

  # 4. Make Targets Template -------------------------------------------------
  counts1 <- allcounts[[1]]
  targets <- as.data.table(counts1$ranges)

  # Add metadata columns
  targets[, `:=`(
    sample = "",
    is_target = TRUE,
    type = "target",
    is_tiled = FALSE,
    chromosome = as.character(seqnames),
    mid = round((end + start) / 2),
    gene = ""
  )]

  # Robust Sorting: Ensure bins are ordered by chromosome (1..22, X, Y) and start position
  # This fixes artifacts where bins might be out of order in input counts (e.g. chr21 bin inside chr1 block)
  
  # Create a factor for sorting standard chromosomes
  chrom_levels <- c(as.character(1:22), "X", "Y")
  # Use clean names for sorting factor
  targets[, clean_chr := clean_chrom_names(chromosome)]
  # Map non-standard to NA (they sort last)
  targets[, chr_factor := factor(clean_chr, levels = chrom_levels)]
  
  # Calculate ordering
  ordering <- order(targets$chr_factor, targets$start)
  
  # Apply ordering to targets
  targets <- targets[ordering]
  targets[, bin := 1:.N] # Re-assign sequential bin IDs
  
  # Cleanup temp columns
  targets[, c("clean_chr", "chr_factor") := NULL]

  # Apply ordering to all count files to maintain consistency
  # And to counts1 which is used below
  counts1$ranges <- counts1$ranges[ordering]
  
  for (i in seq_along(allcounts)) {
       # Sort all relevant count vectors in the object
       if (!is.null(allcounts[[i]]$count)) {
           allcounts[[i]]$count <- allcounts[[i]]$count[ordering]
       }
       if (!is.null(allcounts[[i]]$count_short)) {
           allcounts[[i]]$count_short <- allcounts[[i]]$count_short[ordering]
       }
       # Also ranges need to be sorted to match
       allcounts[[i]]$ranges <- allcounts[[i]]$ranges[ordering]
  }

  # Identify background bins
  if (wgs) {
      # For WGS, assume standard bins. Any deviant width (e.g. ends of chroms) is background/ignored.
      min_w <- min(targets$width)
      targets[width != min_w, is_target := FALSE]
  } else {
      # For targeted panels, widths are variable (exons).
      # However, we must distinguish large background bins (backbone) from targets.
      # generate_counts uses bg_minsize = 300,000 for background bins.
      # So we can safely assume anything larger than e.g. 100kb is background.
      targets[, is_target := TRUE]
      targets[width > 100000, is_target := FALSE]
  }
  
  targets[is_target == FALSE, type := "background"]


  # 5. GC and Mappability ----------------------------------------------------
  message("Calculating GC content and Mappability...")
  ucsc_ranges <- counts1$ranges
  seqlevelsStyle(ucsc_ranges) <- "UCSC"

  # Select Genome
  if (genome == "hg19") {
    bsgenome <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  } else if (genome == "hg38") {
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
      stop("Package 'BSgenome.Hsapiens.UCSC.hg38' is required for hg38 support. Install with: BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')")
    }
    bsgenome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
  } else {
    stop("Unsupported genome for GC/Map calculation: ", genome)
  }

  # GC Content Calculation
  # We calculate the GC content for each bin using the selected  # GC Content Calculation
  # This is crucial for normalization, as coverage often correlates with GC
  # content. Using Biostrings directly instead of Repitools.
  message("Calculating GC content.")
  targets[, gc := as.double(NA)]
  target_ranges <- ucsc_ranges
  seqs <- Biostrings::getSeq(bsgenome, target_ranges)
  gc_values <- Biostrings::letterFrequency(seqs, letters = "GC", as.prob = TRUE)[,1]
  targets[is_target %in% c(TRUE, FALSE)]$gc <- gc_values


  # 6. Create Object ---------------------------------------------------------
  reference <- list(
    target_bed_file = if (wgs) "wgs" else counts1$target_bed_file,
    chromlength = counts1$chromlength,
    ranges = counts1$ranges,
    flag = counts1$flag,
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
    # Default naming reasoning: derived from target_bed_file or generic WGS name
    name <- reference$target_bed_file
    if (is.null(name) || name == "wgs") {
      name <- "jumble.WGS"
    }
    # Remove path from name if present
    name <- stringr::str_remove(name, ".*/")
    # Also remove extension if present (e.g. .bed) to avoid doubling
    name <- stringr::str_remove(name, "\\.[^.]+$")
    
    output_file <- paste0(name, ".reference.RDS")
  }
  
  saveRDS(reference, output_file)
  message("Reference saved to ", output_file)

  invisible(reference)
}
