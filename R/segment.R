#' Run Circular Binary Segmentation (CBS)
#'
#' Prepares the target data by standardizing chromosome names, sorting, 
#' and running the PSCBS segmentation algorithm.
#'
#' @param targets Normalized targets data.table.
#' @param alpha Significance level for segmentation.
#' @return A list containing the `segments` data.table and the sorted `targets` data.table.
#' @importFrom PSCBS segmentByCBS
#' @importFrom data.table as.data.table :=
#' @keywords internal
run_cbs <- function(targets, alpha) {
  # Standardize and sort for PSCBS
  targets[, chromosome := as.character(chromosome)]
  targets[chromosome == "X", chromosome := "23"]
  targets[chromosome == "Y", chromosome := "24"]
  targets[, chromosome := as.numeric(chromosome)]
  
  # Note: This ordering creates a new data.table, breaking reference semantics
  targets <- targets[order(chromosome, start)]
  
  fit <- segmentByCBS(
    y = targets$log2,
    chromosome = targets$chromosome,
    alpha = alpha,
    undo = 1,
    verbose = FALSE
  )
  
  segments <- as.data.table(fit$output)
  
  # Convert numeric representations back to character X/Y
  segments[, chromosome := as.character(chromosome)]
  segments[chromosome == "23", chromosome := "X"]
  segments[chromosome == "24", chromosome := "Y"]
  
  targets[, chromosome := as.character(chromosome)]
  targets[chromosome == "23", chromosome := "X"]
  targets[chromosome == "24", chromosome := "Y"]
  
  return(list(segments = segments, targets = targets))
}

#' Process and Map Segments to Targets
#'
#' Filters out invalid segments, assigns segment IDs to the target bins, 
#' and recalculates the median log-ratio for each segment.
#'
#' @param segments Segmentation data.table from CBS.
#' @param targets Sorted targets data.table.
#' @return A list containing the updated `segments` and `targets` data.tables.
#' @importFrom data.table :=
#' @importFrom stats median
#' @keywords internal
process_segments <- function(segments, targets) {
  # Filter out NAs
  segments <- segments[!is.na(chromosome) & !is.na(mean)]
  
  # Add genomic positions based on target bins
  segments[, start_pos := targets$start[ceiling(start)]]
  segments[, end_pos := targets$end[floor(end)]]
  
  # Initialize new columns
  targets[, segment := as.integer(NA)]
  segments[, genes := ""]
  segments[, relevance := ""]
  
  for (i in seq_len(nrow(segments))) {
    idx_start <- ceiling(segments$start[i])
    idx_end <- floor(segments$end[i])
    
    # Map segment ID back to target bins
    targets[idx_start:idx_end, segment := i]
    
    # Re-calculate mean (median) from data to be robust
    seg_mean <- median(targets[idx_start:idx_end]$log2, na.rm = TRUE)
    segments$mean[i] <- seg_mean
    
    # Store actual bin names/indices instead of row indices
    segments$start[i] <- targets$bin[idx_start]
    segments$end[i] <- targets$bin[idx_end]
  }
  
  return(list(segments = segments, targets = targets))
}

#' Annotate Segments with Cancer Genes and Exons
#'
#' Finds overlaps between derived genomic segments and known cancer genes/exons,
#' appending relevance labels to the segments.
#'
#' @param segments Processed segments data.table.
#' @param cancergenes Data.table of cancer genes with coordinates.
#' @param cancerexons Data.table of cancer exons with coordinates.
#' @return The annotated `segments` data.table.
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom data.table as.data.table :=
#' @keywords internal
annotate_segments <- function(segments, cancergenes, cancerexons) {
  has_coords <- all(c("chromosome", "start", "end") %in% names(cancergenes))
  
  if (!has_coords) {
    warning("cancergenes does not have coordinate columns (chromosome, start, end). Skipping gene annotation. Please rebuild reference with latest version.")
    return(segments)
  }
  
  segments[, start_pos := as.numeric(start_pos)]
  segments[, end_pos := as.numeric(end_pos)]
  
  segranges <- makeGRangesFromDataFrame(
    segments[, .(chromosome, start = start_pos, end = end_pos)],
    seqnames.field = "chromosome",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = FALSE
  )
  
  generanges <- makeGRangesFromDataFrame(
    cancergenes,
    seqnames.field = "chromosome",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )
  
  exons_prep <- cancerexons[, .(
    chromosome = `Chromosome/scaffold name`,
    start = `Exon region start (bp)`,
    end = `Exon region end (bp)`,
    id = `Gene stable ID`
  )]
  
  exonranges <- makeGRangesFromDataFrame(
    exons_prep,
    seqnames.field = "chromosome",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )
  
  gene_overlap <- as.data.table(findOverlaps(segranges, generanges))
  exon_overlap <- as.data.table(findOverlaps(segranges, exonranges))
  
  for (i in seq_len(nrow(segments))) {
    gene_ix <- gene_overlap[queryHits == i]$subjectHits
    
    if (length(gene_ix) > 0) {
      genetable <- cancergenes[gene_ix, .(
        symbol = hugo_symbol,
        id = ensembl_gene_id_version, 
        type = ANNOT
      )]
      genetable[, exonic := FALSE]
      genetable[, label := ""]
      
      segments[i, genes := paste(genetable$symbol, collapse = ",")]
      
      exon_ix <- exon_overlap[queryHits == i]$subjectHits
      
      if (length(exon_ix) > 0) {
        exontable <- exons_prep[exon_ix]
        genetable[id %in% exontable$id, exonic := TRUE]
      }
      
      genetable[exonic == TRUE, label := symbol]
      genetable[exonic == TRUE & type != "AMBI", label := paste0(symbol, "|", type)]
      
      rel_label <- paste(genetable[exonic == TRUE]$label, collapse = ",")
      segments[i, relevance := rel_label]
    }
  }
  
  return(segments)
}

#' Segment Data
#'
#' Main wrapper function that performs segmentation on normalized log-ratios using PSCBS,
#' assigns mapped segments to targets, and optionally annotates the data.
#'
#' @param targets Normalized targets data.table.
#' @param alpha Significance level for segmentation (default: 0.02, or 1e-5 for WGS).
#' @param cancergenes Optional data.table of cancer genes for annotation.
#' @param cancerexons Optional data.table of cancer exons for annotation.
#' @return A list containing the segment table and the updated targets with segment IDs.
#' @export
segment_data <- function(targets, alpha = NULL, cancergenes = NULL, cancerexons = NULL) {
  if (is.null(alpha)) alpha <- 0.02
  
  # 1. Run CBS Segmentation
  cbs_results <- run_cbs(targets, alpha)
  
  # 2. Map and Process Segments
  processed <- process_segments(cbs_results$segments, cbs_results$targets)
  segments <- processed$segments
  targets <- processed$targets
  
  # 3. Annotate Segments (if references are provided)
  if (!is.null(cancergenes) && !is.null(cancerexons) && nrow(cancergenes) > 0) {
    segments <- annotate_segments(segments, cancergenes, cancerexons)
  }
  
  return(list(segments = segments, targets = targets))
}