#' Segment Data
#'
#' Performs segmentation on the normalized log-ratios using PSCBS.
#'
#' @param targets Normalized targets data.table.
#' @param alpha Significance level for segmentation (default: 0.02, or 1e-5 for WGS).
#' @param cancergenes Optional data.table of cancer genes for annotation.
#' @param cancerexons Optional data.table of cancer exons for annotation.
#' @return A list containing the segment table and the updated targets with segment IDs.
#' @importFrom PSCBS segmentByCBS
#' @importFrom data.table as.data.table :=
#' @importFrom stats median
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @export
segment_data <- function(targets, alpha = NULL, cancergenes = NULL,
                         cancerexons = NULL) {
  # 1. Determine Alpha -------------------------------------------------------
  if (is.null(alpha)) alpha <- 0.02

  # 2. Run CBS ---------------------------------------------------------------
  # PSCBS requires numeric chromosomes and sorted data.
  # We temporarily convert X/Y to 23/24 to ensure correct sorting and
  # processing.
  targets[, chromosome := as.character(chromosome)]
  targets[chromosome == "X", chromosome := "23"]
  targets[chromosome == "Y", chromosome := "24"]
  targets[, chromosome := as.numeric(chromosome)]

  targets <- targets[order(chromosome, start)]

  fit <- segmentByCBS(
    y = targets$log2,
    chromosome = targets$chromosome,
    alpha = alpha,
    undo = 1,
    verbose = FALSE
  )

  segments <- as.data.table(fit$output)

  # Convert back to character X/Y
  segments[, chromosome := as.character(chromosome)]
  segments[chromosome == "23", chromosome := "X"]
  segments[chromosome == "24", chromosome := "Y"]

  # Restore targets chromosome to character X/Y
  targets[, chromosome := as.character(chromosome)]
  targets[chromosome == "23", chromosome := "X"]
  targets[chromosome == "24", chromosome := "Y"]

  # 3. Filter Segments -------------------------------------------------------
  segments <- segments[!is.na(chromosome) & !is.na(mean)]

  # Add genomic positions
  segments[, start_pos := targets$start[ceiling(start)]]
  segments[, end_pos := targets$end[floor(end)]]

  # 4. Assign Segment IDs ----------------------------------------------------
  targets[, segment := as.integer(NA)]

  # Initialize annotation columns
  segments[, genes := ""]
  segments[, relevance := ""]

  # Prepare GRanges for annotation if available
  if (!is.null(cancergenes) && !is.null(cancerexons)) {
    segranges <- makeGRangesFromDataFrame(segments[, .(chromosome,
      start = start_pos, end = end_pos
    )])
  }

  for (i in seq_len(nrow(segments))) {
    # Indices
    idx_start <- ceiling(segments$start[i])
    idx_end <- floor(segments$end[i])

    # Assign
    targets[idx_start:idx_end, segment := i]

    # Re-calculate mean (median) from data to be sure
    seg_mean <- median(targets[idx_start:idx_end]$log2, na.rm = TRUE)
    segments$mean[i] <- seg_mean

    # Store bin indices (overwriting start/end to match original script behavior)
    # The original script replaces the index-based start/end with the actual bin
    # numbers
    segments$start[i] <- targets$bin[idx_start]
    segments$end[i] <- targets$bin[idx_end]
  }

  # 5. Annotate Segments -----------------------------------------------------
  if (!is.null(cancergenes) && !is.null(cancerexons) && nrow(cancergenes) > 0) {
    # Check if cancergenes has the required coordinate columns
    # (older references might not have them)
    has_coords <- all(c("chromosome", "start", "end") %in% names(cancergenes))

    if (!has_coords) {
      warning("cancergenes does not have coordinate columns (chromosome, start, end). Skipping gene annotation. Please rebuild reference with latest version.")
    } else {
      # Ensure coordinates are numeric
      segments[, start_pos := as.numeric(start_pos)]
      segments[, end_pos := as.numeric(end_pos)]

      # Handle X/Y for GRanges (makeGRangesFromDataFrame expects standard names
      # or numeric)
      # Our segments have character chromosomes (1..22, X, Y)
      # cancergenes should also have compatible chromosomes.

      segranges <- makeGRangesFromDataFrame(
        segments[, .(chromosome,
          start = start_pos, end = end_pos
        )],
        seqnames.field = "chromosome",
        start.field = "start",
        end.field = "end",
        keep.extra.columns = FALSE
      )

      # Prepare cancergenes GRanges
      # cancergenes has: hugo_symbol, ensembl_gene_id_version, ANNOT,
      # chromosome, start, end
      generanges <- makeGRangesFromDataFrame(cancergenes,
        seqnames.field = "chromosome",
        start.field = "start",
        end.field = "end",
        keep.extra.columns = TRUE
      )

      # Prepare cancerexons GRanges
      # cancerexons has: Gene stable ID, Gene name, Chromosome/scaffold name,
      # Exon region start (bp), Exon region end (bp), Exon rank in transcript
      # We need to rename for makeGRanges
      exons_prep <- cancerexons[, .(
        chromosome = `Chromosome/scaffold name`,
        start = `Exon region start (bp)`,
        end = `Exon region end (bp)`,
        id = `Gene stable ID`
      )]
      exonranges <- makeGRangesFromDataFrame(exons_prep,
        seqnames.field = "chromosome",
        start.field = "start",
        end.field = "end",
        keep.extra.columns = TRUE
      )

      # Find overlaps
      gene_overlap <- as.data.table(findOverlaps(segranges, generanges))
      exon_overlap <- as.data.table(findOverlaps(segranges, exonranges))

      # Iterate segments to annotate
      for (i in seq_len(nrow(segments))) {
        # Genes
        gene_ix <- gene_overlap[queryHits == i]$subjectHits

        if (length(gene_ix) > 0) {
          # Get gene info
          # cancergenes columns: hugo_symbol, ensembl_gene_id_version, ANNOT
          genetable <- cancergenes[gene_ix, .(
            symbol = hugo_symbol,
            id = ensembl_gene_id_version, type = ANNOT
          )]
          genetable[, exonic := FALSE]
          genetable[, label := ""]

          # Basic label: comma separated symbols
          label_str <- paste(genetable$symbol, collapse = ",")
          segments[i, genes := label_str]

          # Exons
          exon_ix <- exon_overlap[queryHits == i]$subjectHits

          if (length(exon_ix) > 0) {
            exontable <- exons_prep[exon_ix]
            # Mark genes as exonic if their ID matches any exon ID
            genetable[id %in% exontable$id, exonic := TRUE]
          }

          # Relevance label
          genetable[exonic == TRUE, label := symbol]
          genetable[
            exonic == TRUE & type != "AMBI",
            label := paste0(symbol, "|", type)
          ]

          rel_label <- paste(genetable[exonic == TRUE]$label, collapse = ",")
          segments[i, relevance := rel_label]
        }
      }
    }
  }
  return(list(segments = segments, targets = targets))
}
