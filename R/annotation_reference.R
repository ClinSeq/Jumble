# Reference annotation compatibility helpers.
#
# Internal helpers for handling legacy references with embedded annotation and
# newer slim references where annotation may be absent or resolved elsewhere.

valid_annotation_modes <- function() {
  c("auto", "none", "embedded", "cached", "bioc", "user")
}

normalize_annotation_mode <- function(annotation_mode = "auto") {
  if (length(annotation_mode) != 1 || is.na(annotation_mode)) {
    stop("annotation_mode must be a single value.", call. = FALSE)
  }
  annotation_mode <- as.character(annotation_mode)
  if (!annotation_mode %in% valid_annotation_modes()) {
    stop(
      "annotation_mode must be one of: ",
      paste(valid_annotation_modes(), collapse = ", "),
      call. = FALSE
    )
  }
  annotation_mode
}

is_nonempty_table <- function(x) {
  !is.null(x) && is.data.frame(x) && nrow(x) > 0
}

empty_reference_annotations <- function(mode = "none") {
  list(
    mode = mode,
    allgenes = NULL,
    allexons = NULL,
    cancergenes_clinseq = NULL,
    cytobands = NULL,
    has_gene_annotation = FALSE,
    has_exon_annotation = FALSE,
    has_cancer_gene_annotation = FALSE,
    has_cytoband_annotation = FALSE
  )
}

reference_has_embedded_annotation <- function(reference) {
  any(vapply(
    c("allgenes", "allexons", "cancergenes_clinseq", "cytobands"),
    function(field) is_nonempty_table(reference[[field]]),
    logical(1)
  ))
}

reference_schema_label <- function(reference) {
  if (!is.null(reference$reference_schema_version)) {
    return(as.character(reference$reference_schema_version))
  }
  if (reference_has_embedded_annotation(reference)) {
    return("legacy-embedded")
  }
  "slim"
}

resolve_reference_annotations <- function(reference, annotation_mode = "auto") {
  annotation_mode <- normalize_annotation_mode(annotation_mode)

  if (annotation_mode == "none") {
    return(empty_reference_annotations(mode = "none"))
  }

  if (annotation_mode %in% c("cached", "bioc", "user")) {
    stop(
      "annotation_mode='", annotation_mode,
      "' is reserved for the external annotation resolver and is not implemented yet. ",
      "Use annotation_mode='auto', 'embedded', or 'none'.",
      call. = FALSE
    )
  }

  has_embedded <- reference_has_embedded_annotation(reference)

  if (annotation_mode == "embedded" && !has_embedded) {
    stop(
      "annotation_mode='embedded' was requested, but the reference contains no embedded annotation resources.",
      call. = FALSE
    )
  }

  if (!has_embedded) {
    return(empty_reference_annotations(mode = "none"))
  }

  annotations <- empty_reference_annotations(mode = "embedded")
  annotations$allgenes <- if (is_nonempty_table(reference$allgenes)) reference$allgenes else NULL
  annotations$allexons <- if (is_nonempty_table(reference$allexons)) reference$allexons else NULL
  annotations$cancergenes_clinseq <- if (is_nonempty_table(reference$cancergenes_clinseq)) reference$cancergenes_clinseq else NULL
  annotations$cytobands <- if (is_nonempty_table(reference$cytobands)) reference$cytobands else NULL
  annotations$has_gene_annotation <- !is.null(annotations$allgenes)
  annotations$has_exon_annotation <- !is.null(annotations$allexons)
  annotations$has_cancer_gene_annotation <- !is.null(annotations$cancergenes_clinseq)
  annotations$has_cytoband_annotation <- !is.null(annotations$cytobands)
  annotations
}

normalize_reference_object <- function(reference, annotation_mode = "auto") {
  if (is.null(reference) || !is.list(reference)) {
    stop("reference must be a reference list or path to a reference RDS file.", call. = FALSE)
  }

  annotation_mode <- normalize_annotation_mode(annotation_mode)
  reference$reference_schema_version <- reference_schema_label(reference)
  reference$annotation_mode <- annotation_mode
  reference$annotation_resources <- resolve_reference_annotations(reference, annotation_mode)
  reference
}

reference_annotations <- function(reference) {
  if (!is.null(reference$annotation_resources)) {
    return(reference$annotation_resources)
  }
  resolve_reference_annotations(reference, annotation_mode = "auto")
}

reference_annotation_table <- function(reference, field) {
  annotations <- reference_annotations(reference)
  value <- annotations[[field]]
  if (is_nonempty_table(value)) value else NULL
}
