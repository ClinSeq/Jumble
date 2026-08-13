#' Annotate Targets for Plotting
#'
#' Adds 'label' and 'selected_genes' columns to targets.
#'
#' @param targets data.table.
#' @return Modifies targets in-place.
#' @keywords internal
annotate_targets <- function(targets) {
  if (!"label" %in% names(targets)) {
    if ("type" %in% names(targets)) {
      targets[, label := type]
      targets[label == "tiled", label := "dense"]
      targets[label == "target", label := "sparse"]
    } else {
      targets[, label := "sparse"] # Fallback
    }
  }

  # Gene Coloring
  label_genes <- c(
    "TP53", "BRCA1", "BRCA2", "EGFR", "KRAS", "PIK3CA", "PTEN", "APC",
    "MYC", "BRAF", "ERBB2", "RB1", "MET", "ALK", "CDH1", "ATM",
    "CDKN2A", "VHL", "FLT3", "IDH1"
  )

  if (!"selected_genes" %in% names(targets)) {
    targets[, selected_genes := as.character(NA)]

    for (g in label_genes) {
      # Match exact: ^g$ OR ,g$ OR ^g, OR ,g,
      targets[
        type != "background" & (
          grepl(paste0("^", g, "$"), gene) |
            grepl(paste0(",", g, "$"), gene) |
            grepl(paste0("^", g, ","), gene) |
            grepl(paste0(",", g, ","), gene)
        ),
        selected_genes := g
      ]
    }

    # Factorize selected_genes
    found_genes <- sort(unique(targets$selected_genes[!is.na(targets$selected_genes)]))
    targets[, selected_genes := factor(selected_genes, levels = found_genes)]
  }
}

#' Setup Plot Theme and Colors
#'
#' @return List with theme settings and color values
#' @keywords internal
setup_plot_theme <- function() {
  ggplot2::theme_set(ggplot2::theme_bw())
  
  list(
    colorvalues = c(
      "background" = "#1010D0", "dense" = "#F8566D",
      "exonic" = "#00BFA4", "sparse" = "#000000", "bin" = "#000000"
    ),
    size = 1,
    size_selected = 1.1,
    pointcolor = "#000000",
    alpha = 0.3
  )
}

#' Prepare Y-Axis Parameters
#'
#' Calculate y-axis limits and breaks for log2 plots
#'
#' @param targets data.table with smooth_log2 column
#' @return List with ylims, ybreaks, yminorbreaks
#' @keywords internal
prepare_yaxis_params <- function(targets) {
  ylims <- c(
    min(0.4, min(2^targets[chromosome != "Y"]$smooth_log2, na.rm = TRUE)),
    max(3, max(2^targets$smooth_log2, na.rm = TRUE))
  )
  
  ybreaks <- c(.5, .75, 1, 1.5, 2, 3, 4, 6, 8)
  ybreaks <- ybreaks[ybreaks >= ylims[1] & ybreaks <= ylims[2]]
  yminorbreaks <- c(1.25, 1.75)
  
  list(ylims = ylims, ybreaks = ybreaks, yminorbreaks = yminorbreaks)
}

#' Prepare Depth (X-axis) Parameters
#'
#' Setup log-scale depth axis parameters
#'
#' @return List with limits, limits_labels, limits_breaks
#' @keywords internal
prepare_depth_axis <- function() {
  list(
    limits = c(1e1, 1e5),
    limits_labels = c("10", "30", "100", "300", "1k", "3k", "10k", "30k", "100k"),
    limits_breaks = c(1e1, 3e1, 1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5)
  )
}

#' Prepare Chromosome Data
#'
#' Create chromosome position mappings for plotting
#'
#' @param targets data.table with chromosome and end columns
#' @param reference Optional reference object with chromlength
#' @return List with chroms and chroms_order data.tables
#' @keywords internal
prepare_chromosome_data <- function(targets, reference = NULL) {
  chrom_levels <- c(as.character(1:22), "X", "Y")
  
  # Get chromosome lengths
  if (!is.null(reference) && !is.null(reference$chromlength)) {
    chroms <- data.table::data.table(
      chromosome = names(reference$chromlength),
      length = as.numeric(reference$chromlength)
    )
  } else {
    chroms <- targets[, .(length = max(end)), by = chromosome]
  }
  
  chroms[, chromosome := clean_chrom_names(as.character(chromosome))]
  chroms <- chroms[match(chrom_levels, chromosome)]
  chroms <- chroms[!is.na(chromosome)]
  
  # Calculate cumulative positions
  chroms[, start := as.double(0)]
  chroms[, stop := as.double(length)]
  chroms[, mid := as.double(round(length / 2))]
  
  if (nrow(chroms) >= 2) {
    for (i in 2:nrow(chroms)) {
      chroms[i, start := chroms$stop[i - 1]]
      chroms[i, stop := chroms$stop[i - 1] + length]
      chroms[i, mid := chroms$stop[i - 1] + round(length / 2)]
    }
  }
  
  # Chromosome positions by order
  chroms_order <- data.table::data.table(
    chromosome = unique(targets$chromosome),
    start = 0,
    end = 0,
    mid = 0
  )
  for (chr in chroms_order$chromosome) {
    chroms_order$start[chr == chroms_order$chromosome] <- targets[chromosome == chr, min(bin)]
    chroms_order$end[chr == chroms_order$chromosome] <- targets[chromosome == chr, max(bin)]
    chroms_order$mid[chr == chroms_order$chromosome] <- targets[chromosome == chr, mean(bin)]
  }
  
  chroms_order <- chroms_order[chromosome %in% chrom_levels]
  
  list(chroms = chroms, chroms_order = chroms_order)
}

#' Add Genomic Positions to Targets and Segments
#'
#' Calculate genomic coordinates across chromosomes
#'
#' @param targets data.table
#' @param segments data.table
#' @param chroms data.table with chromosome lengths
#' @return List with modified targets and segments
#' @keywords internal
add_genomic_positions <- function(targets, segments, chroms) {
  # Add genomic positions to targets
  targets[, gpos := as.numeric(start)]
  for (chr in unique(targets$chromosome)[-1]) {
    chr_idx <- which(chroms$chromosome == chr)
    if (length(chr_idx) > 0) {
      offset <- sum(chroms[1:(chr_idx - 1)]$length)
      targets[chromosome == chr, gpos := gpos + offset]
    }
  }
  
  # Add genomic positions to segments
  segments[, gstart := as.double(start_pos)]
  segments[, gstop := as.double(end_pos)]
  for (chr in unique(segments$chromosome)[-1]) {
    chr_idx <- which(chroms$chromosome == chr)
    if (length(chr_idx) > 0) {
      offset <- sum(chroms[1:(chr_idx - 1)]$length)
      segments[chromosome == chr, gstart := gstart + offset]
      segments[chromosome == chr, gstop := gstop + offset]
    }
  }
  
  list(targets = targets, segments = segments)
}

#' Prepare Cancer Gene Summary for Labeling
#'
#' Create one point per gene for plot annotations
#'
#' @param targets data.table with selected_genes
#' @return data.table with summarized cancer genes
#' @keywords internal
prepare_cancergenes_summary <- function(targets) {
  cancergenes_summary <- targets[!is.na(selected_genes), .(
    bin = median(bin, na.rm = TRUE),
    gpos = median(gpos, na.rm = TRUE),
    count = median(count, na.rm = TRUE),
    log2 = median(smooth_log2, na.rm = TRUE)
  ), by = selected_genes]
  
  if (nrow(cancergenes_summary) > 0) {
    # Calculate median depth
    m_depth <- median(targets[type != "background" & chromosome %in% 1:22]$count, na.rm = TRUE)
    if (is.na(m_depth)) m_depth <- 100
    
    cancergenes_summary[, nudge := 0.1]
    
    # Apply nudge logic
    if (nrow(cancergenes_summary) >= 2) {
      cancergenes_summary[seq(2, .N, 2), nudge := -1]
    }
    cancergenes_summary[count < m_depth * 0.8, nudge := 1]
    cancergenes_summary[count > m_depth * 1.2, nudge := -1]
  }
  
  cancergenes_summary
}

#' Prepare SNP Data for Plotting
#'
#' Filter and map SNP data to targets
#'
#' @param snp_table data.table with SNP data
#' @param targets data.table to map SNPs to
#' @return Modified targets with allele_ratio, DP, and maf columns
#' @keywords internal
prepare_snp_data <- function(snp_table, targets) {
  targets[, allele_ratio := as.double(NA)]
  targets[, DP := as.double(NA)]
  targets[, maf := as.double(NA)]

  if (is.null(snp_table) || nrow(snp_table) == 0) {
    return(targets)
  }

  # Prepare allele ratio column
  snp_table[, allele_ratio_use := allele_ratio]
  if (!is.null(snp_table$allele_ratio_corrected2)) {
    snp_table[, allele_ratio_use := allele_ratio_corrected2]
  }
  
  # Filter SNPs
  snp_table <- snp_table[type != "other"][allele_ratio_use < .99][
    allele_ratio_use > .01
  ][RD > 2][AD > 2]
  
  median_dp <- median(snp_table$DP, na.rm = TRUE)
  snp_table <- snp_table[DP > median_dp / 5][DP < median_dp * 10]
  
  # Map to targets
  idx <- match(targets$bin, snp_table$bin)
  targets[
    !is.na(idx),
    `:=`(
      allele_ratio = snp_table$allele_ratio_use[idx[!is.na(idx)]],
      DP = snp_table$DP[idx[!is.na(idx)]]
    )
  ]
  
  # Calculate MAF and smooth
  targets[, maf := abs(allele_ratio - .5) + .5]
  targets[!is.na(maf), maf := stats::runmed(maf, 9)]
  
  targets
}

#' Prepare Somatic Data for Plotting
#'
#' Map somatic variants to targets and create labels
#'
#' @param somatic_table data.table with somatic variants
#' @param targets data.table to map to
#' @param reference Optional reference object with cancergenes
#' @return List with somatic_plot and somatic_labels data.tables
#' @keywords internal
prepare_somatic_data <- function(somatic_table, targets, reference = NULL) {
  if (is.null(somatic_table) || nrow(somatic_table) == 0) {
    return(list(somatic_plot = NULL, somatic_labels = NULL))
  }
  
  somatic_plot <- data.table::copy(somatic_table)
  somatic_plot <- somatic_plot[!is.na(bin)]
  
  if (nrow(somatic_plot) == 0) {
    return(list(somatic_plot = NULL, somatic_labels = NULL))
  }
  
  # Map to targets
  cols_needed <- c("bin", "gpos", "smooth_log2", "log2", "count")
  cols_exist <- cols_needed %in% names(targets)
  
  if (!all(cols_exist)) {
    warning("Targets missing columns for somatic mapping: ", 
            paste(cols_needed[!cols_exist], collapse = ", "))
    return(list(somatic_plot = NULL, somatic_labels = NULL))
  }
  
  t_lookup <- targets[, cols_needed, with = FALSE]
  data.table::setkey(t_lookup, bin)
  data.table::setkey(somatic_plot, bin)
  somatic_plot <- t_lookup[somatic_plot]
  
  # Classify variant type for shape mapping
  if (all(c("REF", "ALT") %in% names(somatic_plot))) {
    somatic_plot[, var_type := data.table::fifelse(
      nchar(REF) > nchar(ALT), "Deletion",
      data.table::fifelse(nchar(REF) < nchar(ALT), "Insertion", "SNV")
    )]
    somatic_plot[, var_type := factor(var_type, levels = c("SNV", "Insertion", "Deletion"))]
  } else {
    somatic_plot[, var_type := factor("SNV", levels = c("SNV", "Insertion", "Deletion"))]
  }
  
  # Prepare labels
  somatic_labels <- NULL
  cancergenes_list <- NULL
  
  if (!is.null(reference) && !is.null(reference$cancergenes_clinseq)) {
    cancergenes_list <- reference$cancergenes_clinseq$gene
  }
  
  if (is.null(cancergenes_list) && !is.null(targets$selected_genes)) {
    cancergenes_list <- unique(na.omit(as.character(targets$selected_genes)))
  }
  
  if (nrow(somatic_plot) > 0 && !is.null(cancergenes_list)) {
    req_cols <- c("SYMBOL", "CANONICAL", "effect", "Consequence", 
                  "Protein_position", "Amino_acids")
    
    if (all(req_cols %in% names(somatic_plot))) {
      somatic_labels <- somatic_plot[
        SYMBOL %in% cancergenes_list & CANONICAL == "YES" & 
          effect %in% c("hotspot", "high-impact")
      ]
      
      if (nrow(somatic_labels) > 0) {
        somatic_labels[, ppos := stringr::str_extract(Protein_position, "[0-9]*")]
        somatic_labels[, aa := stringr::str_extract(Amino_acids, "[A-Z]$")]
        somatic_labels[is.na(aa), aa := ""]
        somatic_labels[Consequence == "stop_gained", aa := "*"]
        somatic_labels[grepl("frameshift", Consequence), aa := "fs"]
        somatic_labels[grepl("inframe", Consequence), aa := "if"]
        somatic_labels[grepl("splice", Consequence), aa := "sp"]
        
        somatic_labels[, label := paste0(SYMBOL, ":", ppos, aa)]
        somatic_labels[, nudge := 0.1]
        somatic_labels[AF > 0.8, nudge := -0.1]
      } else {
        somatic_labels <- NULL
      }
    }
  }
  
  list(somatic_plot = somatic_plot, somatic_labels = somatic_labels)
}

#' Create GC vs Raw Depth Plot
#'
#' @param targets data.table
#' @param theme_params List from setup_plot_theme()
#' @param axis_params List from prepare_depth_axis()
#' @return ggplot object
#' @keywords internal
plot_gc_rawdepth <- function(targets, theme_params, axis_params) {
  ggplot(targets[count > 0]) +
    xlab("Target GC content") +
    ylab("Fragments") +
    xlim(c(.16, .84)) +
    geom_point(
      mapping = aes(x = gc, y = count),
      col = theme_params$pointcolor, shape = 21, size = 1, 
      alpha = theme_params$alpha / 2
    ) +
    facet_wrap(facets = vars(label), nrow = 1) +
    theme(
      panel.spacing = unit(0, "lines"), 
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    geom_text(
      data = unique(targets[, .(label)])[label == "target", label := "sparse"],
      mapping = aes(label = label, x = 0.16, y = Inf), 
      hjust = 0, vjust = 1.5
    ) +
    scale_y_log10(
      limits = axis_params$limits, 
      breaks = axis_params$limits_breaks,
      minor_breaks = axis_params$limits_breaks,
      labels = axis_params$limits_labels
    )
}

#' Create Order vs Raw Depth Plot
#'
#' @param targets data.table
#' @param cancergenes_summary data.table
#' @param chroms_order data.table
#' @param theme_params List from setup_plot_theme()
#' @param axis_params List from prepare_depth_axis()
#' @return ggplot object
#' @keywords internal
plot_order_rawdepth <- function(targets, cancergenes_summary, chroms_order, 
                                theme_params, axis_params) {
  ggplot(targets) +
    xlab("Order of genomic position") +
    ylab("Fragments") +
    geom_vline(
      data = cancergenes_summary,
      mapping = aes(xintercept = bin, color = selected_genes),
      alpha = 0.5, linewidth = 0.5, show.legend = FALSE
    ) +
    geom_point(
      data = targets[!is.na(label) & is.na(selected_genes) & count > 0],
      mapping = aes(x = bin, y = count),
      fill = "#606060", col = "#202020", shape = 21, 
      size = theme_params$size, alpha = theme_params$alpha
    ) +
    geom_point(
      data = targets[!is.na(label) & !is.na(selected_genes) & count > 0],
      mapping = aes(x = bin, y = count, fill = selected_genes),
      col = "black", shape = 21, size = theme_params$size_selected, alpha = 0.8
    ) +
    geom_text(
      data = cancergenes_summary,
      mapping = aes(x = bin, y = Inf, label = selected_genes, color = selected_genes),
      angle = 45, hjust = 0, size = 2.5, vjust = 0, show.legend = FALSE
    ) +
    scale_fill_hue(l = 70) +
    scale_color_hue(l = 70) +
    scale_y_log10(
      limits = axis_params$limits,
      breaks = axis_params$limits_breaks,
      minor_breaks = axis_params$limits_breaks,
      labels = axis_params$limits_labels
    ) +
    scale_x_continuous(
      breaks = chroms_order$mid, minor_breaks = chroms_order$start[-1],
      expand = c(.01, .01), labels = chroms_order$chromosome
    ) +
    coord_cartesian(clip = "off") +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_line(),
      panel.grid.minor.x = element_line(color = "black"),
      axis.line = element_line(),
      axis.ticks = element_line(),
      plot.margin = unit(c(2, 0.5, 0.5, 0.5), "lines")
    )
}

#' Create Position vs Raw Depth Plot
#'
#' @param targets data.table
#' @param cancergenes_summary data.table
#' @param chroms data.table
#' @param theme_params List from setup_plot_theme()
#' @param axis_params List from prepare_depth_axis()
#' @return ggplot object
#' @keywords internal
plot_pos_rawdepth <- function(targets, cancergenes_summary, chroms, 
                              theme_params, axis_params) {
  ggplot(targets) +
    xlab("Genomic position") +
    ylab("Fragments") +
    geom_vline(
      data = cancergenes_summary,
      mapping = aes(xintercept = gpos, color = selected_genes),
      alpha = 0.5, linewidth = 0.5, show.legend = FALSE
    ) +
    geom_point(
      data = targets[!is.na(label) & is.na(selected_genes) & count > 0],
      mapping = aes(x = gpos, y = count),
      fill = "#606060", col = "#202020", shape = 21,
      size = theme_params$size, alpha = theme_params$alpha
    ) +
    geom_point(
      data = targets[!is.na(label) & !is.na(selected_genes) & count > 0],
      mapping = aes(x = gpos, y = count, fill = selected_genes),
      col = "black", shape = 21, size = theme_params$size_selected, alpha = 0.8
    ) +
    geom_text(
      data = cancergenes_summary,
      mapping = aes(x = gpos, y = Inf, label = selected_genes, color = selected_genes),
      angle = 45, hjust = 0, size = 2.5, vjust = 0, show.legend = FALSE
    ) +
    scale_fill_hue(l = 70) +
    scale_color_hue(l = 70) +
    scale_y_log10(
      limits = axis_params$limits,
      breaks = axis_params$limits_breaks,
      minor_breaks = axis_params$limits_breaks,
      labels = axis_params$limits_labels
    ) +
    scale_x_continuous(
      breaks = chroms$mid, minor_breaks = chroms$start[-1],
      expand = c(.01, .01), labels = chroms$chromosome
    ) +
    coord_cartesian(clip = "off") +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_line(),
      panel.grid.minor.x = element_line(color = "black"),
      axis.line = element_line(),
      axis.ticks = element_line(),
      plot.margin = unit(c(2, 0.5, 0.5, 0.5), "lines")
    )
}

#' Create GC vs Log2 Plot
#'
#' @param targets data.table
#' @param yaxis_params List from prepare_yaxis_params()
#' @param theme_params List from setup_plot_theme()
#' @return ggplot object
#' @keywords internal
plot_gc_log2 <- function(targets, yaxis_params, theme_params) {
  ggplot(targets) +
    xlab("Target GC content") +
    ylab("Corrected depth") +
    xlim(c(.16, .84)) +
    geom_point(
      mapping = aes(x = gc, y = 2^log2),
      col = theme_params$pointcolor, shape = 21, size = 1, 
      alpha = theme_params$alpha / 2
    ) +
    facet_wrap(facets = vars(label), nrow = 1) +
    theme(
      panel.spacing = unit(0, "lines"),
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    scale_y_log10(
      limits = yaxis_params$ylims, 
      breaks = yaxis_params$ybreaks, 
      minor_breaks = yaxis_params$yminorbreaks
    )
}

#' Create Order vs Log2 Plot
#'
#' @param targets data.table
#' @param segments data.table
#' @param cancergenes_summary data.table
#' @param chroms_order data.table
#' @param yaxis_params List from prepare_yaxis_params()
#' @param theme_params List from setup_plot_theme()
#' @return ggplot object
#' @keywords internal
plot_order_log2 <- function(targets, segments, cancergenes_summary, chroms_order,
                            yaxis_params, theme_params) {
  ggplot(targets) +
    xlab("Order of genomic position") +
    ylab("Corrected depth") +
    geom_vline(
      data = cancergenes_summary,
      mapping = aes(xintercept = bin, color = selected_genes),
      alpha = 0.5, linewidth = 0.5, show.legend = FALSE
    ) +
    geom_point(
      data = targets[!is.na(label) & is.na(selected_genes)],
      mapping = aes(x = bin, y = 2^log2),
      fill = "#606060", col = "#202020", shape = 21, 
      size = theme_params$size, alpha = theme_params$alpha
    ) +
    geom_point(
      data = targets[!is.na(label) & !is.na(selected_genes)],
      mapping = aes(x = bin, y = 2^log2, fill = selected_genes),
      col = "black", shape = 21, size = theme_params$size_selected, alpha = 0.8
    ) +
    geom_segment(
      data = segments, col = "green", linewidth = 1,
      mapping = aes(x = start, xend = end, y = 2^mean, yend = 2^mean)
    ) +
    geom_text(
      data = cancergenes_summary,
      mapping = aes(x = bin, y = Inf, label = selected_genes, color = selected_genes),
      angle = 45, hjust = 0, size = 2.5, vjust = 0, show.legend = FALSE
    ) +
    scale_fill_hue(l = 70) +
    scale_color_hue(l = 70) +
    scale_y_log10(
      limits = yaxis_params$ylims, 
      breaks = yaxis_params$ybreaks, 
      minor_breaks = yaxis_params$yminorbreaks
    ) +
    scale_x_continuous(
      breaks = chroms_order$mid, minor_breaks = chroms_order$start[-1],
      expand = c(.01, .01), labels = chroms_order$chromosome
    ) +
    coord_cartesian(clip = "off") +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_line(),
      panel.grid.minor.x = element_line(color = "black"),
      axis.line = element_line(),
      axis.ticks = element_line(),
      plot.margin = unit(c(2, 0.5, 0.5, 0.5), "lines")
    )
}

#' Create Position vs Log2 Plot
#'
#' @param targets data.table
#' @param segments data.table
#' @param cancergenes_summary data.table
#' @param chroms data.table
#' @param yaxis_params List from prepare_yaxis_params()
#' @param theme_params List from setup_plot_theme()
#' @return ggplot object
#' @keywords internal
plot_pos_log2 <- function(targets, segments, cancergenes_summary, chroms,
                          yaxis_params, theme_params) {
  ggplot(targets) +
    xlab("Genomic position") +
    ylab("Corrected depth") +
    geom_vline(
      data = cancergenes_summary,
      mapping = aes(xintercept = gpos, color = selected_genes),
      alpha = 0.5, linewidth = 0.5, show.legend = FALSE
    ) +
    geom_point(
      data = targets[!is.na(label) & is.na(selected_genes)],
      mapping = aes(x = gpos, y = 2^log2),
      fill = "#606060", col = "#202020", shape = 21,
      size = theme_params$size, alpha = theme_params$alpha
    ) +
    geom_point(
      data = targets[!is.na(selected_genes)],
      mapping = aes(x = gpos, y = 2^log2, fill = selected_genes),
      shape = 21, size = theme_params$size_selected, alpha = 0.8, col = "black"
    ) +
    geom_segment(
      data = segments, col = "green", linewidth = 1,
      mapping = aes(x = gstart, xend = gstop, y = 2^mean, yend = 2^mean)
    ) +
    geom_text(
      data = cancergenes_summary,
      mapping = aes(x = gpos, y = Inf, label = selected_genes, color = selected_genes),
      angle = 45, hjust = 0, size = 2.5, vjust = 0, show.legend = FALSE
    ) +
    scale_fill_hue(l = 70, na.value = "grey50") +
    scale_color_hue(l = 70) +
    scale_y_log10(
      limits = yaxis_params$ylims, 
      breaks = yaxis_params$ybreaks, 
      minor_breaks = yaxis_params$yminorbreaks
    ) +
    scale_x_continuous(
      breaks = chroms$mid, minor_breaks = chroms$start[-1],
      expand = c(.01, .01), labels = chroms$chromosome
    ) +
    coord_cartesian(clip = "off") +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_line(),
      panel.grid.minor.x = element_line(color = "black"),
      axis.line = element_line(),
      axis.ticks = element_line(),
      plot.margin = unit(c(2, 0.5, 0.5, 0.5), "lines")
    )
}

#' Create SNP-related Plots
#'
#' Creates depth vs allele ratio, grid, nogrid, order, and position plots
#'
#' @param targets data.table with SNP data
#' @param cancergenes_summary data.table
#' @param chroms data.table
#' @param chroms_order data.table
#' @param somatic_plot data.table or NULL
#' @param somatic_labels data.table or NULL
#' @param yaxis_params List from prepare_yaxis_params()
#' @param theme_params List from setup_plot_theme()
#' @return List of ggplot objects
#' @keywords internal
create_snp_plots <- function(targets, cancergenes_summary, chroms, chroms_order,
                             somatic_plot, somatic_labels, yaxis_params, theme_params) {
  p <- list()
  
  # Calculate depth limits. SNP VCFs provide DP; somatic-only runs use target
  # fragment counts as the x-axis background for mutation overlays.
  valid_dp <- targets[!is.na(DP)]$DP
  if (length(valid_dp) > 0) {
    m <- quantile(valid_dp, c(.01, .50, .99), na.rm = TRUE)
  } else {
    m <- quantile(targets[is_target == TRUE]$count, c(.01, .50, .99), na.rm = TRUE)
  }
  
  if (m[1] > m[2] / 3) m[1] <- m[2] / 3
  if (m[3] < m[2] * 3) m[3] <- m[2] * 3
  m[1] <- max(1, m[1])
  
  if (!is.finite(m[1]) || !is.finite(m[3]) || m[1] <= 0 || m[3] <= 0) {
    m <- c(1, 10, 100)
  }
  
  # Depth vs Allele Ratio
  p$depth_alleleratio <- ggplot(targets[type != "background"]) +
    xlab("Depth") +
    ylab("Allele ratio") +
    geom_point(
      data = targets[type != "background" & !is.na(DP) & !is.na(allele_ratio) & is.na(selected_genes)],
      mapping = aes(x = DP, y = allele_ratio),
      fill = "#606060", col = "#202020",
      shape = 21, size = 1, alpha = theme_params$alpha / 2
    ) +
    geom_point(
      data = targets[type != "background" & !is.na(DP) & !is.na(allele_ratio) & !is.na(selected_genes)],
      mapping = aes(x = DP, y = allele_ratio, fill = selected_genes),
      col = "black", shape = 21, size = theme_params$size_selected, alpha = 0.8
    ) +
    (if (!is.null(somatic_plot) && nrow(somatic_plot) > 0) {
      geom_point(
        data = somatic_plot,
        mapping = aes(x = count, y = AF, shape = var_type),
        col = "red", fill = "red", size = theme_params$size_selected, alpha = 0.8
      )
    }) +
    scale_shape_manual(
      values = c("SNV" = 16, "Insertion" = 24, "Deletion" = 25),
      drop = FALSE, guide = "none"
    ) +
    (if (!is.null(somatic_labels) && nrow(somatic_labels) > 0) {
      ggrepel::geom_label_repel(
        data = somatic_labels,
        mapping = aes(x = count, y = AF, label = label),
        fill = "white", col = "darkred", segment.color = "red",
        nudge_y = somatic_labels$nudge,
        size = 2.5, show.legend = FALSE, min.segment.length = 0,
        box.padding = 0.1, label.padding = 0.1, point.padding = 0.1
      )
    }) +
    theme(
      panel.spacing = unit(0, "lines"), 
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    scale_fill_hue(l = 70) +
    scale_x_log10(limits = c(m[1], m[3])) +
    scale_y_continuous(limits = c(-.05, 1.05), breaks = c(0, .25, .5, .75, 1))
  
  # Grid plot
  p$grid <- ggplot(targets) +
    xlim(c(.2, 1.8)) +
    scale_y_continuous(
      limits = c(.5, 1),
      breaks = seq(0.5, .9, by = 0.1),
      labels = c("0.5", "0.6", "0.7", "0.8", "0.9")
    ) +
    xlab("Corrected depth (smooth)") +
    ylab("Major Allele Ratio") +
    geom_point(
      data = targets[!is.na(smooth_log2) & !is.na(maf), .(smooth_log2, maf)],
      aes(x = 2^smooth_log2, y = maf), col = "lightgrey", alpha = .2
    ) +
    geom_text(
      data = unique(targets[, .(chromosome)]),
      mapping = aes(x = .25, y = .95, label = chromosome)
    ) +
    geom_point(
      data = targets[label != "background" & !is.na(smooth_log2) & !is.na(maf) & is.na(selected_genes)],
      aes(x = 2^smooth_log2, y = maf),
      fill = "#606060", col = "#202020", shape = 21, 
      size = theme_params$size, alpha = theme_params$alpha
    ) +
    geom_point(
      data = targets[label != "background" & !is.na(smooth_log2) & !is.na(maf) & !is.na(selected_genes)],
      aes(x = 2^smooth_log2, y = maf, fill = selected_genes),
      shape = 21, size = theme_params$size_selected, alpha = 0.8, col = "black"
    ) +
    facet_wrap(facets = vars(factor(chromosome,
      levels = unique(chromosome),
      ordered = TRUE
    )), ncol = 8) +
    scale_fill_hue(l = 70) +
    theme(
      panel.spacing = unit(0, "lines"),
      strip.background = element_blank(), 
      strip.text.x = element_blank()
    )
  
  # No-grid plot
  p$nogrid <- ggplot(targets) +
    xlim(c(0.2, 1.8)) +
    xlab("Corrected depth (smooth)") +
    ylab("Major allele ratio (smooth)") +
    scale_fill_hue(l = 70) +
    scale_y_continuous(
      limits = c(.5, 1),
      sec.axis = ggplot2::sec_axis(
        ~ 2 * (. - 0.5),
        name = "Homozygous DNA fraction",
        breaks = seq(0, 1, by = 0.1)
      )
    ) +
    geom_point(
      data = targets[!is.na(smooth_log2) & !is.na(maf), .(smooth_log2, maf)],
      aes(x = 2^smooth_log2, y = maf), col = "lightgrey", alpha = .2
    ) +
    geom_point(
      data = targets[label != "background" & !is.na(smooth_log2) & !is.na(maf) & is.na(selected_genes)],
      aes(x = 2^smooth_log2, y = maf),
      fill = "#606060", col = "#202020", shape = 21, 
      size = theme_params$size, alpha = theme_params$alpha
    ) +
    geom_point(
      data = targets[label != "background" & !is.na(smooth_log2) & !is.na(maf) & !is.na(selected_genes)],
      aes(x = 2^smooth_log2, y = maf, fill = selected_genes),
      shape = 21, size = theme_params$size_selected, alpha = 0.8, col = "black"
    )
  
  # Order vs Allele Ratio
  p$order_alleleratio <- ggplot(targets) +
    xlab("Order of genomic position") +
    ylab("Allele ratio") +
    geom_vline(
      data = cancergenes_summary,
      mapping = aes(xintercept = bin, color = selected_genes),
      alpha = 0.5, linewidth = 0.5, show.legend = FALSE
    ) +
    geom_point(
      data = targets[!is.na(label) & is.na(selected_genes)],
      mapping = aes(x = bin, y = allele_ratio),
      fill = "#606060", col = "#202020", shape = 21, 
      size = theme_params$size, alpha = theme_params$alpha
    ) +
    geom_point(
      data = targets[!is.na(label) & !is.na(selected_genes)],
      mapping = aes(x = bin, y = allele_ratio, fill = selected_genes),
      col = "black", shape = 21, size = theme_params$size_selected, alpha = 0.8
    ) +
    (if (!is.null(somatic_plot) && nrow(somatic_plot) > 0) {
      geom_point(
        data = somatic_plot,
        mapping = aes(x = bin, y = AF, shape = var_type),
        col = "red", fill = "red", size = theme_params$size_selected, alpha = 0.8
      )
    }) +
    scale_shape_manual(
      values = c("SNV" = 16, "Insertion" = 24, "Deletion" = 25),
      drop = FALSE, guide = "none"
    ) +
    (if (!is.null(somatic_labels) && nrow(somatic_labels) > 0) {
      ggrepel::geom_label_repel(
        data = somatic_labels,
        mapping = aes(x = bin, y = AF, label = label),
        fill = "white", col = "darkred", segment.color = "red",
        nudge_y = somatic_labels$nudge,
        size = 2.5, show.legend = FALSE, min.segment.length = 0,
        box.padding = 0.1, label.padding = 0.1, point.padding = 0.1
      )
    }) +
    geom_text(
      data = cancergenes_summary,
      mapping = aes(x = bin, y = Inf, label = selected_genes, color = selected_genes),
      angle = 45, hjust = 0, size = 2.5, vjust = 0, show.legend = FALSE
    ) +
    scale_fill_hue(l = 70) +
    scale_color_hue(l = 70) +
    scale_x_continuous(
      breaks = chroms_order$mid, minor_breaks = chroms_order$start[-1],
      expand = c(.01, .01), labels = chroms_order$chromosome
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, .25, .5, .75, 1)) +
    coord_cartesian(clip = "off") +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_line(),
      panel.grid.minor.x = element_line(color = "black"),
      axis.line = element_line(),
      axis.ticks = element_line(),
      plot.margin = unit(c(2, 0.5, 0.5, 0.5), "lines")
    )
  
  # Position vs Allele Ratio
  p$pos_alleleratio <- ggplot(targets) +
    xlab("Genomic position") +
    ylab("Allele ratio") +
    geom_vline(
      data = cancergenes_summary,
      mapping = aes(xintercept = gpos, color = selected_genes),
      alpha = 0.5, linewidth = 0.5, show.legend = FALSE
    ) +
    geom_point(
      data = targets[!is.na(label) & is.na(selected_genes)],
      mapping = aes(x = gpos, y = allele_ratio),
      fill = "#606060", col = "#202020", shape = 21, 
      size = theme_params$size, alpha = theme_params$alpha
    ) +
    geom_point(
      data = targets[!is.na(label) & !is.na(selected_genes)],
      mapping = aes(x = gpos, y = allele_ratio, fill = selected_genes),
      col = "black", shape = 21, size = theme_params$size_selected, alpha = 0.8
    ) +
    (if (!is.null(somatic_plot) && nrow(somatic_plot) > 0) {
      geom_point(
        data = somatic_plot,
        mapping = aes(x = gpos, y = AF, shape = var_type),
        col = "red", fill = "red", size = 1.2, alpha = 0.8
      )
    }) +
    scale_shape_manual(
      values = c("SNV" = 16, "Insertion" = 24, "Deletion" = 25),
      drop = FALSE, guide = "none"
    ) +
    (if (!is.null(somatic_labels) && nrow(somatic_labels) > 0) {
      ggrepel::geom_label_repel(
        data = somatic_labels,
        mapping = aes(x = gpos, y = AF, label = label),
        fill = "white", col = "darkred", segment.color = "red",
        nudge_y = somatic_labels$nudge,
        size = 2.5, show.legend = FALSE, min.segment.length = 0,
        box.padding = 0.1, label.padding = 0.1, point.padding = 0.1
      )
    }) +
    geom_text(
      data = cancergenes_summary,
      mapping = aes(x = gpos, y = Inf, label = selected_genes, color = selected_genes),
      angle = 45, hjust = 0, size = 2.5, vjust = 0, show.legend = FALSE
    ) +
    scale_fill_hue(l = 70) +
    scale_color_hue(l = 70) +
    scale_x_continuous(
      breaks = chroms$mid, minor_breaks = chroms$start[-1],
      expand = c(.01, .01), labels = chroms$chromosome
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, .25, .5, .75, 1)) +
    coord_cartesian(clip = "off") +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_line(),
      panel.grid.minor.x = element_line(color = "black"),
      axis.line = element_line(),
      axis.ticks = element_line(),
      plot.margin = unit(c(2, 0.5, 0.5, 0.5), "lines")
    )
  
  p
}

#' Create GIS Plot
#'
#' @param gis_table data.table or NULL
#' @return ggplot object
#' @keywords internal
plot_gis <- function(gis_table) {
  if (!is.null(gis_table)) {
    p <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = gis_table, 
        aes(x = fraction, y = predicted_gis, col = "predicted_gis", group = 1), 
        linewidth = 2
      ) +
      ggplot2::geom_line(
        data = gis_table, 
        aes(x = fraction, y = local_cnv * 10, col = "local_cnv"), 
        lty = 1, linewidth = 1
      ) +
      ggplot2::geom_line(
        data = gis_table, 
        aes(x = fraction, y = focal_gain * 10, col = "focal_gain"), 
        lty = 1, linewidth = 1
      ) +
      ggplot2::geom_line(
        data = gis_table, 
        aes(x = fraction, y = loh * 10, col = "loh"), 
        lty = 1, linewidth = 1
      )
      
    if ("custom_HRD" %in% names(gis_table) && any(!is.na(gis_table$custom_HRD))) {
      p <- p + ggplot2::geom_line(
        data = gis_table,
        aes(x = fraction, y = custom_HRD, col = "custom_HRD"),
        linetype = 2,
        linewidth = 1.5
      )
    }
    
    p <- p + ggplot2::scale_color_manual(
        name = "Feature (x10)",
        values = c("predicted_gis" = "black", "local_cnv" = "#F8766D", "focal_gain" = "#00BA38", "loh" = "#619CFF", "custom_HRD" = "#8B0000"),
        breaks = c("predicted_gis", "local_cnv", "focal_gain", "loh", "custom_HRD"),
        labels = c("predicted_gis" = "GIS score", "local_cnv" = "local_cnv", "focal_gain" = "focal_gain", "loh" = "loh", "Custom HRD")
      ) +
      scale_x_continuous(
        name = "Tumor DNA fraction",
        breaks = seq(0, 1, .05), minor_breaks = seq(0, 1, .01),
        labels = seq(0, 100, 5)
      ) +
      scale_y_continuous(
        name = "GMCK GIS score", limits = c(0, 100),
        breaks = seq(0, 100, 10),
        sec.axis = sec_axis(~., breaks = seq(0, 100, 10))
      ) +
      ggplot2::geom_hline(yintercept = 42, lty = 2, col = "#00000050", linewidth = 2) +
      theme(legend.position = "bottom")
    
    return(p)
  } else {
    ggplot2::ggplot() +
      ggplot2::geom_text(aes(x = 0.5, y = 0.5, label = "No GIS Data")) +
      theme(
        axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank()
      )
  }
}

#' Add Legend Formatting to Plots
#'
#' @param plot_list List of ggplot objects
#' @return Modified list with consistent legend formatting
#' @keywords internal
format_plot_legends <- function(plot_list) {
  for (i in 1:length(plot_list)) {
    plot_list[[i]] <- plot_list[[i]] + 
      guides(fill = guide_legend(override.aes = list(shape = 21, size = 3))) + 
      labs(fill = NULL)
  }
  plot_list
}

#' Arrange and Save Final Plot
#'
#' @param plot_list List of ggplot objects
#' @param use_snp Logical: whether SNP plots are included
#' @param title Character: plot title
#' @param output_file Character path or NULL
#' @return ggplot object (invisibly if saved)
#' @keywords internal
arrange_and_save_plots <- function(plot_list, use_snp, title, output_file) {
  if (use_snp) {
    layout <- "AABBBBBB
               CCDDDDDD
               EEFFFFFF
               GGHHHHHH"
    fig <- plot_list$gc_rawdepth + plot_list$order_rawdepth +
      plot_list$gc_log2 + plot_list$order_log2 +
      plot_list$depth_alleleratio + plot_list$order_alleleratio +
      plot_list$nogrid + plot_list$grid +
      patchwork::plot_layout(design = layout, guides = "collect") &
      theme(legend.position = "none")
  } else {
    layout <- "ABBB
               CDDD
               EEEE"
    fig <- plot_list$gc_rawdepth + plot_list$order_rawdepth +
      plot_list$gc_log2 + plot_list$order_log2 +
      plot_list$pos_log2 +
      patchwork::plot_layout(design = layout, guides = "collect") &
      theme(legend.position = "none")
  }
  
  # Add title and footer
  pa <- patchwork::plot_annotation(
    title = title,
    caption = paste0("Jumble ", get_jumble_version(), " on ", format(Sys.time(), "%a %b %e %Y, %H:%M")),
    theme = ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 1, size = 8, color = "grey50")
    )
  )
  fig <- fig + pa
  
  # Save if specified
  if (!is.null(output_file)) {
    png(file = output_file, width = 1600, height = 1300, res = 100)
    suppressWarnings(print(fig))
    dev.off()
  }
  
  invisible(fig)
}

#' Plot Jumble Results (Refactored)
#'
#' Generates plots for Jumble analysis results.
#'
#' @param targets Normalized targets data.table (must have label, selected_genes).
#' @param segments Segments data.table.
#' @param reference Reference object (for chromlength).
#' @param snp_table Optional SNP table.
#' @param output_file Path to save the plot (PNG).
#' @param title Plot title.
#' @param gis_table Optional GIS results table.
#' @param somatic_table Optional table of somatic mutations.
#' @return The ggplot object (invisibly).
#' @keywords internal
plot_results <- function(targets, segments, reference = NULL, snp_table = NULL,
                         output_file = NULL, title = "", gis_table = NULL, 
                         somatic_table = NULL) {
  # 1. Setup
  annotate_targets(targets)
  theme_params <- setup_plot_theme()
  axis_params <- prepare_depth_axis()
  yaxis_params <- prepare_yaxis_params(targets)
  
  # 2. Filter and prepare data
  chrom_levels <- c(as.character(1:22), "X", "Y")
  targets <- targets[chromosome %in% chrom_levels]
  
  # 3. Prepare chromosome data
  chr_data <- prepare_chromosome_data(targets, reference)
  chroms <- chr_data$chroms
  chroms_order <- chr_data$chroms_order
  
  # 4. Add genomic positions
  pos_data <- add_genomic_positions(targets, segments, chroms)
  targets <- pos_data$targets
  segments <- pos_data$segments
  
  # 5. Prepare cancer gene summaries
  cancergenes_summary <- prepare_cancergenes_summary(targets)
  
  # 6. Prepare somatic data
  somatic_data <- prepare_somatic_data(somatic_table, targets, reference)
  somatic_plot <- somatic_data$somatic_plot
  somatic_labels <- somatic_data$somatic_labels
  
  # 7. Initialize plot list
  p <- list()
  
  # 8. Create basic plots
  p$gc_rawdepth <- plot_gc_rawdepth(targets, theme_params, axis_params)
  p$order_rawdepth <- plot_order_rawdepth(targets, cancergenes_summary, chroms_order,
                                           theme_params, axis_params)
  p$pos_rawdepth <- plot_pos_rawdepth(targets, cancergenes_summary, chroms,
                                       theme_params, axis_params)
  p$gc_log2 <- plot_gc_log2(targets, yaxis_params, theme_params)
  p$order_log2 <- plot_order_log2(targets, segments, cancergenes_summary, chroms_order,
                                   yaxis_params, theme_params)
  p$pos_log2 <- plot_pos_log2(targets, segments, cancergenes_summary, chroms,
                               yaxis_params, theme_params)
  
  # 9. Create variant panels if SNPs or somatic mutations are available.
  use_variant_panels <- !is.null(snp_table) || (!is.null(somatic_plot) && nrow(somatic_plot) > 0)
  if (use_variant_panels) {
    targets <- prepare_snp_data(snp_table, targets)
    snp_plots <- create_snp_plots(targets, cancergenes_summary, chroms, chroms_order,
                                   somatic_plot, somatic_labels, yaxis_params, theme_params)
    p <- c(p, snp_plots)
    p$gis_plot <- plot_gis(gis_table)
  }
  
  # 10. Format legends
  p <- format_plot_legends(p)
  
  # 11. Arrange and save
  arrange_and_save_plots(p, use_variant_panels, title, output_file)
}

#' Plot GIS Score
#'
#' Plots the Genomic Instability Score analysis:
#' 1. Smooth corrected depth vs Allele Ratio (top panel)
#' 2. GIS metrics vs Tumor Fraction (bottom panel)
#'
#' @param gis_table Data.table result from compute_gis_table.
#' @param targets Processed targets (must have smooth_log2, maf).
#' @param output_file Path to save the plot (PNG).
#' @param title Plot title.
#' @return NULL (saves file)
#' @import ggplot2 patchwork
#' @keywords internal
plot_gis_score <- function(gis_table, targets, output_file, title = "GIS Analysis") {
  # Ensure valid data
  if (is.null(gis_table) || nrow(gis_table) == 0) {
    warning("No GIS table provided for plotting.")
    return(NULL)
  }

  # Melt table to long format for plotting (if not already)
  if (!"feature" %in% names(gis_table)) {
    # Check if we have the expected columns
    cols <- c("fraction", "predicted_gis", "local_cnv", "focal_gain", "loh")
    if (all(cols %in% names(gis_table))) {
      # Ensure all columns are numeric
      gis_table[, `:=`(
        predicted_gis = as.numeric(predicted_gis),
        local_cnv = as.numeric(local_cnv),
        focal_gain = as.numeric(focal_gain),
        loh = as.numeric(loh),
        fraction = as.numeric(fraction)
      )]
      measure_vars <- c("predicted_gis", "local_cnv", "focal_gain", "loh")
      has_custom <- "custom_HRD" %in% names(gis_table) && any(!is.na(gis_table$custom_HRD))
      if (has_custom) {
        measure_vars <- c(measure_vars, "custom_HRD")
      }
      
      plot_data <- data.table::melt(
        gis_table,
        id.vars = "fraction",
        measure.vars = measure_vars,
        variable.name = "feature",
        value.name = "value"
      )
    } else {
      warning("GIS table missing expected columns for plotting.")
      return(NULL)
    }
  } else {
    plot_data <- gis_table
  }

  # 1. Prepare Data ----------------------------------------------------------
  # --- Top Panel: p (Smooth Depth vs MAF) ---
  # Needs: targets with smooth_log2, maf, label, selected_genes

  # Basic plot setup
  size <- ifelse(nrow(targets) > 50000, 0.4, 0.8) # Adjust point size based on data density
  alpha <- 0.6

  p1 <- ggplot(targets) +
    xlim(c(0.2, 1.8)) +
    xlab("Corrected depth (smooth)") +
    ylab("Major allele ratio (smooth)") +
    scale_fill_hue(l = 70) +
    scale_y_continuous(
      limits = c(.5, 1),
      sec.axis = ggplot2::sec_axis(
        ~ 2 * (. - 0.5),
        name = "Homozygous DNA fraction",
        breaks = seq(0, 1, by = 0.1)
      )
    ) +
    geom_point(
      data = targets[, .(smooth_log2, maf)],
      aes(x = 2^smooth_log2, y = maf), col = "lightgrey", alpha = .2
    ) +
    # Background
    geom_point(
      data = targets[label != "background" & is.na(selected_genes)],
      aes(x = 2^smooth_log2, y = maf),
      fill = "#606060", col = "#202020", shape = 21, size = size, alpha = alpha
    ) +
    # Selected
    geom_point(
      data = targets[label != "background" & !is.na(selected_genes)],
      aes(x = 2^smooth_log2, y = maf, fill = selected_genes),
      shape = 21, size = size + 0.5, alpha = 0.8, col = "black"
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )


  # 3. Bottom Panel (Curves) -------------------------------------------------
  # Minimal cancer fraction to display (matching frankenscript)
  plot_data <- plot_data[fraction >= 0.2]

  # If estimated fraction exists in metadata, use it? For now assume it's just the table data.

  p2 <- ggplot() +
    geom_line(data = plot_data[feature == "predicted_gis"], aes(x = fraction, y = value, col = "predicted_gis"), linewidth = 2) +
    geom_line(
      data = plot_data[feature %in% c("local_cnv", "focal_gain", "loh")],
      aes(x = fraction, y = value, col = feature), linetype = 1, linewidth = 1
    )
    
  if ("custom_HRD" %in% plot_data$feature) {
    p2 <- p2 + geom_line(
      data = plot_data[feature == "custom_HRD"],
      aes(x = fraction, y = value, col = "custom_HRD"),
      linetype = 2, linewidth = 1.5
    )
  }
  
  p2 <- p2 + scale_color_manual(
      name = NULL,
      values = c("predicted_gis" = "black", "local_cnv" = "#F8766D", "focal_gain" = "#00BA38", "loh" = "#619CFF", "custom_HRD" = "#8B0000"),
      breaks = c("predicted_gis", "local_cnv", "focal_gain", "loh", "custom_HRD"),
      labels = c("predicted_gis" = "GIS score", "local_cnv" = "local_cnv", "focal_gain" = "focal_gain", "loh" = "loh", "custom_HRD" = "Custom HRD")
    ) +
    scale_x_continuous(
      name = "Tumor DNA fraction",
      breaks = seq(0, 1, .05), minor_breaks = seq(0, 1, .01),
      labels = seq(0, 100, 5)
    ) +
    scale_y_continuous(
      name = "GMCK GIS score", limits = c(0, 100),
      breaks = seq(0, 100, 5), minor_breaks = seq(0, 100, 1),
      sec.axis = sec_axis(~., breaks = seq(0, 100, 5))
    ) +
    geom_hline(yintercept = 42, linetype = 2, col = "#00000050", linewidth = 2) +
    theme_bw() +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      legend.position = "bottom",
      legend.title = element_blank()
    )

  # 4. Combine & Save --------------------------------------------------------
  final_plot <- p1 + p2 +
    plot_annotation(
      title = title,
      caption = paste0("Jumble ", utils::packageVersion("Jumble"), " on ", format(Sys.time(), "%a %b %e %Y, %H:%M")),
      theme = theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
                    plot.caption = ggplot2::element_text(hjust = 1, size = 8, color = "grey50"))
    )

  # Save if requested, otherwise return the plot for report embedding.
  if (!is.null(output_file)) {
    suppressWarnings(ggsave(output_file, plot = final_plot, width = 16, height = 12, dpi = 300))
  }

  invisible(final_plot)
}


#' Plot MSI Indel Count vs Minimum VAF Threshold
#'
#' Creates a diagnostic plot showing how MSI-like indel counts change
#' across minimum VAF thresholds. Helps distinguish MSI signal from noise.
#'
#' @param somatic data.table with AF and MSI columns.
#' @param output_file Path to save the plot (PNG).
#' @param title Plot title.
#' @return NULL (saves file).
#' @import ggplot2
#' @keywords internal
plot_msi_vaf <- function(somatic, output_file, title = "MSI VAF Analysis") {
  if (is.null(somatic) || nrow(somatic) == 0 || !"MSI" %in% names(somatic)) {
    return(NULL)
  }

  # Only indels (MSI is NA for SNVs)
  indels <- somatic[!is.na(MSI)]
  if (nrow(indels) == 0) return(NULL)

  # Define VAF thresholds to sweep
  vaf_thresholds <- seq(0.01, 0.50, by = 0.01)

  # Build summary data for each threshold
  plot_data <- data.table::rbindlist(lapply(vaf_thresholds, function(vaf_min) {
    sub <- indels[AF >= vaf_min]
    data.table::data.table(
      vaf_threshold = vaf_min,
      category = c("Total indels", "Mono repeats", "Di repeats", "Tri repeats"),
      count = c(
        nrow(sub),
        sum(sub$MSI == 1),
        sum(sub$MSI == 2),
        sum(sub$MSI == 3)
      )
    )
  }))

  if (nrow(plot_data) == 0) return(NULL)

  # Always show all categories (including zero-count ones) for consistent legend
  cat_order <- c("Total indels", "Mono repeats", "Di repeats", "Tri repeats")
  plot_data[, category := factor(category, levels = cat_order)]

  # Colors
  cat_colors <- c(
    "Total indels"  = "grey50",
    "Mono repeats"  = "#E41A1C",
    "Di repeats"    = "#377EB8",
    "Tri repeats"   = "#4DAF4A"
  )

  # Log Y axis with minimum range of 400
  y_max <- max(400, max(plot_data$count, na.rm = TRUE) * 1.1)

  # Replace 0 with 0.5 for log scale display
  plot_data[, plot_count := ifelse(count == 0, 0.5, count)]

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = vaf_threshold, y = plot_count,
    color = category, linewidth = category
  )) +
    ggplot2::geom_line() +
    ggplot2::scale_x_continuous(
      name = "Minimum VAF threshold",
      breaks = seq(0, 0.5, by = 0.05),
      labels = function(x) paste0(x * 100, "%"),
      expand = c(0.01, 0.01)
    ) +
    ggplot2::scale_y_log10(
      name = "Number of variants",
      breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 400, 1000, 2000, 5000)
    ) +
    ggplot2::coord_cartesian(ylim = c(0.5, y_max)) +
    ggplot2::scale_color_manual(values = cat_colors, drop = FALSE) +
    ggplot2::scale_linewidth_manual(
      values = c("Total indels" = 0.8,
                 "Mono repeats" = 1, "Di repeats" = 1, "Tri repeats" = 1),
      drop = FALSE
    ) +
    ggplot2::labs(title = title, color = NULL, linewidth = NULL,
                  caption = paste0("Jumble ", get_jumble_version(), " on ", format(Sys.time(), "%a %b %e %Y, %H:%M"))) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey90"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5),
      plot.caption = ggplot2::element_text(hjust = 1, size = 8, color = "grey50")
    ) +
    ggplot2::guides(linewidth = "none")

  suppressWarnings(ggplot2::ggsave(output_file, plot = p, width = 10, height = 7, dpi = 300))
}

#' Plot TMB Estimate vs Minimum VAF Threshold
#'
#' @param somatic data.table with AF and variants mapping to valid bins.
#' @param targets targets data.table.
#' @param output_file Path to save the plot (PNG).
#' @param title Plot title.
#' @return NULL (saves file).
#' @import ggplot2
#' @keywords internal
plot_tmb_vaf <- function(somatic, targets, output_file, title = "TMB VAF Analysis") {
  if (is.null(somatic) || nrow(somatic) == 0) return(NULL)

  # Mask poorly covered / background
  valid_bin_indices <- which(targets$is_target == TRUE & targets$count > (0.2 * median(targets$count[targets$is_target == TRUE], na.rm = TRUE)) & targets$count > 50)
  if (length(valid_bin_indices) == 0) return(NULL)
  
  target_mb <- sum(targets[valid_bin_indices, end - start], na.rm = TRUE) / 1e6
  if (target_mb <= 0) return(NULL)

  if ("is_rare_snp" %in% names(somatic)) {
    tmb_pool <- somatic[!is.na(bin) & bin %in% valid_bin_indices & (is_rare_snp == FALSE | is.na(is_rare_snp))]
  } else {
    tmb_pool <- somatic[!is.na(bin) & bin %in% valid_bin_indices]
  }
  
  if (nrow(tmb_pool) == 0) return(NULL)

  tmb_pool[, is_indel := nchar(REF) != nchar(ALT)]
  
  vaf_thresholds <- seq(0.01, 0.50, by = 0.01)
  
  plot_data <- data.table::rbindlist(lapply(vaf_thresholds, function(vaf_min) {
    sub <- tmb_pool[AF >= vaf_min]
    tmb_snv <- sum(!sub$is_indel)
    tmb_indel <- sum(sub$is_indel)
    n_estim <- round(tmb_indel + (tmb_snv * 0.70))
    est <- n_estim / target_mb
    
    data.table::data.table(
      vaf_threshold = vaf_min,
      category = c("Substitutions", "Indels", "TMB Estimate (/Mb)"),
      value = c(tmb_snv, tmb_indel, est)
    )
  }))
  
  if (nrow(plot_data) == 0) return(NULL)
  
  cat_order <- c("Substitutions", "Indels", "TMB Estimate (/Mb)")
  plot_data[, category := factor(category, levels = cat_order)]
  
  cat_colors <- c(
    "Substitutions" = "#FF7F00",
    "Indels"        = "#984EA3",
    "TMB Estimate (/Mb)" = "#000000"
  )
  
  y_max <- max(10, max(plot_data$value, na.rm = TRUE) * 1.1)
  
  plot_data[, plot_value := ifelse(value == 0, 0.5, value)]

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = vaf_threshold, y = plot_value,
    color = category, linewidth = category, linetype = category
  )) +
    ggplot2::geom_line() +
    ggplot2::scale_x_continuous(
      name = "Minimum VAF threshold",
      breaks = seq(0, 0.5, by = 0.05),
      labels = function(x) paste0(x * 100, "%"),
      expand = c(0.01, 0.01)
    ) +
    ggplot2::scale_y_log10(
      name = "Count / Score",
      breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)
    ) +
    ggplot2::coord_cartesian(ylim = c(0.5, y_max)) +
    ggplot2::scale_color_manual(values = cat_colors, drop = FALSE) +
    ggplot2::scale_linewidth_manual(
      values = c("Substitutions" = 1, "Indels" = 1, "TMB Estimate (/Mb)" = 2.0),
      drop = FALSE
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Substitutions" = "solid", "Indels" = "solid", "TMB Estimate (/Mb)" = "solid"),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = title, color = NULL, linewidth = NULL, linetype = NULL,
      caption = paste0("Jumble ", get_jumble_version(), " on ", format(Sys.time(), "%a %b %e %Y, %H:%M"))
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey90"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5),
      plot.caption = ggplot2::element_text(hjust = 1, size = 8, color = "grey50")
    ) +
    ggplot2::guides(linewidth = "none", linetype = "none")

  suppressWarnings(ggplot2::ggsave(output_file, plot = p, width = 10, height = 7, dpi = 300))
}