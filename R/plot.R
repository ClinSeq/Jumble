
#' Annotate Targets for Plotting
#'
#' Adds 'label' and 'selected_genes' columns to targets.
#'
#' @param targets data.table.
#' @return Modifies targets in-place.
#' @export
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

#' Plot Jumble Results
#'
#' Generates plots for Jumble analysis results matching original script.
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
#' @importFrom ggplot2 ggplot aes geom_point scale_fill_manual scale_y_log10 scale_x_log10 scale_x_continuous scale_y_continuous sec_axis theme element_blank element_line labs guides guide_legend geom_segment facet_wrap vars unit geom_text xlim xlab ylab theme_set theme_bw geom_line geom_hline geom_vline scale_fill_hue scale_color_hue coord_cartesian
#' @importFrom ggrepel geom_label_repel
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom grDevices png dev.off
#' @importFrom stats quantile
#' @export
plot_results <- function(targets, segments, reference = NULL, snp_table = NULL,
                         output_file = NULL, title = "", gis_table = NULL, somatic_table = NULL) {
  # 1. Setup & Theme ---------------------------------------------------------
  # Set theme to match original script (line 68)
  ggplot2::theme_set(ggplot2::theme_bw())

  # Define colors matching original
  colorvalues <- c(
    "background" = "#1010D0", "dense" = "#F8566D",
    "exonic" = "#00BFA4", "sparse" = "#000000", "bin" = "#000000"
  )

  size <- 1
  size_selected <- 1.1
  pointcolor <- "#000000"
  alpha <- 0.3

  # 2. Prepare Data ----------------------------------------------------------
  # Ensure targets are annotated (if not already)
  annotate_targets(targets)

  # Filter to standard chromosomes only
  # This removes alternative contigs (1_KI...) which cause plotting artifacts
  chrom_levels <- c(as.character(1:22), "X", "Y")
  targets <- targets[chromosome %in% chrom_levels]

  # 3. Plot Setup (Axes) -----------------------------------------------------
  ylims <- c(
    min(0.4, min(2^targets[chromosome != "Y"]$smooth_log2, na.rm = TRUE)),
    max(3, max(2^targets$smooth_log2, na.rm = TRUE))
  )

  ybreaks <- c(.5, .75, 1, 1.5, 2, 3, 4, 6, 8)
  ybreaks <- ybreaks[ybreaks >= ylims[1] & ybreaks <= ylims[2]]
  yminorbreaks <- c(1.25, 1.75)

  # Raw depth limits
  limits <- c(1e1, 1e5)
  limits_labels <- c("10", "30", "100", "300", "1k", "3k", "10k", "30k", "100k")
  limits_breaks <- c(1e1, 3e1, 1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5)

  # Chromosome positions (by genomic position)
  if (!is.null(reference) && !is.null(reference$chromlength)) {
    chroms <- data.table::data.table(
      chromosome = names(reference$chromlength),
      length = as.numeric(reference$chromlength)
    )
  } else {
    chroms <- targets[, .(length = max(end)), by = chromosome]
  }

  chroms[, chromosome := clean_chrom_names(as.character(chromosome))]
  chrom_levels <- c(as.character(1:22), "X", "Y")
  chroms <- chroms[match(chrom_levels, chromosome)]
  chroms <- chroms[!is.na(chromosome)]

  # Cumulative positions for genomic plots
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

  # Add genomic positions to targets
  targets[, gpos := as.numeric(start)]
  for (chr in unique(targets$chromosome)[-1]) {
    targets[chromosome == chr, gpos := gpos + sum(chroms[1:(which(chroms$chromosome == chr) - 1)]$length)]
  }

  # Add genomic positions to segments
  segments[, gstart := as.double(start_pos)]
  segments[, gstop := as.double(end_pos)]
  for (chr in unique(segments$chromosome)[-1]) {
    segments[chromosome == chr, gstart := gstart + sum(chroms[1:(which(chroms$chromosome == chr) - 1)]$length)]
    segments[chromosome == chr, gstop := gstop + sum(chroms[1:(which(chroms$chromosome == chr) - 1)]$length)]
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

  # Filter to standard chromosomes only
  chroms_order <- chroms_order[chromosome %in% chrom_levels]

  # Prepare cancergenes_summary for labeling (one point per gene)
  cancergenes_summary <- targets[!is.na(selected_genes), .(
    bin = median(bin, na.rm = TRUE),
    gpos = median(gpos, na.rm = TRUE),
    count = median(count, na.rm = TRUE),
    log2 = median(smooth_log2, na.rm = TRUE)
  ), by = selected_genes]

  if (nrow(cancergenes_summary) > 0) {
    # Calculate median count for nudge logic
    m_depth <- median(targets[type != "background" & chromosome %in% 1:22]$count, na.rm = TRUE)
    if (is.na(m_depth)) m_depth <- 100 # Fallback

    cancergenes_summary[, nudge := 0.1]

    # Apply nudge logic similar to frankenscript (in log units)
    if (nrow(cancergenes_summary) >= 2) {
      cancergenes_summary[seq(2, .N, 2), nudge := -1]
    }
    cancergenes_summary[count < m_depth * 0.8, nudge := 1]
    cancergenes_summary[count > m_depth * 1.2, nudge := -1]
  }

  # Prepare Somatic Data for Plotting
  somatic_plot <- NULL
  somatic_labels <- NULL
  if (!is.null(somatic_table) && nrow(somatic_table) > 0) {
    # We need to map bins to gpos/log2/count if not already present
    # somatic_table has chromosome, start, AF, bin

    # Copy to avoid modifying original by reference if it's used elsewhere
    somatic_plot <- data.table::copy(somatic_table)

    # Only keep somatic variants that fall into valid bins
    somatic_plot <- somatic_plot[!is.na(bin)]

    if (nrow(somatic_plot) > 0) {
      # Match bins to targets
      # targets should be keyed by bin or we assume unique bins
      # Use character subsetting to avoid closure/function name conflicts
      cols_needed <- c("bin", "gpos", "smooth_log2", "log2", "count")
      # Ensure columns exist
      cols_exist <- cols_needed %in% names(targets)
      if (all(cols_exist)) {
        t_lookup <- targets[, cols_needed, with = FALSE]
        data.table::setkey(t_lookup, bin)
        data.table::setkey(somatic_plot, bin)

        # Left join somatic_plot with t_lookup
        somatic_plot <- t_lookup[somatic_plot]
        
        # Prepare Somatic Labels (matching frankenscript logic)
        somatic_labels <- NULL
        cancergenes_list <- NULL
        
        if (!is.null(reference) && !is.null(reference$cancergenes_clinseq)) {
            cancergenes_list <- reference$cancergenes_clinseq$gene
        }
        
        if (is.null(cancergenes_list) && !is.null(targets$selected_genes)) {
            cancergenes_list <- unique(na.omit(as.character(targets$selected_genes)))
        }
    
        if (nrow(somatic_plot) > 0 && !is.null(cancergenes_list)) { 
          req_cols <- c("SYMBOL", "CANONICAL", "effect", "Consequence", "Protein_position", "Amino_acids")
          if (all(req_cols %in% names(somatic_plot))) {
            
            somatic_labels <- somatic_plot[SYMBOL %in% cancergenes_list & CANONICAL == "YES" & effect %in% c("hotspot", "high-impact")]
            
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

      } else {
        warning("Targets missing columns for somatic mapping: ", paste(cols_needed[!cols_exist], collapse = ", "))
        somatic_plot <- NULL
      }

      # For AF plots: y = AF
      # For Log2 plots: y = 2^log2
      # For Grid plots: y = folded MAF = abs(AF - 0.5) + 0.5
    }
  }

  # Initialize plot list
  p <- list()

  # 5. Plot: GC vs Depth -----------------------------------------------------
  # 1. GC vs Raw Depth (TOP LEFT)
  # Filter count > 0 to avoid log10(0) warnings
  p$gc_rawdepth <- ggplot(targets[count > 0]) +
    xlab("Target GC content") +
    ylab("Fragments") +
    xlim(c(.16, .84)) +
    # Layer 1: Background only (no selected genes)
    geom_point(
      mapping = aes(x = gc, y = count),
      col = pointcolor, shape = 21, size = 1, alpha = alpha / 2
    ) +
    facet_wrap(facets = vars(label), nrow = 1) +
    theme(
      panel.spacing = unit(0, "lines"), strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    geom_text(
      data = unique(targets[, .(label)])[label == "target", label := "sparse"],
      mapping = aes(label = label, x = 0.16, y = Inf), hjust = 0, vjust = 1.5
    ) +
    scale_y_log10(
      limits = limits, breaks = limits_breaks,
      minor_breaks = yminorbreaks, labels = limits_labels
    )

  # 6. Plot: Order vs Depth --------------------------------------------------
  # 2. Order vs Raw Depth (TOP RIGHT)
  p$order_rawdepth <- ggplot(targets) +
    xlab("Order of genomic position") +
    ylab("Fragments") +
    # Layer 0: Vertical lines (behind points)
    geom_vline(
      data = cancergenes_summary,
      mapping = aes(xintercept = bin, color = selected_genes),
      alpha = 0.5, linewidth = 0.5, show.legend = FALSE
    ) +
    # Layer 1: Background
    geom_point(
      data = targets[!is.na(label) & is.na(selected_genes) & count > 0],
      mapping = aes(x = bin, y = count),
      fill = "#606060", col = "#202020", shape = 21, size = size, alpha = alpha
    ) +
    # Layer 2: Selected
    geom_point(
      data = targets[!is.na(label) & !is.na(selected_genes) & count > 0],
      mapping = aes(x = bin, y = count, fill = selected_genes),
      col = "black", shape = 21, size = size_selected, alpha = 0.8
    ) +
    # Layer 3: Text labels (on top)
    geom_text(
      data = cancergenes_summary,
      mapping = aes(x = bin, y = Inf, label = selected_genes, color = selected_genes),
      angle = 45, hjust = 0, size = 2.5, vjust = 0, show.legend = FALSE
    ) +
    scale_fill_hue(l = 70) +
    scale_color_hue(l = 70) +
    scale_y_log10(
      limits = limits, breaks = limits_breaks,
      minor_breaks = yminorbreaks, labels = limits_labels
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

  # 7. Plot: Position vs Depth -----------------------------------------------
  # Genomic Position vs Raw Depth
  p$pos_rawdepth <- ggplot(targets) +
    xlab("Genomic position") +
    ylab("Fragments") +
    # Layer 0: Vertical lines
    geom_vline(
      data = cancergenes_summary,
      mapping = aes(xintercept = gpos, color = selected_genes),
      alpha = 0.5, linewidth = 0.5, show.legend = FALSE
    ) +
    # Layer 1: Background
    geom_point(
      data = targets[!is.na(label) & is.na(selected_genes) & count > 0],
      mapping = aes(x = gpos, y = count),
      fill = "#606060", col = "#202020", shape = 21, size = size, alpha = alpha
    ) +
    # Layer 2: Selected
    geom_point(
      data = targets[!is.na(label) & !is.na(selected_genes) & count > 0],
      mapping = aes(x = gpos, y = count, fill = selected_genes),
      col = "black", shape = 21, size = size_selected, alpha = 0.8
    ) +
    # Layer 3: Labels
    geom_text(
      data = cancergenes_summary,
      mapping = aes(x = gpos, y = Inf, label = selected_genes, color = selected_genes),
      angle = 45, hjust = 0, size = 2.5, vjust = 0, show.legend = FALSE
    ) +
    scale_fill_hue(l = 70) +
    scale_color_hue(l = 70) +
    scale_y_log10(
      limits = limits, breaks = limits_breaks,
      minor_breaks = yminorbreaks, labels = limits_labels
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

  # 3. GC vs Log2 (MIDDLE LEFT)
  # Updated to 2-layer coloring
  # 7. Plot: GC vs Log2 ------------------------------------------------------
  # 3. GC vs Log2 (MIDDLE LEFT)
  p$gc_log2 <- ggplot(targets) +
    xlab("Target GC content") +
    ylab("Corrected depth") +
    xlim(c(.16, .84)) +
    # Layer 1: Background only
    geom_point(
      mapping = aes(x = gc, y = 2^log2),
      col = pointcolor, shape = 21, size = 1, alpha = alpha / 2
    ) +
    facet_wrap(facets = vars(label), nrow = 1) +
    theme(
      panel.spacing = unit(0, "lines"), strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    scale_y_log10(limits = ylims, breaks = ybreaks, minor_breaks = yminorbreaks)

  # 8. Plot: Order vs Log2 ---------------------------------------------------
  # 4. Order vs Log2 (MIDDLE RIGHT)
  p$order_log2 <- ggplot(targets) +
    xlab("Order of genomic position") +
    ylab("Corrected depth") +
    # Layer 0: Vertical lines
    geom_vline(
      data = cancergenes_summary,
      mapping = aes(xintercept = bin, color = selected_genes),
      alpha = 0.5, linewidth = 0.5, show.legend = FALSE
    ) +
    # Layer 1: Background
    geom_point(
      data = targets[!is.na(label) & is.na(selected_genes)],
      mapping = aes(x = bin, y = 2^log2),
      fill = "#606060", col = "#202020", shape = 21, size = size, alpha = alpha
    ) +
    # Layer 2: Selected
    geom_point(
      data = targets[!is.na(label) & !is.na(selected_genes)],
      mapping = aes(x = bin, y = 2^log2, fill = selected_genes),
      col = "black", shape = 21, size = size_selected, alpha = 0.8
    ) +

    # Layer 3: Segments
    geom_segment(
      data = segments, col = "green", linewidth = 1,
      mapping = aes(x = start, xend = end, y = 2^mean, yend = 2^mean)
    ) +
    # Layer 4: Labels
    geom_text(
      data = cancergenes_summary,
      mapping = aes(x = bin, y = Inf, label = selected_genes, color = selected_genes),
      angle = 45, hjust = 0, size = 2.5, vjust = 0, show.legend = FALSE
    ) +
    scale_fill_hue(l = 70) +
    scale_color_hue(l = 70) +
    scale_y_log10(limits = ylims, breaks = ybreaks, minor_breaks = yminorbreaks) +
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


  # 9. Plot: Position vs Log2 ------------------------------------------------
  # 5. Genomic Position vs Log2 (BOTTOM - full width)
  p$pos_log2 <- ggplot(targets) +
    xlab("Genomic position") +
    ylab("Corrected depth") +
    # Layer 0: Vertical lines
    geom_vline(
      data = cancergenes_summary,
      mapping = aes(xintercept = gpos, color = selected_genes),
      alpha = 0.5, linewidth = 0.5, show.legend = FALSE
    ) +
    # Layer 1: Background points (grey)
    geom_point(
      data = targets[!is.na(label) & is.na(selected_genes)],
      mapping = aes(x = gpos, y = 2^log2),
      fill = "#606060", col = "#202020", shape = 21, size = size, alpha = alpha
    ) +
    # Layer 2: Selected cancer genes (colored by gene)
    geom_point(
      data = targets[!is.na(selected_genes)],
      mapping = aes(x = gpos, y = 2^log2, fill = selected_genes),
      shape = 21, size = size_selected, alpha = 0.8, col = "black"
    ) +

    # Layer 3: Segments
    geom_segment(
      data = segments, col = "green", linewidth = 1,
      mapping = aes(x = gstart, xend = gstop, y = 2^mean, yend = 2^mean)
    ) +
    # Layer 4: Labels
    geom_text(
      data = cancergenes_summary,
      mapping = aes(x = gpos, y = Inf, label = selected_genes, color = selected_genes),
      angle = 45, hjust = 0, size = 2.5, vjust = 0, show.legend = FALSE
    ) +
    scale_fill_hue(l = 70, na.value = "grey50") +
    scale_color_hue(l = 70) +
    scale_y_log10(limits = ylims, breaks = ybreaks, minor_breaks = yminorbreaks) +
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

  # 10. SNP Plots ------------------------------------------------------------
  if (!is.null(snp_table)) {
    # Prepare SNP data
    # defaults to using raw allele ratio:
    snp_table[, allele_ratio_use := allele_ratio]
    # but if there is a corrected allele ratio, use it (not implemented yet but
    # keeping logic):
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
    targets[, allele_ratio := as.double(NA)]
    targets[, DP := as.double(NA)] # Add DP
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

    # 5. Depth vs Allele Ratio (E)
    # Calculate quantiles for depth-based limits (using DP)
    valid_dp <- targets[!is.na(DP)]$DP
    if (length(valid_dp) > 0) {
      m <- quantile(valid_dp, c(.01, .50, .99), na.rm = TRUE)
    } else {
      # Fallback to counts if DP is missing (should verify why)
      m <- quantile(targets[is_target == TRUE]$count, c(.01, .50, .99), na.rm = TRUE)
    }
    
    if (m[1] > m[2] / 3) m[1] <- m[2] / 3
    if (m[3] < m[2] * 3) m[3] <- m[2] * 3
    m[1] <- max(1, m[1])

    # Robustness checks
    if (!is.finite(m[1]) || !is.finite(m[3]) || m[1] <= 0 || m[3] <= 0) {
       m <- c(1, 10, 100)
    }

    p$depth_alleleratio <- ggplot(targets[type != "background"]) +
      xlab("Depth") + # Changed from Fragments
      ylab("Allele ratio") +
      # Layer 1: Background (grey)
      geom_point(
        data = targets[type != "background" & is.na(selected_genes)],
        mapping = aes(x = DP, y = allele_ratio), # Changed count to DP
        fill = "#606060", col = "#202020",
        shape = 21, size = 1, alpha = alpha / 2
      ) +
      # Layer 2: Selected Genes (colored)
      geom_point(
        data = targets[type != "background" & !is.na(selected_genes)],
        mapping = aes(x = DP, y = allele_ratio, fill = selected_genes), # Changed count to DP
        col = "black", shape = 21, size = size_selected, alpha = 0.8
      ) +
      # Layer 2.5: Somatic Mutations (Red) - AF vs Depth
      (if (!is.null(somatic_plot) && nrow(somatic_plot) > 0) {
        geom_point(
          data = somatic_plot,
          mapping = aes(x = DP, y = AF), # Changed count to DP
          col = "red", size = size_selected, alpha = 0.8, shape = 16
        )
      }) +
      (if (!is.null(somatic_labels) && nrow(somatic_labels) > 0) {
        geom_label_repel(
          data = somatic_labels,
          mapping = aes(x = DP, y = AF, label = label), # Changed count to DP
          fill = "white", col = "darkred", segment.color = "red",
          nudge_y = somatic_labels$nudge,
          size = 2.5, show.legend = FALSE, min.segment.length = 0,
          box.padding = 0.1, label.padding = 0.1, point.padding = 0.1
        )
      }) +

      theme(
        panel.spacing = unit(0, "lines"), strip.background = element_blank(),
        strip.text.x = element_blank()
      ) +
      scale_fill_hue(l = 70) +
      scale_x_log10(limits = c(m[1], m[3])) +
      scale_y_continuous(limits = c(-.05, 1.05), breaks = c(0, .25, .5, .75, 1))

    # snp (grid) smooth-to-allele-ratio plot
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
        data = targets[, .(smooth_log2, maf)],
        aes(x = 2^smooth_log2, y = maf), col = "lightgrey", alpha = .2
      ) +
      geom_text(
        data = unique(targets[, .(chromosome)]),
        mapping = aes(x = .25, y = .95, label = chromosome)
      ) +
      # Background layer
      geom_point(
        data = targets[label != "background" & is.na(selected_genes)],
        aes(x = 2^smooth_log2, y = maf),
        fill = "#606060", col = "#202020", shape = 21, size = size, alpha = alpha
      ) +
      # Selected Genes layer
      geom_point(
        data = targets[label != "background" & !is.na(selected_genes)],
        aes(x = 2^smooth_log2, y = maf, fill = selected_genes),
        shape = 21, size = size_selected, alpha = 0.8, col = "black"
      ) +

      facet_wrap(facets = vars(factor(chromosome,
        levels = unique(chromosome),
        ordered = TRUE
      )), ncol = 8) +
      scale_fill_hue(l = 70) +
      theme(
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank()
      )

    # snp (all) smooth-to-allele-ratio plot
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
        shape = 21, size = size_selected, alpha = 0.8, col = "black"
      )


    # allele ratio by order
    p$order_alleleratio <- ggplot(targets) +
      xlab("Order of genomic position") +
      ylab("Allele ratio") +
      # Layer 0: Vertical lines
      geom_vline(
        data = cancergenes_summary,
        mapping = aes(xintercept = bin, color = selected_genes),
        alpha = 0.5, linewidth = 0.5, show.legend = FALSE
      ) +
      # Background
      geom_point(
        data = targets[!is.na(label) & is.na(selected_genes)],
        mapping = aes(x = bin, y = allele_ratio),
        fill = "#606060", col = "#202020", shape = 21, size = size, alpha = alpha
      ) +
      # Selected
      geom_point(
        data = targets[!is.na(label) & !is.na(selected_genes)],
        mapping = aes(x = bin, y = allele_ratio, fill = selected_genes),
        col = "black", shape = 21, size = size_selected, alpha = 0.8
      ) +
      # Layer 2.5: Somatic Mutations (Red)
      (if (!is.null(somatic_plot) && nrow(somatic_plot) > 0) {
        geom_point(
          data = somatic_plot,
          mapping = aes(x = bin, y = AF), # Plot raw AF here (0-1)
          col = "red", size = size_selected, alpha = 0.8, shape = 16
        )
      }) +
      (if (!is.null(somatic_labels) && nrow(somatic_labels) > 0) {
        geom_label_repel(
            data = somatic_labels,
            mapping = aes(x = bin, y = AF, label = label),
            fill = "white", col = "darkred", segment.color = "red",
            nudge_y = somatic_labels$nudge,
            size = 2.5, show.legend = FALSE, min.segment.length = 0,
            box.padding = 0.1, label.padding = 0.1, point.padding = 0.1
        )
      }) +
      # Labels
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

    # Remove redundant quantile calculation here as it is done for depth_alleleratio
    # But for safety, keep m definition if depth_alleleratio wasn't called?
    # Actually depth_alleleratio code block is always run (it's not conditional on snp_table being present in my previous edit?
    # Wait, in original script, depth_alleleratio logic WAS conditional on snp_allele_ratio.


    # 7. GIS Plot (replaces depth_alleleratio in old layout, now moved to K)
    if (!is.null(gis_table)) {
      p$gis_plot <- ggplot2::ggplot() +
        ggplot2::geom_line(data = gis_table, aes(x = fraction, y = predicted_gis, group = 1), linewidth = 2) +
        ggplot2::geom_line(data = gis_table, aes(x = fraction, y = local_cnv * 10, col = "local_cnv"), lty = 1, linewidth = 1) +
        ggplot2::geom_line(data = gis_table, aes(x = fraction, y = focal_gain * 10, col = "focal_gain"), lty = 1, linewidth = 1) +
        ggplot2::geom_line(data = gis_table, aes(x = fraction, y = loh * 10, col = "loh"), lty = 1, linewidth = 1) +
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
        labs(color = "Feature (x10)") +
        theme(legend.position = "bottom")
    } else {
      p$gis_plot <- ggplot2::ggplot() +
        ggplot2::geom_text(aes(x = 0.5, y = 0.5, label = "No GIS Data")) +
        theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
    }


    # allele ratio by genomic position
    # 10. Plot: Position vs Allele Ratio ---------------------------------------
    p$pos_alleleratio <- ggplot(targets) +
      xlab("Genomic position") +
      ylab("Allele ratio") +
      # Layer 0: Vertical lines
      geom_vline(
        data = cancergenes_summary,
        mapping = aes(xintercept = gpos, color = selected_genes),
        alpha = 0.5, linewidth = 0.5, show.legend = FALSE
      ) +
      # Layer 1: Background
      geom_point(
        data = targets[!is.na(label) & is.na(selected_genes)],
        mapping = aes(x = gpos, y = allele_ratio),
        fill = "#606060", col = "#202020", shape = 21, size = size, alpha = alpha
      ) +
      # Layer 2: Selected
      geom_point(
        data = targets[!is.na(label) & !is.na(selected_genes)],
        mapping = aes(x = gpos, y = allele_ratio, fill = selected_genes),
        col = "black", shape = 21, size = size_selected, alpha = 0.8
      ) +
      # Layer 2.5: Somatic Mutations (Red)
      (if (!is.null(somatic_plot) && nrow(somatic_plot) > 0) {
        geom_point(
          data = somatic_plot,
          mapping = aes(x = gpos, y = AF),
          col = "red", size = 1.2, alpha = 0.8, shape = 16
        )
      }) + 
      # Layer 2.6: Labels
      (if (!is.null(somatic_labels) && nrow(somatic_labels) > 0) {
        geom_label_repel(
            data = somatic_labels,
            mapping = aes(x = gpos, y = AF, label = label),
            fill = "white", col = "darkred", segment.color = "red",
            nudge_y = somatic_labels$nudge,
            size = 2.5, show.legend = FALSE, min.segment.length = 0,
            box.padding = 0.1, label.padding = 0.1, point.padding = 0.1
        )
      }) +
      # Layer 3: Labels
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
  }

  # Apply legend formatting
  for (i in 1:length(p)) {
    p[[i]] <- p[[i]] + guides(fill = guide_legend(override.aes = list(
      shape = 21, size = 3
    ))) + labs(fill = NULL)
  }

  # 11. Layout & Save --------------------------------------------------------
  # We use the 'patchwork' package to arrange the individual plots into a
  # composite figure.
  # The layout string defines the grid, where letters correspond to the plot
  # objects below.
  if (!is.null(snp_table)) {
    layout <- "AABBBBBB
                CCDDDDDD
                EEFFFFFF
                GGHHHHHH"
    fig <- p$gc_rawdepth + p$order_rawdepth +
      p$gc_log2 + p$order_log2 +
      p$depth_alleleratio + p$order_alleleratio +
      # p$pos_log2 +
      # p$pos_alleleratio +
      p$nogrid + p$grid +
      plot_layout(design = layout, guides = "collect") &
      theme(legend.position = "none")
  } else {
    layout <- "ABBB
               CDDD
               EEEE"

    fig <- p$gc_rawdepth + p$order_rawdepth +
      p$gc_log2 + p$order_log2 +
      p$pos_log2 +
      plot_layout(design = layout, guides = "collect") &
      theme(legend.position = "none")
  }

  # Add title/caption matching original format
  pa <- plot_annotation(
    title = title,
    caption = paste("Jumble on", format(Sys.time(), "%a %b %e %Y, %H:%M"))
  )
  fig <- fig + pa

  # Save if output file specified
  if (!is.null(output_file)) {
    png(file = output_file, width = 1600, height = 1300, res = 100)
    suppressWarnings(print(fig))
    dev.off()
  }

  invisible(fig)
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
#' @export
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
      plot_data <- data.table::melt(
        gis_table,
        id.vars = "fraction",
        measure.vars = c("predicted_gis", "local_cnv", "focal_gain", "loh"),
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
    geom_line(data = plot_data[feature == "predicted_gis"], aes(x = fraction, y = value), linewidth = 2) +
    geom_line(
      data = plot_data[feature %in% c("local_cnv", "focal_gain", "loh")],
      aes(x = fraction, y = value, col = feature), linetype = 1, linewidth = 1
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
      theme = theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5))
    )

  # Save
  suppressWarnings(ggsave(output_file, plot = final_plot, width = 16, height = 12, dpi = 300))
}
