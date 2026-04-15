#' GIS Model
#'
#' Hard-coded prediction function for Genomic Instability Score (GIS).
#'
#' @param feats List or data.table containing features: focal_gain, local_cnv, loh.
#' @return Predicted GIS score.
#' @keywords internal
gis_model <- function(feats = NULL) {
    # 1. Coefficients --------------------------------------------------------
    coefficients <- list(
        `(Intercept)` = 2.63187808497402,
        focal_gain = -0.851910778439913,
        local_cnv = 1.79873839268124,
        loh = 2.72722833695396
    )

    if (is.null(feats)) {
        return(coefficients)
    }

    # 2. Calculate Prediction ------------------------------------------------
    predicted_GIS <-
        coefficients$`(Intercept)` +
        coefficients$focal_gain * feats$focal_gain +
        coefficients$local_cnv * feats$local_cnv +
        coefficients$loh * feats$loh

    return(predicted_GIS)
}

#' Compute GIS Table
#'
#' Computes GIS scores across a range of purity fractions.
#'
#' @param targets Data.table of target bins (must have segment, log2, chromosome, start, end).
#' @param snps Data.table of SNPs (must have chromosome, start, end, allele_ratio).
#' @param fractions Vector of purity fractions to test (default: 0.01 to 1.00 by 0.01).
#' @param genome Genome version ("hg19" or "hg38"). Default "hg19".
#' @return Data.table with fraction, predicted_gis, and feature counts.
#' @importFrom data.table data.table copy as.data.table setnames := rbindlist
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom stats median quantile runmed
#' @importFrom stringr str_detect
#' @keywords internal
compute_gis_table <- function(targets, snps = NULL, fractions = seq(0.01, 1.00, by = 0.01), genome = "hg19") {
    # 1. 1MB Binning of Targets ----------------------------------------------
    # parsing targets
    targets_sub <- targets[chromosome %in% c(1:22, "X", "Y") & !is.na(log2)]

    # parsing snps (if available) - map to targets first
    if (!is.null(snps) && nrow(snps) > 0) {
        snps_sub <- snps[chromosome %in% c(1:22, "X", "Y")][allele_ratio > .02 & allele_ratio < .98]
        snps_sub[, maf := 0.5 + abs(allele_ratio - 0.5)]

        # Map MAF to targets
        targetranges <- makeGRangesFromDataFrame(targets_sub,
            seqnames.field = "chromosome",
            start.field = "start",
            end.field = "end",
            ignore.strand = TRUE
        )
        snpranges <- makeGRangesFromDataFrame(snps_sub, keep.extra.columns = TRUE)

        overlap <- as.data.table(findOverlaps(snpranges, targetranges))

        # Initialize maf column in targets
        targets_sub[, maf := as.numeric(NA)]

        if (nrow(overlap) > 0) {
            targets_sub[overlap$subjectHits, maf := snps_sub[overlap$queryHits]$maf]
        }
    } else {
        targets_sub[, maf := as.numeric(NA)]
    }

    # "Let all positions be megabase"
    targets_binned <- copy(targets_sub)
    targets_binned[, start_mb := round((start + end) / 2e6) * 1e6]
    targets_binned[, pos := paste0(chromosome, ":", start_mb / 1e6)]

    # Aggregate per 1MB bin (targets_mb)
    # Median log2, maf, and segment per 1MB bin
    targets_mb <- targets_binned[, .(
        log2 = median(log2, na.rm = TRUE),
        maf = median(maf, na.rm = TRUE),
        segment = if (all(is.na(segment))) as.numeric(NA) else round(median(segment, na.rm = TRUE)),
        start = unique(start_mb),
        end = unique(start_mb)
    ), by = .(chromosome, pos)]

    # "Set segmented log2 + maf value" - re-averaging by segment ID
    # Note: original script calculates 'log2' as median(log2) by segment.
    # This means the 1MB bin value becomes the SEGMENT mean.
    targets_mb[, log2 := median(log2, na.rm = TRUE), by = segment]
    targets_mb[, maf := median(maf, na.rm = TRUE), by = segment]

    targets_mb <- targets_mb[!is.na(log2)]


    # 2. 5MB Binning (Template) ----------------------------------------------
    # Get arm definitions
    chroms <- get_chrom_arms(genome)

    binlength <- 5e6 # 5MB

    binning_list <- list()
    for (i in 1:nrow(chroms)) {
        bins <- data.table(
            chromosome = chroms[i]$chromosome,
            arm = chroms[i]$arm,
            start = seq(chroms[i]$start, chroms[i]$end - binlength, binlength)
        )
        bins[, end := start + binlength]
        bins[, telomeric := 0][, centromeric := 0]

        # First and last bin of arm are telomeric/centromeric boundaries,
        # but for 'telomeric' feature only actual telomeres matter.
        # Script logic (matched to frankenscript.R):
        # if p-arm: bin[1] is telomeric, bin[.N] is centromeric
        # if q-arm: bin[.N] is telomeric, bin[1] is centromeric
        if (str_detect(chroms[i]$arm, "p")) bins[1, telomeric := 1][.N, centromeric := 1]
        if (str_detect(chroms[i]$arm, "q")) bins[.N, telomeric := 1][1, centromeric := 1]
        binning_list[[i]] <- bins
    }
    binning <- rbindlist(binning_list)[centromeric == 0]

    # 3. Map 1MB -> 5MB ------------------------------------------------------
    # Note: targets_mb has 'start' and 'end' as start_mb (point).
    # makeGRangesFromDataFrame needs width.
    targets_mb[, end := start] # Point ranges
    tr <- makeGRangesFromDataFrame(targets_mb,
        seqnames.field = "chromosome",
        start.field = "start", end.field = "end", ignore.strand = TRUE
    )

    br <- makeGRangesFromDataFrame(binning,
        seqnames.field = "chromosome",
        start.field = "start", end.field = "end", ignore.strand = TRUE
    )

    overlap <- as.data.table(findOverlaps(br, tr))

    # Aggregate values to bins
    # For each bin (queryHit), get median of mapped targets (subjectHit)
    # Using data.table join for speed

    # Add indices
    binning[, bid := 1:.N]
    targets_mb[, tid := 1:.N]

    # Join overlap with targets
    ov_vals <- merge(overlap, targets_mb[, .(tid, log2, maf)], by.x = "subjectHits", by.y = "tid")

    # Compute medians per bin
    bin_stats <- ov_vals[, .(
        log2 = median(log2, na.rm = TRUE),
        maf = median(maf, na.rm = TRUE)
    ), by = .(queryHits)]

    # Merge back to binning
    bins_final <- merge(binning, bin_stats, by.x = "bid", by.y = "queryHits", all.x = TRUE)

    # Filter valid bins (Must have both log2 AND maf, per frankenscript.R)
    bins_final <- bins_final[!is.na(log2) & !is.na(maf)]

    # 4. Compute Features ----------------------------------------------------
    bins_final[, arm_median := median(log2, na.rm = TRUE), by = arm]

    # Long median (smoothed per chromosome, k=5 -> 25MB)
    bins_final[, long_median := safe_runmed(log2, k = 5), by = chromosome]

    results_list <- lapply(fractions, function(f) {
        comp_gis_for_fraction(bins_final, f)
    })

    results_dt <- rbindlist(results_list)
    attr(results_dt, "bins_final") <- bins_final

    return(results_dt)
}

#' Get Chromosome Arms
#' @keywords internal
get_chrom_arms <- function(genome = "hg19") {
    if (genome == "hg19") {
        chroms <- data.table(
            chromosome = c(
                "1", "1", "2", "2", "3", "3", "4", "4", "5", "5", "6", "6", "7", "7", "8", "8", "9", "9", "10", "10",
                "11", "11", "12", "12", "13", "13", "14", "14", "15", "15", "16", "16", "17", "17", "18", "18", "19", "19", "20", "20",
                "21", "21", "22", "22", "X", "X"
            ),
            arm = c(
                "1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q", "5p", "5q", "6p", "6q", "7p", "7q", "8p", "8q", "9p", "9q", "10p", "10q",
                "11p", "11q", "12p", "12q", "13p", "13q", "14p", "14q", "15p", "15q", "16p", "16q", "17p", "17q", "18p", "18q", "19p", "19q", "20p", "20q",
                "21p", "21q", "22p", "22q", "Xp", "Xq"
            ),
            start = c(
                0, 1.25e+08, 0, 93300000, 0, 9.1e+07, 0, 50400000, 0, 48400000, 0, 6.1e+07, 0, 59900000, 0, 45600000, 0, 4.9e+07, 0, 40200000,
                0, 53700000, 0, 35800000, 0, 17900000, 0, 17600000, 0, 19000000, 0, 36600000, 0, 24000000, 0, 17200000, 0, 26500000, 0, 27500000,
                0, 13200000, 0, 14700000, 0, 60600000
            ),
            end = c(
                1.25e+08, 249250621, 93300000, 243199373, 9.1e+07, 198022430, 50400000, 191154276, 48400000, 180915260, 6.1e+07, 171115067, 59900000, 159138663, 45600000, 146364022, 4.9e+07, 141213431, 40200000, 135534747,
                53700000, 135006516, 35800000, 133851895, 17900000, 115169878, 17600000, 107349540, 19000000, 102531392, 36600000, 90354753, 24000000, 81195210, 17200000, 78077248, 26500000, 59128983, 27500000, 63025520,
                13200000, 48129895, 14700000, 51304566, 60600000, 155270560
            )
        )
    } else if (genome == "hg38") {
        # Approx hg38 centromeres (splits)
        chroms <- data.table(
            chromosome = c(
                "1", "1", "2", "2", "3", "3", "4", "4", "5", "5", "6", "6", "7", "7", "8", "8", "9", "9", "10", "10",
                "11", "11", "12", "12", "13", "13", "14", "14", "15", "15", "16", "16", "17", "17", "18", "18", "19", "19", "20", "20",
                "21", "21", "22", "22", "X", "X"
            ),
            arm = c(
                "1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q", "5p", "5q", "6p", "6q", "7p", "7q", "8p", "8q", "9p", "9q", "10p", "10q",
                "11p", "11q", "12p", "12q", "13p", "13q", "14p", "14q", "15p", "15q", "16p", "16q", "17p", "17q", "18p", "18q", "19p", "19q", "20p", "20q",
                "21p", "21q", "22p", "22q", "Xp", "Xq"
            ),
            start = c(
                0, 123.4e6, 0, 93.9e6, 0, 90.9e6, 0, 49.8e6, 0, 48.8e6, 0, 58.6e6, 0, 60.1e6, 0, 45.2e6, 0, 43e6, 0, 39.8e6,
                0, 53.4e6, 0, 35.5e6, 0, 17.7e6, 0, 17.2e6, 0, 19e6, 0, 36.8e6, 0, 25.1e6, 0, 18.5e6, 0, 26.2e6, 0, 28.1e6,
                0, 12e6, 0, 15e6, 0, 61e6 # X approx
            ),
            end = c(
                123.4e6, 248956422, 93.9e6, 242193529, 90.9e6, 198295559, 49.8e6, 190214555, 48.8e6, 181538259, 58.6e6, 170805979, 60.1e6, 159345973, 45.2e6, 145138636, 43e6, 138394717, 39.8e6, 133797422,
                53.4e6, 135086622, 35.5e6, 133275309, 17.7e6, 114364328, 17.2e6, 107043718, 19e6, 101991189, 36.8e6, 90338345, 25.1e6, 83257441, 18.5e6, 80373285, 26.2e6, 58617616, 28.1e6, 64444167,
                12e6, 46709983, 15e6, 50818468, 61e6, 156040895
            )
        )
        # Remove Y from definition if not in hg19 list?
        # hg19 list has X but not Y.
        # Let's keep it consistent.
        if ("Y" %in% chroms$chromosome) chroms <- chroms[chromosome != "Y"]
    } else {
        stop("Unsupported genome: ", genome)
    }
    return(chroms)
}

#' Helper to compute GIS for a single fraction
#' @keywords internal
comp_gis_for_fraction <- function(bins_in, tfr) {
    # Work on a copy
    bins <- copy(bins_in)

    # 1. Transitions ---------------------------------------------------------
    threshold <- ifelse(is.na(tfr), 0.05, tfr / 3)
    if (threshold < 0.05) threshold <- 0.05
    threshold_log2 <- log2(1 + threshold)

    # Copy number transitions
    # func: abs(diff(data)) > threshold
    bins[, transitions_val := 0]

    # diff based on adjacent bins in the 5MB template
    calc_trans <- function(x, thr) {
        if (length(x) < 2) {
            return(rep(0, length(x)))
        }
        d <- which(abs(diff(x)) > thr)
        res <- rep(0, length(x))
        res[d] <- 1
        res[d + 1] <- 1
        return(res)
    }

    bins[, transitions_val := calc_trans(log2, threshold_log2), by = chromosome]

    # 2. LOH -----------------------------------------------------------------
    calc_loh <- function(log2_vec, maf_vec, tfr_val) {
        res <- rep(0, length(log2_vec))
        if (all(is.na(maf_vec))) {
            return(res)
        }

        copyratio <- 1 + 1 * ((2^log2_vec - 1) / tfr_val)
        dnaratio <- copyratio * tfr_val / (copyratio * tfr_val + 1 * (1 - tfr_val))

        thr <- dnaratio * 0.8
        thr_maf <- 0.5 + thr * 0.5
        thr_maf[thr_maf < 0.55] <- 0.55

        # Allow NA maf to not trigger LOH
        valid <- !is.na(maf_vec)
        is_loh <- which(valid & maf_vec > thr_maf)
        res[is_loh] <- 1
        return(res)
    }

    bins[, loh_raw := calc_loh(log2, maf, tfr)]

    # Smooth LOH (k=5 bins = 25MB)
    bins[, loh_smooth := safe_runmed(loh_raw, 5), by = arm]

    # Arm stats
    bins[, arm_loh := sum(loh_smooth == 1), by = arm]
    bins[, arm_noloh := sum(loh_smooth == 0), by = arm]

    # 3. Allelic Imbalance ---------------------------------------------------
    bins[, imbalance := 0]
    bins[, imbalance := {
        m_clean <- maf
        m_clean[is.na(m_clean)] <- 0.5
        as.integer(safe_runmed(m_clean, 5) > 0.6)
    }, by = arm]

    # 4. Features ------------------------------------------------------------
    # Note: long_median and arm_median already computed in parent function
    bins[, long_cnv := as.integer(abs(long_median - arm_median) > threshold_log2)]

    bins[, local_gain := as.integer(log2 - arm_median > threshold_log2)]
    bins[, local_loss := as.integer(log2 - arm_median < -threshold_log2)]

    bins[, focal_gain := as.integer(log2 - long_median > threshold_log2 & log2 - arm_median > threshold_log2)]
    bins[, focal_loss := as.integer(log2 - long_median < -threshold_log2 & log2 - arm_median < -threshold_log2)]

    # Aggregating results
    res <- list()
    res$fraction <- tfr
    res$transitions <- round(sum(bins$transitions_val))
    res$long_cnv <- length(unique(bins[long_cnv == 1]$arm))
    res$local_cnv <- length(unique(bins[local_gain == 1 | local_loss == 1]$arm))
    res$local_gain <- length(unique(bins[local_gain == 1]$arm))
    res$local_loss <- length(unique(bins[local_loss == 1]$arm))
    res$focal_gain <- length(unique(bins[focal_gain == 1]$arm))
    res$focal_loss <- length(unique(bins[focal_loss == 1]$arm))

    # loh: length(unique(bins[arm_loh>3 & arm_noloh>3]$arm))
    res$loh <- length(unique(bins[arm_loh > 3 & arm_noloh > 3]$arm))

    # tai: telomeric imbalance
    # uses 'telomeric' column set in parent
    res$tai <- length(unique(bins[imbalance == 1 & telomeric == 1]$arm))

    # 5. Predict GIS ---------------------------------------------------------
    res$predicted_gis <- round(gis_model(res))

    return(res)
}

#' Safe runmed helper
#' @keywords internal
safe_runmed <- function(x, k) {
    if (length(x) < k) {
        return(median(x, na.rm = TRUE))
    }
    stats::runmed(x, k)
}
