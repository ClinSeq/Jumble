#' GIS Model
#'
#' Hard-coded linear prediction function for Genomic Instability Score (GIS).
#'
#' @details
#' The GIS score is predicted from three genomic instability features using a
#' linear model with hard-coded coefficients:
#'
#' \deqn{GIS = 2.632 + (-0.852 \times focal\_gain) + (1.799 \times local\_cnv) + (2.727 \times loh)}
#'
#' The coefficients were derived from a training set of samples with known HRD
#' status (Myriad myChoice GIS scores). The model uses only three of the nine
#' available features:
#'
#' \itemize{
#'   \item \strong{loh} (coefficient 2.727): Strongest predictor. Loss of
#'     heterozygosity is a hallmark of homologous recombination deficiency.
#'   \item \strong{local_cnv} (coefficient 1.799): Captures arm-level copy
#'     number changes characteristic of HRD genomes.
#'   \item \strong{focal_gain} (coefficient -0.852): Negative coefficient.
#'     Focal gains are more characteristic of non-HRD tumors (e.g., oncogene
#'     amplification) and thus reduce the predicted GIS.
#' }
#'
#' When called with \code{feats = NULL}, returns the coefficient list for
#' inspection. Users who wish to use a different model (e.g., randomForest
#' trained on their own cohort) can supply it via the \code{hrd_model_file}
#' parameter in \code{\link{run_jumble}}.
#'
#' See also \code{docs/METHODS.md} Section 6.7 for full details.
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
#' @details
#' This function orchestrates the full Genomic Instability Score computation.
#' The algorithm operates through a multi-scale binning pipeline:
#'
#' \strong{1. 1 Mb Binning:} Target-level data is aggregated into 1 Mb bins
#' by taking the median log2 and MAF per bin. Values are then re-averaged by
#' segment ID to use the segmented (smoothed) signal rather than raw bin
#' values.
#'
#' \strong{2. 5 Mb Binning:} A genome-wide 5 Mb bin template is constructed
#' from chromosome arm definitions (centromeric bins excluded). The 1 Mb
#' values are mapped to 5 Mb bins via genomic overlap, taking the median per
#' bin. Bins must have both log2 AND MAF data to be included.
#'
#' \strong{3. Arm-level and long-range statistics:} Per-arm median log2 and
#' a 25 Mb running median (k=5 bins of 5 Mb) are computed as baselines for
#' feature extraction.
#'
#' \strong{4. Tumor fraction sweep:} Nine genomic instability features are
#' computed at each of 100 assumed tumor fractions (0.01 to 1.00), producing
#' a GIS score profile. This allows identification of the most likely tumor
#' fraction and assessment of GIS sensitivity to purity uncertainty.
#'
#' The multi-scale approach ensures that GIS features capture arm-level and
#' sub-arm-level events rather than bin-level noise, making the scoring
#' applicable to both sparse gene panel data and dense WGS data.
#'
#' See also \code{\link{comp_gis_for_fraction}} for per-fraction feature
#' computation and \code{docs/METHODS.md} Section 6 for the full description.
#'
#' @param targets Data.table of target bins (must have segment, log2, chromosome, start, end).
#' @param snps Data.table of SNPs (must have chromosome, start, end, allele_ratio).
#' @param fractions Vector of purity fractions to test (default: 0.01 to 1.00 by 0.01).
#' @param genome Genome version ("hg19" or "hg38"). Default "hg19".
#' @param hrd_model Optional HRD model.
#' @return Data.table with fraction, predicted_gis, and feature counts.
#' @importFrom data.table data.table copy as.data.table setnames := rbindlist
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom stats median quantile runmed
#' @importFrom stringr str_detect
#' @keywords internal
compute_gis_table <- function(targets, snps = NULL, fractions = seq(0.01, 1.00, by = 0.01), genome = "hg19", hrd_model = NULL) {
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
        comp_gis_for_fraction(bins_final, f, hrd_model = hrd_model)
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

#' Compute GIS Features for a Single Tumor Fraction
#'
#' Computes nine genomic instability features and the predicted GIS score
#' for a single assumed tumor fraction.
#'
#' @details
#' All feature thresholds are scaled by the assumed tumor fraction:
#' \code{threshold = max(tumor_fraction / 3, 0.05)} and
#' \code{threshold_log2 = log2(1 + threshold)}. This ensures that at low
#' tumor fractions, thresholds are relaxed to detect diluted signals.
#'
#' \strong{Nine features are computed:}
#'
#' \emph{Copy number features} (counted as number of distinct chromosome arms):
#' \itemize{
#'   \item \strong{transitions}: 5 Mb bins where |Δlog2| > threshold between
#'     adjacent bins (sum of bins, not arms).
#'   \item \strong{long_cnv}: Arms where the 25 Mb running median deviates
#'     from the arm median by > threshold_log2.
#'   \item \strong{local_cnv}: Arms with any bin showing local gain OR loss.
#'   \item \strong{local_gain / local_loss}: Arms with bins deviating from
#'     arm median by > threshold_log2.
#'   \item \strong{focal_gain / focal_loss}: Arms with bins deviating from
#'     BOTH the long-range median AND the arm median. This dual requirement
#'     distinguishes true focal events from arm-level shifts.
#' }
#'
#' \emph{Allele-based features:}
#' \itemize{
#'   \item \strong{loh}: Arms with > 3 bins showing LOH AND > 3 bins without
#'     LOH (partial LOH). LOH is called using a purity-aware threshold:
#'     \code{copyratio = 1 + (2^log2 - 1) / tfr}, then
#'     \code{dnaratio = copyratio * tfr / (copyratio * tfr + 1 * (1 - tfr))},
#'     and \code{thr_maf = max(0.5 + dnaratio * 0.8 * 0.5, 0.55)}.
#'     A bin is LOH if MAF > thr_maf. LOH calls are smoothed with a running
#'     median (k=5 = 25 Mb) per arm.
#'   \item \strong{tai}: Arms with allelic imbalance (smoothed MAF > 0.6) at
#'     a telomeric bin.
#' }
#'
#' See \code{docs/METHODS.md} Sections 6.3–6.6 for mathematical details.
#'
#' @param bins_in Data.table of 5 Mb bins with log2, maf, arm, chromosome,
#'   telomeric, arm_median, and long_median columns.
#' @param tfr Numeric tumor fraction (0 to 1).
#' @param hrd_model Optional custom HRD model object.
#' @return Named list with fraction, 9 feature counts, predicted_gis, and
#'   optionally custom_HRD.
#' @keywords internal
comp_gis_for_fraction <- function(bins_in, tfr, hrd_model = NULL) {
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

    custom_score <- apply_hrd_model(hrd_model, res)
    if (!is.null(custom_score)) {
        res$custom_HRD <- custom_score
    }

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

#' Apply Custom HRD Model
#'
#' Applies a user-supplied HRD model to the GIS feature set for a single
#' fraction. The model receives a data frame with all 9 standard GIS
#' features and selects its own columns by name.
#'
#' For regression models the raw numeric prediction is returned.
#' For classification models the positive-class probability × 100 is returned.
#'
#' @param model Model object, function, or NULL.
#' @param feats List of feature values from comp_gis_for_fraction.
#' @return Numeric score, or NULL if no model supplied.
#' @keywords internal
apply_hrd_model <- function(model, feats) {
    if (is.null(model)) return(NULL)

    # Check randomForest availability when model requires it
    if (inherits(model, "randomForest") &&
        !require("randomForest", character.only = TRUE, quietly = TRUE)) {
        stop("Package 'randomForest' is required to use the supplied hrd_model. ",
             "Install it with: install.packages('randomForest')")
    }

    # Build full feature data frame — model picks its columns by name
    feats_df <- as.data.frame(feats[c(
        "focal_gain", "focal_loss", "local_cnv",
        "local_gain", "local_loss", "loh",
        "transitions", "long_cnv", "tai"
    )], check.names = FALSE)

    # Plain function interface
    if (is.function(model)) {
        return(as.numeric(model(feats_df)))
    }

    # Try classification interface (type = "prob")
    pred <- tryCatch(
        predict(model, newdata = feats_df, type = "prob"),
        error = function(e) NULL
    )

    if (!is.null(pred) && (is.data.frame(pred) || is.matrix(pred))) {
        col_names <- colnames(pred)
        priority  <- c("GIS", "HRD", "Pos", "1")
        match_idx <- which(toupper(col_names) %in% toupper(priority))

        if (length(match_idx) > 0) {
            best <- match_idx[which.min(match(
                toupper(col_names[match_idx]), toupper(priority)
            ))]
            prob_col <- as.numeric(pred[, best])
        } else {
            prob_col <- as.numeric(pred[, min(2, ncol(pred))])
        }
        return(prob_col * 100)
    }

    # Regression fallback
    as.numeric(predict(model, newdata = feats_df))
}
