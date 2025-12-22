#' @title Label Cells Based on Binary Prediction and Significance
#' @description
#' Assigns categorical labels (\code{"Positive"}, \code{"Negative"}, or \code{"Neutral"})
#' to cells using a precomputed binary prediction (0/1) and FDR-adjusted significance.
#' Only cells passing the FDR threshold (\code{PIPET_FDR < 0.05}) receive non-neutral labels.
#'
#' Designed as a companion to probabilistic or classifier-based annotation pipelines
#' (e.g., logistic regression, SVM, or random forest outputs) where \code{PIPET_prediction}
#' encodes the hard class assignment.
#'
#' @param res_df A data frame containing at least:
#'   \itemize{
#'     \item \code{PIPET_prediction}: numeric (0 or 1), where 1 = predicted positive class.
#'     \item \code{PIPET_FDR}: FDR-adjusted p-values for the prediction confidence.
#'   }
#'
#' @return A copy of \code{res_df} with a new column \code{PIPET} (factor-like character vector)
#'   containing the final labels:
#'   \describe{
#'     \item{\code{"Positive"}}{Predicted class = 1 AND FDR < 0.05}
#'     \item{\code{"Negative"}}{Predicted class = 0 AND FDR < 0.05}
#'     \item{\code{"Neutral"}}{Otherwise (non-significant or prediction unclear)}
#'   }
#'
#' @note The FDR threshold is currently hardcoded at 0.05. To support configurable thresholds,
#'   consider adding an \code{alpha} argument in future versions.
#'
#' @examples
#' df <- data.frame(
#'   PIPET_prediction = c(1, 1, 0, 0, 1),
#'   PIPET_FDR = c(0.01, 0.08, 0.03, 0.1, 0.04)
#' )
#' BinaryLabelCell(df)
#' #>   PIPET_prediction PIPET_FDR   PIPET
#' #> 1                1      0.01 Positive
#' #> 2                1      0.08 Neutral
#' #> 3                0      0.03 Negative
#' #> 4                0      0.10 Neutral
#' #> 5                1      0.04 Positive
#'
#' @export
BinaryLabelCell <- function(res_df) {
  res_df$PIPET <- "Neutral"
  res_df$PIPET[
    res_df$PIPET_prediction == 1 & res_df$PIPET_FDR < 0.05
  ] <- "Positive"
  res_df$PIPET[
    res_df$PIPET_prediction == 0 & res_df$PIPET_FDR < 0.05
  ] <- "Negative"

  res_df
}

#' @title Label Cells Using Continuous Phenotype Scoring from PIPET
#' @description
#' Assigns each cell a categorical label (`"Positive"`, `"Negative"`, or `"Neutral"`)
#' based on a continuous phenotype score derived from weighted barycentric projection
#' onto reference group centroids (i.e. PIPET score). The score is smoothed via kernel density estimation (KDE),
#' and thresholds are adaptively inferred from bimodality (valleys) or spread statistics.
#' Significance (via FDR-corrected p-values) gates confident assignments.
#'
#' Intended for outputs from proximity-inference tools (e.g., PIPET-style pipelines),
#' where multiple distance columns encode dissimilarity to pre-defined reference groups.
#'
#' @param res_df A data frame containing:
#'   \itemize{
#'     \item Distance columns named like \code{paste0(group_prefix, <id>)}, e.g.
#'           \code{"PIPET_dist_1"}, \code{"PIPET_dist_2"}, ...
#'     \item A p-value column (e.g., \code{"PIPET_Pvalue"})
#'     \item An FDR column (e.g., \code{"PIPET_FDR"})
#'   }
#' @param group_prefix Character string prefix used to identify distance columns.
#'   Default: \code{"PIPET_dist_"}.
#' @param pval_col Name of the column containing raw p-values (used only for significance flagging).
#' @param fdr_col Name of the column containing FDR-adjusted p-values. Cells with
#'   \code{res_df\[\[fdr_col\]\] <= alpha} are considered statistically significant.
#' @param alpha Significance threshold for FDR. Default: \code{0.05}.
#' @param lambda Smoothing parameter for distance-to-weight conversion:
#'   \eqn{w_{ik} = \exp(-\lambda \cdot d_{ik})}. Controls how sharply distances decay
#'   into weights. If \code{NULL} (default), \code{lambda} is auto-tuned so that the
#'   median ratio of best-to-second-best weight is ~5 (empirically balances specificity
#'   and robustness). Clamped to \code{[2, 50]}.
#'
#' @return A copy of \code{res_df} with three new columns added:
#'   \describe{
#'     \item{\code{PIPET_phenotype_score}}{Continuous barycentric score \eqn{s_i \in [\min(k), \max(k)]},
#'       where higher values indicate stronger alignment with higher-indexed groups.}
#'     \item{\code{PIPET_match_confidence}}{Per-cell confidence:
#'       \eqn{\mathrm{plogis}(\lambda \cdot \Delta d_i) \cdot \mathbb{I}(FDR \le \alpha)},
#'       where \eqn{\Delta d_i} = gap between 1st and 2nd smallest distances. Ranges in \eqn{[0, 1]}.}
#'     \item{\code{PIPET}}{Categorical label: \code{"Positive"}, \code{"Negative"}, or \code{"Neutral"}.}
#'   }
#'
#' @section Labeling Logic:
#' \itemize{
#'   \item \strong{Positive}: Significant (\code{FDR <= alpha}) AND score > \code{thresh_high}
#'   \item \strong{Negative}: Significant AND score < \code{thresh_low}
#'   \item \strong{Neutral}: Otherwise (non-significant or intermediate score)
#' }
#'
#' Thresholds \code{thresh_low} / \code{thresh_high} are inferred hierarchically:
#' \enumerate{
#'   \item If KDE of significant scores is clearly bimodal (\eqn{\ge 2} valleys),
#'         use leftmost and rightmost valley positions.
#'   \item If unimodal (1 valley), use \eqn{\text{valley} \pm \mathrm{MAD}}.
#'   \item If no valleys detected, fall back to 25th/75th percentiles of significant scores.
#' }
#'
#'
#' @examples
#' # Simulate 3 reference groups (e.g., anterior/mid/posterior)
#' set.seed(42)
#' n_cells <- 200
#' dist_mat <- matrix(rexp(n_cells * 3, rate = 0.5), ncol = 3)
#' colnames(dist_mat) <- c("PIPET_dist_1", "PIPET_dist_2", "PIPET_dist_3")
#' pvals <- runif(n_cells, 0, 0.1)
#' fdrs <- p.adjust(pvals, "BH")
#' res_df <- data.frame(dist_mat, PIPET_Pvalue = pvals, PIPET_FDR = fdrs)
#'
#' labeled <- ContinuousLabelCell(res_df)
#' table(labeled$PIPET)
#' head(labeled[, c("PIPET_phenotype_score", "PIPET_match_confidence", "PIPET")])
#'
#' @export
ContinuousLabelCell <- function(
  res_df,
  group_prefix = "PIPET_dist_",
  pval_col = "PIPET_Pvalue",
  fdr_col = "PIPET_FDR",
  alpha = 0.05,
  lambda = NULL
) {
  dist_cols <- grep(paste0("^", group_prefix), names(res_df), value = TRUE)
  # value tag
  group_ids <- as.numeric(gsub(paste0(".*", group_prefix), "", dist_cols))

  ord <- order(group_ids)
  dist_cols <- dist_cols[ord]
  group_ids <- group_ids[ord]
  K <- length(group_ids)

  #  s = weighted barycenter
  dist_mat <- as.matrix(res_df[dist_cols]) # n × K
  n <- nrow(dist_mat)

  # let median(max_weight / second_max_weight) ≈ 5
  if (is.null(lambda)) {
    margins <- apply(dist_mat, 1, function(d) {
      ord_d <- order(d)
      d[ord_d[1]] - d[ord_d[2]] # best - second
    })
    # lambda ≈ log(5) / median(margin) for ratio=5
    lambda <- ifelse(
      stats::median(margins) > 0,
      log(5) / stats::median(margins),
      10
    )
    lambda <- min(max(lambda, 2), 50) # clamp
  }

  weights <- exp(-lambda * dist_mat) # n × K, larger weight for smaller dist
  s <- rowSums(sweep(weights, 2, group_ids, `*`)) / rowSums(weights) # weighted average of z_k = k

  # Weighted KDE on s
  # compute margin & weight_raw always so later code can use it
  margin <- apply(dist_mat, 1, function(d) {
    d_sorted <- sort(d)
    d_sorted[2] - d_sorted[1] # gap between 1st and 2nd best
  })
  significant <- res_df[[fdr_col]] <= alpha
  weight_raw <- stats::plogis(margin * lambda) * significant

  if (rlang::is_installed("KernSmooth")) {
    # compute numeric bandwidth using KernSmooth::dpik; fall back to density() on error
    h <- tryCatch(
      KernSmooth::dpik(s, kernel = "normal"),
      error = function(e) NA_real_
    )
    if (is.na(h) || !is.finite(h) || h <= 0) {
      s_sub <- s[weight_raw > 0]
      if (length(s_sub) < 2) {
        s_sub <- s
      }
      dens0 <- stats::density(s_sub, n = 512)
      dens <- list(x = dens0$x, y = dens0$y)
    } else {
      dens <- KernSmooth::bkde(
        x = s,
        bandwidth = h,
        kernel = "normal",
        range.x = range(s, na.rm = TRUE) + c(-0.1, 0.1)
      )
    }
  } else {
    # original fallback: weighted selection + density()
    w <- (weight_raw - min(weight_raw)) /
      (max(weight_raw) - min(weight_raw) + .Machine$double.eps)
    w[w < 0.001] <- 0

    s_sub <- s[w > 0]
    dens0 <- stats::density(if (length(s_sub) >= 2) s_sub else s, n = 512)
    dens <- list(x = dens0$x, y = dens0$y)
  }

  x_grid <- dens$x
  y_grid <- dens$y

  # Find valleys (local minima)
  d1 <- diff(y_grid)
  d2 <- diff(sign(d1))
  valleys_idx <- which(d2 == 2) + 1

  # Default: use extremes of significant cells
  s_sig <- s[significant]
  if (length(s_sig) < 10) {
    s_sig <- s
  }

  if (length(valleys_idx) >= 2) {
    # 双峰：取最左/最右 valley
    left_v <- min(valleys_idx)
    right_v <- max(valleys_idx)
    thresh_low <- x_grid[left_v]
    thresh_high <- x_grid[right_v]
  } else if (length(valleys_idx) == 1) {
    center <- x_grid[valleys_idx]
    spread <- stats::mad(s_sig)
    thresh_low <- center - spread
    thresh_high <- center + spread
  } else {
    thresh_low <- stats::quantile(s_sig, 0.25)
    thresh_high <- stats::quantile(s_sig, 0.75)
  }

  class_vec <- dplyr::case_when(
    significant & s > thresh_high ~ "Positive",
    significant & s < thresh_low ~ "Negative",
    TRUE ~ "Neutral"
  )

  dplyr::mutate(
    res_df,
    PIPET_phenotype_score = s,
    PIPET_match_confidence = weight_raw,
    PIPET = class_vec
  )
}
