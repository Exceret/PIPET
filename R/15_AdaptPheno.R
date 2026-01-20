#' @title Adapt Phenotype Vector for Downstream Analysis
#' @description
#' Converts a named numeric phenotype vector into a standardized data frame with
#' discrete class labels. Supports both binary, continuous, and survival phenotypes.
#'
#' For continuous phenotypes, performs automatic discretization based on the
#' specified method (e.g., quantile-based, k-means clustering, or custom cutoffs).
#' Internally applies log2(x+1) transformation and z-scoring before discretization
#' to stabilize variance and improve cluster separation.
#'
#'
#' @param phenotype A named numeric vector (for \code{"binary"} or \code{"continuous"})
#'   or a data frame with columns \code{time} and \code{status} (for \code{"survival"}).
#'   - For \code{"binary"}: values should be 0/1 (or equivalent); exactly two unique values required.
#'   - For \code{"continuous"}: numeric values; must have >2 unique values unless explicitly declared binary.
#'   - For \code{"survival"}: a data frame where \code{rownames} are sample IDs, and \code{status} is 0/1.
#' @param discretize_method \code{c("median", "kmeans", "custom")}. Discretization
#'   strategy for continuous phenotypes. Note: `"median"` is mapped internally to
#'   `"quantile"` (2-group quantile split). Default: `"kmeans"`.
#' @param cutoff Numeric vector of length `n_group - 1`. Required only when
#'   \code{discretize_method = "custom"}. Defines interior breakpoints on the
#'   *normalized, log2-transformed scale* (i.e., after `scale(log2(x + 1))`).
#'   Must be sorted in ascending order.
#' @param phenotype_type \code{c("binary", "continuous", "survival")}. If missing or
#'   length > 1, auto-detection is attempted:
#' @param ... Additional arguments (currently unused, reserved for future extension).
#'
#' @return A two-column \code{data.frame}:
#'   \describe{
#'     \item{sample}{Sample IDs (from \code{names(phenotype)})}
#'     \item{class}{Discrete class labels as a factor. For binary phenotypes:
#'       `"group_1"` (value == 1) and `"group_0"` (value != 1). For continuous:
#'       `"group_1"`, `"group_2"`, ..., ordered by increasing mean of normalized
#'       phenotype within group (for \code{"kmeans"}), or by quantile/custom interval.}
#'   }
#'
#' @section Notes:
#' \itemize{
#'   \item Binary phenotype validation: exactly two unique values required; if more,
#'     use \code{phenotype_type = "continuous"}.
#'   \item Continuous phenotype with only two unique values will error — explicitly
#'     declare as binary in such edge cases (e.g., dichotomized lab values).
#'   \item Discretization is performed on \code{scale(log2(phenotype + 1))} to
#'     ensure robustness for skewed genomic/clinical data (e.g., gene expression,
#'     biomarker levels). Cutoffs must be specified on this transformed scale.
#' }
#'
#' @examples
#' \dontrun{
#' # Binary phenotype (auto-detected)
#' pheno_bin <- c(A = 1, B = 0, C = 1, D = 0)
#' AdaptPheno(pheno_bin)
#'
#' # Continuous phenotype: 3-group quantile split
#' set.seed(123)
#' pheno_cont <- rlnorm(20, meanlog = 5, sdlog = 1)
#' names(pheno_cont) <- paste0("S", 1:20)
#' AdaptPheno(pheno_cont, phenotype_type = "continuous", n_group = 3)
#'
#' # Custom cutoffs (e.g., normalized expression > 1.5 = high)
#' # Suppose we want 2 groups with cutoff at z = 0.5 on normalized scale
#' AdaptPheno(pheno_cont,
#'            phenotype_type = "continuous",
#'            discretize_method = "custom",
#'            cutoff = 0.5,
#'            n_group = 2)
#' }
#' @export
AdaptPheno <- function(
  phenotype,
  discretize_method = c("kmeans", "median", "custom"),
  cutoff = NULL,
  phenotype_type = c("binary", "continuous", "survival"),
  ...
) {
  # * Auto-detect phenotype type
  if (is.vector(phenotype)) {
    phenotype_type <- ifelse(
      length(table(phenotype)) == 2,
      "binary",
      "continuous"
    )
  } else if (is.data.frame(phenotype)) {
    phenotype_type <- "survival"
  } else {
    cli::cli_abort("Unable to auto-detect phenotype type, please specify")
  }

  switch(
    phenotype_type,
    "continuous" = {
      n_group <- length(table(phenotype))

      if (n_group == 2) {
        cli::cli_abort(c(
          "x" = "Only 2 classes in continuous phenotype, use {.arg phenotype_type} = {.val binary} instead"
        ))
      }

      # Discretization
      group <- DiscretizePheno(
        phenotype = phenotype,
        method = discretize_method,
        cutoff = cutoff,
        ngroup = n_group
      )

      data.frame(sample = names(phenotype), class = factor(group))
    },
    "binary" = {
      n_group <- length(table(phenotype))

      if (n_group > 2) {
        cli::cli_abort(c(
          "x" = "Only 2 classes allowed in binary phenotype, use {.arg phenotype_type} = {.val continuous} instead when {.arg n_group} > 2"
        ))
      }
      data.frame(
        sample = names(phenotype),
        class = factor(ifelse(phenotype == 1, "group_1", "group_0"))
      )
    },
    "survival" = {
      data.frame(
        sample = rownames(phenotype),
        class = paste0("group_", phenotype$status)
      )
    }
  )
}

#' @keywords internal
DiscretizePheno <- function(
  phenotype,
  method = c("kmeans", "quantile", "custom"),
  cutoff = NULL,
  n_group = NULL,
  ...
) {
  # * default method is "kmeans"
  method <- SigBridgeRUtils::MatchArg(
    method,
    c("kmeans", "quantile", "custom")
  )

  if (method == "custom") {
    if (is.null(cutoff)) {
      cli::cli_abort(c(
        "x" = "{.arg cutoff} must be provided when {.arg method} is {.val custom}"
      ))
    }
    if (!is.numeric(cutoff) || length(cutoff) != n_group - 1) {
      cli::cli_abort(c(
        "x" = "{.arg cutoff} must be a vector of {.val {n_group-1}} numeric value"
      ))
    }
  }

  pheno_norm <- as.vector(scale(log2(phenotype + 1)))

  switch(
    method,
    "quantile" = {
      probs <- seq(0, 1, length.out = n_group + 1)
      breaks <- stats::quantile(
        pheno_norm,
        probs = probs,
        na.rm = TRUE,
        names = FALSE
      )

      cut(
        pheno_norm,
        breaks = breaks,
        labels = paste0("group_", seq_len(n_group)),
        include.lowest = TRUE,
        right = FALSE
      )
    },
    "kmeans" = {
      km <- stats::kmeans(matrix(pheno_norm, ncol = 1L), centers = n_group)
      group_temp <- paste0("group_", km$cluster)
      group_order <- names(sort(tapply(pheno_norm, group_temp, mean)))
      factor(group_temp, levels = group_order)
    },
    "custom" = {
      if (!is.numeric(cutoff) || length(cutoff) != n_group - 1) {
        cli::cli_abort(c(
          "x" = "{.arg cutoff} must be a vector of {.val {n_group-1}} numeric value"
        ))
      }
      breaks <- c(min(pheno_norm), cutoff, max(pheno_norm))
      cut(
        pheno_norm,
        breaks = breaks,
        labels = paste0("group_", seq_len(n_group)),
        include.lowest = TRUE,
        right = FALSE
      )
    },
    cli::cli_abort(c(
      "x" = "Unsupported {.arg method}: {.val {method}}",
      ">" = "Supported methods: median, kmeans, custom"
    ))
  )
}
