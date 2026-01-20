#' @title Run PIPET Screening
#' @description
#' Predicts cell subpopulations in single-cell data by matching expression profiles
#' to predefined marker gene templates using various distance/similarity metrics.
#'
#'
#' @param sc_data A Seurat object of single cell data.
#' @param markers A data frame of phenotypic information from bulk data.
#' @param group A character, name of one metadata column to group cells by (for example, orig.ident).
#' @param freq_counts An integer, keep genes expressed in more than a certain number of cells.
#' @param normalize Select whether to perform normalization of count data. The default value is TRUE.
#' @param scale Select whether to scale and center features in the dataset. The default value is TRUE.
#' @param nPerm An integer, number of permutations to do. The default value is 1000.
#' @param distance A character, the distance algorithm must be included in "cosine", "pearson", "spearman", "kendall","euclidean","maximum".
#' @param ... Additional arguments to be passed to \code{\link{PIPET_SingleAnalysis}}.
#' - seed: Random seed for reproducibility
#' - verbose: Whether to show progress messages
#' - parallel: Whether to use parallel processing
#'
#'
#' @export
#' @seealso [PIPET_SingleAnalysis()], [PIPET_GroupAnalysis()]
#'
PIPET <- function(
  sc_data,
  markers,
  group = NULL,
  freq_counts = NULL,
  normalize = TRUE,
  scale = TRUE,
  nPerm = 1000L,
  distance = c(
    "cosine",
    "pearson",
    "spearman",
    "kendall",
    "euclidean",
    "maximum"
  ),
  ...
) {
  # Handle group-wise analysis
  if (!is.null(group)) {
    return(PIPET_GroupAnalysis(
      Seurat_data = sc_data,
      markers = markers,
      group = group,
      freq_counts = freq_counts,
      normalize = normalize,
      scale = scale,
      nPerm = nPerm,
      distance = distance,
      ...
    ))
  }
  PIPET_SingleAnalysis(
    sc_data = sc_data,
    markers = markers,
    freq_counts = freq_counts,
    normalize = normalize,
    scale = scale,
    nPerm = nPerm,
    distance = distance,
    ...
  )
}
