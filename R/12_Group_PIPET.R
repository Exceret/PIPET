#' @title Grouping to predict subpopulations in single-cell data
#' @description Make group predictions for relevant subpopulations in single-cell data from phenotypic information in bulk data.
#'
#' @param Seurat_data A Seurat object of single cell data.
#' @param markers A data frame of phenotypic information from bulk data.
#' @param group A character, name of one metadata column to group cells by (for example, orig.ident).
#' @param rm_NA Select Whether to remove NA values. The default value is TRUE.
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
#' @return  This function returns a \code{data.frame} with rows are cells and the columns contain the following attributes:
#'     \item{prediction}{Subpopulation labels corresponding to single cell data determined based on distance.}
#'     \item{dist_}{Distance or similarity for each subclasses (determined by the chosen distance algorithm).}
#'     \item{Pvalue}{A nominal p-value estimated based on a null distribution for the distance  to the feature vectors.}
#'     \item{FDR}{Adjusted p values with false discovery rate.}
#'
#' @export
#'
#' @references Hoshida, Y. (2010). Nearest Template Prediction: A Single-Sample-Based Flexible Class Prediction with Confidence Assessment. PLoS ONE 5, e15543.
#'
PIPET_GroupAnalysis <- function(
    Seurat_data,
    markers,
    group,
    rm_NA = TRUE,
    freq_counts = NULL,
    normalize = TRUE,
    scale = TRUE,
    nPerm = 1000,
    distance = "cosine",
    ...
) {
    dots <- rlang::list2(...)
    verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption("verbose")
    dots$verbose <- FALSE # suppress verbose in single process

    g <- unique(Seurat_data[[]][[group]])
    if (is.null(g)) {
        cli::cli_abort(c(
            "x" = "Group column {.arg {group}} not found in Seurat metadata",
            ">" = "Available groups are: {colnames(Seurat_data[[]])}"
        ))
    }
    if (verbose) {
        ts_cli$cli_alert_info(
            "Performing group-wise analysis on {length(g)} groups"
        )
    }

    result <- data.frame()
    if (verbose) {
        cli::cli_progress_bar(
            name = "Processing groups",
            type = 'tasks',
            total = length(g)
        )
    }
    for (group in g) {
        sc_sub <- subset(Seurat_data, idents = group)
        res <- do.call(
            PIPET_SingleAnalysis,
            c(
                list(
                    sc_data = sc_sub,
                    markers = markers,
                    rm_NA = rm_NA,
                    freq_counts = freq_counts,
                    normalize = normalize,
                    scale = scale,
                    nPerm = nPerm,
                    distance = distance
                ),
                dots
            )
        )
        result <- rbind(result, res)
        if (verbose) {
            cli::cli_progress_update()
        }
        gc(verbose = FALSE)
    }
    if (verbose) {
        cli::cli_progress_done()
    }

    rownames(result) <- make.names(g)

    result
}
