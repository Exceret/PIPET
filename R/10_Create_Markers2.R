#' @title Create Markers on Log2-transformed Bulk Data
#' @description Based on differential expression analysis using limma (suitable for log2-transformed normalized expression data, e.g., log2(TPM+1), log2(FPKM+1), or log2(counts+1)),
#'     the feature vector of each subclass is established. The input bulk_data must be **log2-transformed continuous values**, not raw counts.
#'     At least two subclasses are required. Subclasses will be converted to factors if not already.
#'
#' @param bulk_data A numeric matrix of log2-transformed expression values (e.g., log2(x + 1)), with rows = genes and columns = samples.
#'     Must contain **non-integer (decimal) values**, if all integers, function aborts and suggests using \code{Create_Markers()}.
#' @param colData A \code{data.frame} with rows corresponding to columns of \code{bulk_data}.
#' @param class_col The name of the column in \code{colData} indicating subclass labels.
#' @param lg2FC The absolute log2 fold change cutoff (default: 1). Genes with |log2FC| >= lg2FC are retained.
#' @param p.adjust Adjusted p-value cutoff (default: 0.05).
#' @param show_lg2FC Whether to include log2 fold change in output (default: TRUE).
#' @param ... Additional arguments: \code{verbose} (default: \code{TRUE}), \code{seed} (ignored but accepted), \code{parallel} (passed to \code{limma::eBayes} if \code{trend = TRUE} uses \code{parallel = TRUE} via \code{BiocParallel}).
#'
#' @return A \code{data.frame} with columns:
#'   \item{genes}{Differential genes passing cutoffs.}
#'   \item{class}{The subclass to which the marker gene is upregulated.}
#'   \item{log2FoldChange}{(if \code{show_lg2FC = TRUE}) Signed log2 fold change.}
#'
#' @export
#'
Create_Markers2 <- function(
  bulk_data,
  colData,
  class_col = NULL,
  lg2FC = 1L,
  p.adjust = 0.05,
  show_lg2FC = TRUE,
  ...
) {
  #   all_integer <- sum(bulk_data == floor(bulk_data)) >
  #     ncol(bulk_data) * nrow(bulk_data) / 10

  #   if (all_integer) {
  #     cli::cli_abort(c(
  #       "x" = "Detected integer, non-negative {.arg bulk_data} probably like raw counts.",
  #       "!" = "limma (voom/lmFit) expects log2-transformed, continuous expression data (e.g., log2(x + 1)).",
  #       ">" = "For raw count data, please use {.fn Create_Markers} (based on DESeq2)."
  #     ))
  #   }

  colData_dt <- data.table::as.data.table(colData)

  if (!class_col %chin% colnames(colData)) {
    cli::cli_abort(c(
      "x" = "{.val {class_col}} must be a column in {.arg colData}",
      ">" = "Available columns are: {colnames(colData)}"
    ))
  }

  if (lg2FC < 0) {
    cli::cli_abort(c(
      "x" = "{.val {lg2FC}} must be >= 0",
      ">" = "Current: {.val {lg2FC}}"
    ))
  }
  if (p.adjust < 0 || p.adjust > 1) {
    cli::cli_abort(c(
      "x" = "{.val {p.adjust}} must be in [0, 1]",
      ">" = "Current: {.val {p.adjust}}"
    ))
  }

  # --- Parse ... ---
  dots <- rlang::list2(...)
  verbose <- dots$verbose %||%
    SigBridgeRUtils::getFuncOption("verbose") %||%
    TRUE
  seed <- dots$seed %||% SigBridgeRUtils::getFuncOption("seed") %||% 123L

  set.seed(seed)

  # Factor conversion
  if (is.numeric(colData_dt[[class_col]])) {
    colData_dt[[class_col]] <- paste0("group_", colData_dt[[class_col]])
  }
  if (!is.factor(colData_dt[[class_col]])) {
    data.table::set(
      colData_dt,
      j = class_col,
      value = as.factor(colData_dt[[class_col]])
    )
  }
  colData_dt$class <- colData_dt[[class_col]]
  class_levels <- levels(colData_dt$class)

  if (length(class_levels) < 2) {
    cli::cli_abort(c(
      "x" = "At least two classes required",
      ">" = "Current classes: {.val {class_levels}}"
    ))
  }

  # Design matrix setup
  design <- stats::model.matrix(~ 0 + class, data = colData_dt)
  colnames(design) <- class_levels
  # limma pipeline
  if (verbose) {
    ts_cli$cli_alert_info("Fitting model with limma")
  }

  fit <- limma::lmFit(bulk_data, design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

  if (length(class_levels) == 2) {
    # Two-class: contrast = class2 - class1
    contrast_mat <- limma::makeContrasts(
      contrasts = paste0(class_levels[2], "-", class_levels[1]),
      levels = design
    )

    fit2 <- limma::contrasts.fit(fit, contrast_mat)
    fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)

    res <- as.data.frame(limma::topTable(fit2, number = Inf, sort.by = "P"))
    res_dt <- data.table::as.data.table(res, keep.rownames = "genes")

    if ("ID" %in% colnames(res_dt)) {
      UP <- res_dt[
        adj.P.Val < p.adjust & logFC >= lg2FC,
        .(ID, class = class_levels[2], log2FoldChange = logFC)
      ]
      DOWN <- res_dt[
        adj.P.Val < p.adjust & logFC <= -lg2FC,
        .(ID, class = class_levels[1], log2FoldChange = logFC)
      ]
    } else {
      UP <- res_dt[
        adj.P.Val < p.adjust & logFC >= lg2FC,
        .(genes, class = class_levels[2], log2FoldChange = logFC)
      ]
      DOWN <- res_dt[
        adj.P.Val < p.adjust & logFC <= -lg2FC,
        .(genes, class = class_levels[1], log2FoldChange = logFC)
      ]
    }

    if (!show_lg2FC) {
      UP[, log2FoldChange := NULL]
      DOWN[, log2FoldChange := NULL]
    }
    markers <- data.table::rbindlist(list(DOWN, UP), use.names = TRUE)

    if (nrow(markers) == 0) {
      cli::cli_warn(
        "No markers found, please try different {.arg lg2FC} and {.arg p.adjust}"
      )
    } else if (verbose) {
      ts_cli$cli_alert_success("Created {.val {nrow(markers)}} marker genes")
    }

    return(as.data.frame(markers))
  }

  # Multi-class: one-vs-rest
  if (verbose) {
    cli::cli_progress_bar(
      name = "limma one-vs-rest",
      type = "tasks",
      total = length(class_levels)
    )
  }

  markers_list <- vector("list", length(class_levels))
  for (i in seq_along(class_levels)) {
    # Contrast: class_i vs mean of others
    others <- setdiff(class_levels, class_levels[i])
    contrast_vec <- rep(0, length(class_levels))
    names(contrast_vec) <- class_levels

    contrast_vec[class_levels[i]] <- 1
    contrast_vec[others] <- -1 / length(others)

    contrast_mat <- matrix(
      contrast_vec,
      nrow = length(class_levels),
      ncol = 1
    )
    rownames(contrast_mat) <- class_levels
    colnames(contrast_mat) <- paste0(class_levels[i], "_vs_rest")

    fit_contrast <- limma::contrasts.fit(fit, contrast_mat)
    fit_contrast <- limma::eBayes(fit_contrast, trend = TRUE, robust = TRUE)

    res <- as.data.frame(limma::topTable(
      fit_contrast,
      number = Inf,
      sort.by = "P"
    ))
    res_dt <- data.table::as.data.table(res, keep.rownames = "genes")

    UP <- res_dt[
      adj.P.Val < p.adjust & logFC >= lg2FC
    ]

    if (!"ID" %in% colnames(res_dt)) {
      if (show_lg2FC) {
        UP <- UP[, .(genes, class = class_levels[i], log2FoldChange = logFC)]
      } else {
        UP <- UP[, .(genes, class = class_levels[i])]
      }
    } else {
      if (show_lg2FC) {
        UP <- UP[, .(ID, class = class_levels[i], log2FoldChange = logFC)]
      } else {
        UP <- UP[, .(ID, class = class_levels[i])]
      }
    }

    markers_list[[i]] <- UP
    if (verbose) cli::cli_progress_update()
  }
  if (verbose) {
    cli::cli_progress_done()
  }

  markers <- data.frame(do.call(rbind, markers_list))

  if (nrow(markers) == 0) {
    cli::cli_warn(
      "No markers found, please try different {.arg lg2FC} and {.arg p.adjust}"
    )
  } else if (verbose) {
    ts_cli$cli_alert_success(
      "Created {.val {nrow(markers)}} marker genes for {.val {length(class_levels)}} classes"
    )
    markers$comparison <- paste0(seq_along(class_levels), "_vs_rest")
  }

  markers
}
