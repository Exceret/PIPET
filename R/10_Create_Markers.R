#' @title Create Markers on Raw Count Bulk Data
#' @description Based on the differential genes of bulk data (differential expression analysis was performed by DESeq2),
#'     the feature vector of each subclass is established. The bulk data to be analyzed should have at least two subclasses.
#'     You could customize the levels of subclasses, otherwise the subclasses will be automatically converted into a factor variable.
#'
#' @param bulk_data A matrix of non-negative integers with row features and sample columns.
#' @param colData A data.frame with at least a single column. Rows of colData correspond to columns of bulk_data.
#' @param class_col The variable name of subclasses, and must be contained in the colData columns.
#' @param lg2FC In the DESeq differential expression analysis results, the cutoff value of lg2FC. The default value is 1.
#' @param p.adjust In the DESeq differential expression analysis results, the cutoff value of adjust P. The default value is 0.05.
#' @param show_lg2FC Select whether to show log2 fold changes. The default value is TRUE.
#' @param ... Additional parameters like verbose (default `TRUE`), seed (default `123L`), and parallel (default `FALSE`).
#'
#' @return This function returns a \code{data.frame} with rows are samples and the columns contain the following attributes:
#'     \item{genes}{Differential genes screened out based on lg2FC and p.adjust.}
#'     \item{class}{The subclass of bulk data to which the marker gene belongs.}
#'
#' @export
#'
#'
Create_Markers <- function(
  bulk_data,
  colData,
  class_col = NULL,
  lg2FC = 1L,
  p.adjust = 0.05,
  show_lg2FC = TRUE,
  ...
) {
  if (any(bulk_data != floor(bulk_data), na.rm = TRUE)) {
    cli::cli_abort(c(
      "x" = "{.arg bulk_data} must contain only non-negative integers (raw counts)",
      ">" = "Detected non-integer values (e.g., decimals). DESeq2 requires raw integer counts",
      ">" = "Please use {.fun Create_Markers2} for log-transformed data"
    ))
  }
  colData_dt <- data.table::as.data.table(colData)

  if (!class_col %chin% colnames(colData)) {
    cli::cli_abort(c(
      "x" = "{.val {class_col}} must be a column in {.arg colData}",
      ">" = "Available columns are: {colnames(colData)}"
    ))
  }
  if (!is.numeric(colData_dt[[class_col]])) {
    cli::cli_abort(c(
      "x" = "{.val {class_col}} must be a numeric column"
    ))
  }
  if (lg2FC < 0) {
    cli::cli_abort(c(
      "x" = "{.val {lg2FC}} must be greater than or equal to 0",
      ">" = "Current value: {.val {lg2FC}}"
    ))
  }
  if (p.adjust < 0 || p.adjust > 1) {
    cli::cli_abort(c(
      "x" = "{.val {p.adjust}} must be between 0 and 1",
      ">" = "Current value: {.val {p.adjust}}"
    ))
  }

  # * default value
  dots <- rlang::list2(...)
  verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption('verbose')
  seed <- dots$seed %||% SigBridgeRUtils::getFuncOption('seed')
  parallel <- dots$seed %||% FALSE

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
  c <- levels(colData_dt$class)

  if (length(c) < 2) {
    cli::cli_abort(c(
      "x" = "At least two classes required",
      ">" = "Current class in column {.val {class_col}}: {.val {c}}"
    ))
  }

  if (length(c) == 2) {
    if (verbose) {
      ts_cli$cli_alert_info(
        "Two classes detected, using two-class comparison"
      )
    }
    # 二分类比较
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = bulk_data,
      colData = as.data.frame(colData_dt),
      design = ~class
    )
    dds <- DESeq2::DESeq(dds, quiet = !verbose)
    res <- as.data.frame(DESeq2::results(dds, parallel = parallel))

    res_dt <- data.table::as.data.table(res, keep.rownames = "genes")

    UP <- res_dt[
      padj < p.adjust & log2FoldChange >= lg2FC,
      .(genes, class = c[2], log2FoldChange)
    ]

    DOWN <- res_dt[
      padj < p.adjust & log2FoldChange <= -lg2FC,
      .(genes, class = c[1], log2FoldChange)
    ]

    if (!show_lg2FC) {
      UP[, log2FoldChange := NULL]
      DOWN[, log2FoldChange := NULL]
    }

    markers <- data.table::rbindlist(list(DOWN, UP), use.names = TRUE)
  } else {
    # 多分类比较
    if (verbose) {
      ts_cli$cli_alert_info(
        "More than two classes detected, using multi-class comparison"
      )
    }
    markers_list <- vector("list", length(c))

    if (verbose) {
      cli::cli_progress_bar(
        name = "Creating markers",
        type = "tasks",
        total = length(c)
      )
    }
    for (i in seq_along(c)) {
      colData_dt[, class_new := ifelse(class == c[i], c[i], "Others")]
      colData_dt[,
        class_new := factor(
          class_new,
          levels = c("Others", c[i])
        )
      ]

      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = bulk_data,
        colData = as.data.frame(colData_dt),
        design = ~class_new
      )
      dds <- DESeq2::DESeq(dds, quiet = !verbose)
      res <- as.data.frame(DESeq2::results(dds, parallel = parallel))

      res_dt <- data.table::as.data.table(res, keep.rownames = "genes")

      UP <- res_dt[
        padj < p.adjust & log2FoldChange >= lg2FC
      ]

      if (show_lg2FC) {
        UP <- UP[, .(genes, class = c[i], log2FoldChange)]
      } else {
        UP <- UP[, .(genes, class = c[i])]
      }

      markers_list[[i]] <- UP
      if (verbose) {
        cli::cli_progress_update()
      }
    }
    if (verbose) {
      cli::cli_progress_done()
    }

    markers <- data.table::rbindlist(markers_list, use.names = TRUE)
  }

  if (verbose) {
    ts_cli$cli_alert_success("Created {.val {nrow(markers)}} marker genes")
  }

  as.data.frame(markers)
}
