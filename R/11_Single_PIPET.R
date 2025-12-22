#' @title Single-Cell Subpopulation Prediction
#' @description
#' Predicts cell subpopulations in single-cell data by matching expression profiles
#' to predefined marker gene templates using various distance/similarity metrics.
#' This function implements a template-based classification approach with permutation
#' testing for significance assessment.
#'
#' @param sc_data A Seurat object containing single-cell RNA-seq data
#' @param markers A data frame containing marker gene information with columns:
#' - genes: Gene identifiers matching rownames of sc_data
#' - class: Class labels for each gene
#' @param freq_counts Integer threshold for minimum number of cells expressing a gene.
#' Genes expressed in fewer cells will be removed (default: NULL, no filtering)
#' @param normalize Logical indicating whether to normalize count data using CPM
#' and log1p transformation (default: TRUE)
#' @param scale Logical indicating whether to z-score normalize features (default: TRUE)
#' @param nPerm Number of permutations for null distribution estimation (default: 1000)
#' @param distance Distance/similarity metric to use. Options include:
#' - "cosine": Cosine similarity
#' - "pearson": Pearson correlation
#' - "spearman": Spearman correlation
#' - "kendall": Kendall correlation
#' - "euclidean": Euclidean distance
#' - "maximum": Maximum distance
#' @param ... Additional parameters including:
#' - seed: Random seed for reproducibility
#' - verbose: Whether to show progress messages
#' - parallel: Whether to use parallel processing
#'
#'
#' @return A data.frame with rows representing cells and columns containing:
#' \item{prediction}{Predicted subpopulation labels based on template matching}
#' \item{dist_*}{Distance/similarity scores for each class template}
#' \item{Pvalue}{Nominal p-value estimated from permutation testing}
#' \item{FDR}{False discovery rate adjusted p-values}
#'
#'
#' @examples
#' \dontrun{
#' # Run template matching
#' results <- PIPET_SingleAnalysis(
#' sc_data = seurat_obj,
#' markers = markers,
#' distance = "cosine",
#' nPerm = 1000
#' )
#'
#' # View results
#' head(results)
#' table(results$prediction)
#' }
#'
#' @references
#' Hoshida, Y. (2010). Nearest Template Prediction: A Single-Sample-Based Flexible
#' Class Prediction with Confidence Assessment. PLoS ONE 5, e15543.
#' @export
#'
#' @seealso
#' [PIPET_GroupAnalysis()] for group-wise analysis,
#' [corCosine()] for cosine similarity calculation,
#' [disFun()] for distance calculations
PIPET_SingleAnalysis <- function(
  sc_data,
  markers,
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
  distance <- SigBridgeRUtils::MatchArg(
    distance,
    c(
      "cosine",
      "pearson",
      "spearman",
      "kendall",
      "euclidean",
      "maximum"
    )
  )
  if (min(table(markers$class)) < 2) {
    cli::cli_warn(c(
      "x" = "Some classes have fewer than 2 marker genes, returning NULL"
    ))
    return(NULL)
  }

  dots <- rlang::list2(...)
  seed <- dots$seed %||% SigBridgeRUtils::getFuncOption("seed")
  verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption("verbose")
  parallel <- (dots$verbose %||% FALSE) &
    !inherits(future::plan("list")[[1]], "sequential")

  SC <- SeuratObject::LayerData(sc_data, assay = "RNA", layer = "counts")

  # 过滤markers以匹配单细胞数据
  keep_gene <- markers$genes %chin% rownames(SC)
  if (sum(keep_gene) == 0) {
    cli::cli_warn(c(
      "x" = "No overlapping genes between markers and single-cell data, returning NULL"
    ))
    return(NULL)
  }

  markers <- markers[keep_gene, ]

  if (verbose) {
    ts_cli$cli_alert_info(
      "The classification of markers is: {purrr::imap(table(markers$class), ~ paste0(.y, ': ', .x))}"
    )
    table(markers$class)
  }
  if (min(table(markers$class)) < 5) {
    cli::cli_warn(
      "The number of features in a class is less than 5, the predictions may be unstable."
    )
  }

  #  Match vector for SC and markers
  mm <- match(markers$genes, rownames(SC), nomatch = 0)
  if (!all(rownames(SC)[mm] == markers$genes)) {
    cli::cli_abort(c("x" = "Gene matching failed"))
  }

  set.seed(seed)

  # keep those genes expressed in more than 'freq_counts' cells
  if (!is.null(freq_counts)) {
    if (verbose) {
      ts_cli$cli_alert_info(
        "Filter cells with {.arg freq_counts} = {.val {freq_counts}}"
      )
    }
    nonzero <- SC > 0
    keep_genes <- Matrix::rowSums(nonzero) >= freq_counts
    SC <- SC[keep_genes, ]
  }

  # Normalizing and scaling SC data
  if (normalize) {
    if (verbose) {
      ts_cli$cli_alert_info(
        "Normalize count data with CPM and log1p"
      )
    }
    col_sums <- Matrix::colSums(SC)
    SC <- log1p(Matrix::t(
      Matrix::t(SC) / (col_sums * 1e4 + .Machine$double.eps)
    ))
  }
  if (scale) {
    if (verbose) {
      ts_cli$cli_alert_info(
        "Scale features with z-score normalization"
      )
    }
    # # z-score normalization
    # SC <- Matrix::t(scale(Matrix::t(SC), center=TRUE, scale=TRUE))

    n <- ncol(SC)
    mu <- Matrix::rowMeans(SC)
    sd_s <- sqrt(pmax(
      (Matrix::rowMeans(SC^2) - mu^2) * n / (n - 1),
      .Machine$double.eps
    ))
    SC <- (SC - mu) / sd_s
  }
  SC[is.na(SC)] <- 0

  # Prepare templates
  markers$class <- as.numeric(sub("group_", "", markers$class))
  class_names <- unique(markers$class)
  n_levels <- length(class_names)
  markers_num <- markers$class

  # markers matrix
  M_mat <- matrix(
    as.numeric(
      col(matrix(0, length(markers_num), n_levels)) == markers_num
    ),
    ncol = n_levels
  )
  if (n_levels == 2) {
    M_mat[M_mat == 0] <- -1
  }

  # 定义预测函数
  pred_fun <- function(
    n,
    distance,
    SC,
    mm = mm,
    nPerm = nPerm,
    M_mat = M_mat,
    n_levels = n_levels
  ) {
    # 距离和相关性转换函数
    if (distance %in% c("cosine", "pearson", "spearman", "kendall")) {
      # 计算相关性
      cor <- as.vector(corFun(
        x = SC[mm, n, drop = FALSE],
        y = M_mat,
        distance = distance
      ))

      # 置换检验
      perm_mat <- matrix(
        SC[, n][sample.int(
          nrow(SC),
          length(mm) * nPerm,
          replace = TRUE
        )],
        ncol = nPerm
      )
      cor_perm_max <- SigBridgeRUtils::rowMaxs3(
        x = corFun(x = perm_mat, y = M_mat, distance = distance)
      )

      pred <- which.max(cor)
      cor_ranks <- rank(-c(cor[pred], cor_perm_max))
      pval <- cor_ranks[1] / length(cor_ranks)
      dist <- CorToDist(cor)
    }

    if (distance %in% c("euclidean", "maximum")) {
      # A vector
      dist <- disFun(
        x = SC[mm, n, drop = FALSE],
        y = M_mat,
        distance = distance,
        n_levels = n_levels
      )
      cor <- DistToCor(dist)

      # 置换检验
      perm_mat <- matrix(
        SC[, n][sample.int(
          nrow(SC),
          length(mm) * nPerm,
          replace = TRUE
        )],
        ncol = nPerm
      )
      cor_perm_max <- SigBridgeRUtils::rowMaxs3(DistToCor(vapply(
        X = seq_len(ncol(perm_mat)),
        FUN = function(i) {
          disFun(
            x = perm_mat[, i],
            y = M_mat,
            distance = distance,
            n_levels = 2L
          )
        },
        FUN.VALUE = numeric(nrow(M_mat))
      )))

      pred <- which.min(dist)
      cor_ranks <- rank(-c(cor[pred], cor_perm_max))
      pval <- cor_ranks[1] / length(cor_ranks)
    }

    return(c(
      pred, # prediction
      dist, # distance
      pval # p-value
    ))
  }

  res <- if (parallel) {
    rlang::check_installed("furrr")
    if (verbose) {
      ts_cli$cli_alert_info(
        "Running parallel Prediction"
      )
    }
    furrr::future_map(
      .x = seq_len(ncol(SC)),
      .f = ~ pred_fun(
        n = .x,
        distance = distance,
        SC = SC,
        mm = mm,
        nPerm = nPerm,
        M_mat = M_mat,
        n_levels = n_levels
      ),
      .options = furrr::furrr_options(
        seed = seed,
        packages = c('SigBridgeRUtils', 'PIPET')
      ),
      .progress = verbose
    )
  } else {
    purrr::map(
      .x = seq_len(ncol(SC)),
      .f = ~ pred_fun(
        n = .x,
        distance = distance,
        SC = SC,
        mm = mm,
        nPerm = nPerm,
        M_mat = M_mat,
        n_levels = n_levels
      ),
      .progress = 'Prediction'
    )
  }

  if (verbose) {
    ts_cli$cli_alert_info(
      "Organize the computed results"
    )
  }

  # 格式化结果
  res_df <- data.frame(do.call(rbind, res))
  colnames(res_df) <- c(
    "PIPET_prediction",
    paste0("PIPET_dist_", class_names),
    "PIPET_Pvalue"
  )
  res_df$PIPET_prediction <- factor(
    class_names[res_df$PIPET_prediction],
    levels = class_names
  )
  rownames(res_df) <- colnames(SC)
  res_df$PIPET_FDR <- stats::p.adjust(res_df$PIPET_Pvalue, "fdr")

  # Add label
  if (n_levels == 2) {
    return(BinaryLabelCell(res_df))
  }

  ContinuousLabelCell(res_df)
}
