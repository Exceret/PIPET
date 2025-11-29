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
#' @param rm_NA Logical indicating whether to remove NA values (default: TRUE)
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
#' - parallel.type: Type of parallel backend
#' - workers: Number of parallel workers
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
#'
#' @seealso
#' [PIPET_GroupAnalysis()] for group-wise analysis,
#' [corCosine()] for cosine similarity calculation,
#' [disFun()] for distance calculations
PIPET_SingleAnalysis <- function(
    sc_data,
    markers,
    rm_NA = TRUE,
    freq_counts = NULL,
    normalize = TRUE,
    scale = TRUE,
    nPerm = 1000L,
    distance = "cosine",
    ...
) {
    if (min(table(markers$class)) < 2) {
        cli::cli_warn(c(
            "x" = "Some classes have fewer than 2 marker genes, returning NULL"
        ))
        return(NULL)
    }
    SC <- SeuratObject::LayerData(sc_data, assay = "RNA", layer = "counts")

    # 过滤markers以匹配单细胞数据
    keep_gene <- markers$genes %chin% rownames(SC)
    if (length(keep_gene) == 0) {
        cli::cli_warn(c(
            "x" = "No overlapping genes between markers and single-cell data, returning NULL"
        ))
        return(NULL)
    }
    markers <- markers[keep_gene, ]

    #  Match vector for SC and markers
    mm <- match(markers$genes, rownames(SC), nomatch = 0)
    if (!all(rownames(SC)[mm] == markers$genes)) {
        cli::cli_abort(c("x" = "Gene matching failed"))
    }

    dots <- rlang::list2(...)
    seed <- dots$seed %||% SigBridgeRUtils::getFuncOption("seed")
    verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption("verbose")
    parallel <- dots$parallel %||% SigBridgeRUtils::getFuncOption("parallel")
    parallel_type <- dots$parallel.type %||%
        SigBridgeRUtils::getFuncOption("parallel.type")
    workers <- dots$workers %||% SigBridgeRUtils::getFuncOption("workers")

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
        # z-score normalization
        gene_means <- Matrix::rowMeans(SC)
        gene_sds <- sqrt(Matrix::rowMeans(SC^2) - gene_means^2)
        SC <- (SC - gene_means) / (gene_sds + .Machine$double.eps)
    }

    # Prepare templates
    class_names <- unique(markers$class)
    n_levels <- length(class_names)
    markers_num <- as.numeric(markers$class)

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
    pred_fun <- function(n, distance) {
        call <- rlang::caller_env()
        distance <- SigBridgeRUtils::MatchArg(
            distance,
            c(
                "cosine",
                "pearson",
                "spearman",
                "kendall",
                "euclidean",
                "maximum"
            ),
            NULL,
            call = call
        )

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
            cor_perm_max <- SigBridgeRUtils::rowMaxs(
                x = corFun(x = perm_mat, y = M_mat)
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
            cor_perm_max <- SigBridgeRUtils::rowMaxs(DistToCor(vapply(
                seq_len(ncol(perm_mat)),
                function(i) {
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

    if (verbose) {
        ts_cli$cli_alert_info(
            "Running parallel computation with {workers} workers"
        )
    }

    res <- if (parallel) {
        SigBridgeRUtils::plan(parallel_type, workers = workers)
        on.exit(SigBridgeRUtils::plan('sequential'), add = TRUE)

        SigBridgeRUtils::future_map(
            .x = seq_len(ncol(SC)),
            .f = ~ pred_fun(n = .x, distance = distance),
            .options = furrr::furrr_options(
                seed = seed,
                packages = c('SigBridgeRUtils', 'PIPET'),
                globals = list(
                    SC = SC,
                    nPerm = nPerm,
                    mm = mm,
                    M_mat = M_mat
                )
            ),
            .progress = verbose
        )
    } else {
        purrr::map(
            .x = seq_len(ncol(SC)),
            .f = ~ pred_fun(n = .x, distance = distance),
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
    colnames(res_df) <- c("prediction", paste0("dist_", class_names), "Pvalue")
    res_df$prediction <- factor(
        class_names[res_df$prediction],
        levels = class_names
    )
    rownames(res_df) <- colnames(SC)
    res_df$FDR <- stats::p.adjust(res_df$Pvalue, "fdr")

    res_df
}
