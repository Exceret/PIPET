#' @title Single-Cell Subpopulation Prediction using Template Matching
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
        cli::cli_abort(c("x" = "Some classes have fewer than 2 marker genes"))
    }
    SC <- SeuratObject::LayerData(sc_data)

    # 过滤markers以匹配单细胞数据
    keep_gene <- markers$genes %chin% rownames(SC)
    if (length(keep_gene) == 0) {
        cli::cli_abort(c(
            "x" = "No overlapping genes between markers and single-cell data"
        ))
    }
    markers <- markers[keep_gene, ]

    # 匹配基因
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

    # 设置并行计划
    set.seed(seed)
    SigBridgeRUtils::plan(parallel_type, workers = workers)
    on.exit(SigBridgeRUtils::plan('sequential'), add = TRUE)

    if (!is.null(freq_counts)) {
        nonzero <- SC > 0
        keep_genes <- Matrix::rowSums(nonzero) >= freq_counts
        SC <- SC[keep_genes, ]
    }

    # 标准化和缩放 - 使用 Matrix 操作
    if (normalize) {
        # 计算CPM并使用log1p
        col_sums <- Matrix::colSums(SC)
        SC <- log1p(Matrix::t(
            Matrix::t(SC) / (col_sums * 1e4 + .Machine$double.eps)
        ))
    }
    if (scale) {
        # 对每个基因进行z-score标准化
        gene_means <- Matrix::rowMeans(SC)
        gene_sds <- sqrt(Matrix::rowMeans(SC^2) - gene_means^2)
        SC <- (SC - gene_means) / (gene_sds + .Machine$double.eps)
    }

    # 准备模板矩阵
    class_names <- unique(markers$class)
    n_levels <- length(class_names)
    markers_num <- as.numeric(markers$class) # ! confusing

    M_mat <- matrix(rep(markers_num, n_levels), ncol = n_levels)
    for (i in seq_len(n_levels)) {
        M_mat[, i] <- as.numeric(M_mat[, i] == i)
    }
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

        corFun <- if (distance == "cosine") {
            corCosine
        } else {
            function(x, y) stats::cor(x, y, method = distance)
        }

        # 距离和相关性转换函数
        if (distance %in% c("cosine", "pearson", "spearman", "kendall")) {
            # 计算相关性
            cor <- as.vector(corFun(SC[mm, n, drop = FALSE], M_mat))

            # 置换检验
            perm_mat <- matrix(
                SC[, n][sample.int(
                    nrow(SC),
                    length(mm) * nPerm,
                    replace = TRUE
                )],
                ncol = nPerm
            )
            cor_perm_max <- matrixStats::rowMaxs(corFun(perm_mat, M_mat))

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
            cor_perm_max <- matrixStats::rowMaxs(DistToCor(t(apply(
                perm_mat,
                2,
                disFun,
                y = M_mat
            ))))

            pred <- which.min(dist)
            cor_ranks <- rank(-c(cor[pred], cor_perm_max))
            pval <- cor_ranks[1] / length(cor_ranks)
        }

        return(c(pred, dist, pval))
    }
    # 使用 furrr 进行并行计算
    if (verbose) {
        ts_cli$cli_alert_info(
            "Running parallel computation with {workers} workers"
        )
    }

    # 使用 furrr 的 future_map
    res <- furrr::future_map(
        .x = seq_len(ncol(SC)),
        .f = ~ pred_fun(.x),
        .options = furrr::furrr_options(seed = seed),
        .progress = verbose
    )

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
