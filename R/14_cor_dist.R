#' @title Distance and Correlation Conversion Functions
#' @name dist-cor-conversion
#' @description
#' Utility functions for converting between distance metrics and correlation measures.
#' These functions are useful for data preprocessing and similarity calculations.
NULL

#' @rdname dist-cor-conversion
#' @description
#' DistToCor converts squared Euclidean distances to correlation coefficients.
#' Assumes data has been normalized to unit variance.
#'
#' @param x Numeric vector or matrix of squared Euclidean distances
#'
#' @return Correlation coefficients ranging from -1 to 1
#'
#' @details
#' The conversion formula is: correlation = 1 - 2 * distance^2
#' This assumes that the input distances are squared Euclidean distances
#' computed from standardized data (mean 0, variance 1).
#'
#' @examples
#' # Convert squared distances to correlations
#' distances <- c(0, 0.25, 0.5, 1)
#' correlations <- DistToCor(distances)
#' print(correlations)
#'
#' @export
DistToCor <- function(x) 1 - 2 * x^2

#' @rdname dist-cor-conversion
#' @description
#' CorToDist converts correlation coefficients to Euclidean distances.
#' Useful for converting correlation matrices to distance matrices for clustering.
#'
#' @param x Numeric vector or matrix of correlation coefficients (-1 to 1)
#'
#' @return Euclidean distances (non-negative values)
#'
#' @details
#' The conversion formula is: distance = sqrt(0.5 * (1 - correlation))
#' This produces Euclidean distances that are compatible with most
#' distance-based algorithms.
#'
#' @examples
#' # Convert correlations to distances
#' correlations <- c(1, 0.5, 0, -0.5, -1)
#' distances <- CorToDist(correlations)
#' print(distances)
#'
#' @export
CorToDist <- function(x) sqrt(1 / 2 * (1 - x))

#' @title Cosine Similarity
#' @description
#' corCosine computes cosine similarity between vectors or matrices.
#' Cosine similarity measures the cosine of the angle between vectors,
#' providing a scale-invariant measure of similarity.
#'
#' @param x Numeric matrix where columns represent vectors
#' @param y Optional second numeric matrix for cross-correlation. If NULL,
#' computes self-correlation of x
#'
#' @return Cosine similarity matrix. For self-correlation, returns a symmetric
#' matrix; for cross-correlation, returns a matrix with dimensions
#' ncol(x) × ncol(y)
#'
#' @details
#' Cosine similarity is calculated as:
#' cos(θ) = (x·y) / (||x|| * ||y||)
#'
#' This measure is particularly useful for:
#' - Text analysis (TF-IDF vectors)
#' - Recommendation systems
#' - High-dimensional data where magnitude is less important than direction
#'
#' The function handles zero vectors gracefully by setting their norm to 1,
#' preventing division by zero.
#'
#' @examples
#' # Self-correlation (symmetric matrix)
#' x <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
#' cos_sim_self <- corCosine(x)
#' print(cos_sim_self)
#'
#' # Cross-correlation
#' y <- matrix(c(2, 4, 6, 1, 3, 5), nrow = 3, ncol = 2)
#' cos_sim_cross <- corCosine(x, y)
#' print(cos_sim_cross)
#'
#' # With zero vectors
#' z <- matrix(c(1, 0, 0, 0, 0, 0), nrow = 3, ncol = 2)
#' cos_sim_zero <- corCosine(z)
#' print(cos_sim_zero)
#'
#' @export
corCosine <- function(x, y = NULL) {
    x <- as.matrix(x)

    if (is.null(y)) {
        # 自相关情况,只需计算一次
        x_norm <- sqrt(colSums(x^2))
        result <- crossprod(x) / outer(x_norm, x_norm)
        return(result)
    }

    y <- as.matrix(y)
    x_norm <- sqrt(colSums(x^2))
    y_norm <- sqrt(colSums(y^2))

    crossprod(x, y) / outer(x_norm, y_norm)
}

#' @title Distance Calculation Between Matrices
#' @description
#' Computes pairwise distances between columns of two matrices using specified
#' distance metrics. Useful for comparing feature distributions or sample similarities.
#'
#' @param x First numeric matrix (features × samples)
#' @param y Second numeric matrix (features × samples)
#' @param distance Distance metric to use. Common options include:
#' - "euclidean": Euclidean distance
#' - "manhattan": Manhattan distance
#' - "maximum": Maximum distance
#' - "canberra": Canberra distance
#' - "binary": Binary distance
#' - "minkowski": Minkowski distance
#' @param n_levels Number of distance values to return. Typically the number of
#' columns in the smaller matrix or a specified subset.
#'
#' @return Numeric vector of distances between corresponding columns of x and y
#'
#' @details
#' This function combines the columns of both matrices and computes a distance matrix
#' using R's built-in dist function, then extracts the relevant pairwise distances
#' between columns of x and y.
#'
#' The function is particularly useful for:
#' - Comparing gene expression profiles between different conditions
#' - Measuring similarity between single-cell and bulk expression patterns
#' - Feature selection based on distribution differences
#'
#' @examples
#' # Create example matrices
#' x <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
#' y <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
#'
#' # Compute Euclidean distances between corresponding columns
#' distances <- disFun(x, y, distance = "euclidean", n_levels = 5)
#' print(distances)
#'
#' # Compute Manhattan distances
#' manhattan_dists <- disFun(x, y, distance = "manhattan", n_levels = 5)
#' print(manhattan_dists)
#'
#' @seealso
#' [stats::dist()] for the underlying distance calculation
#' @export
disFun <- function(x, y, distance, n_levels = 2L) {
    tmp <- stats::dist(t(cbind(x, y)), method = distance)
    as.vector(tmp)[seq_len(n_levels)]
}

#' @title Correlation Calculation Function
#' @description
#' Creates a correlation function for computing similarities between matrices using
#' either cosine similarity or standard correlation methods.
#'
#' @param x First numeric matrix (features × samples)
#' @param y Second numeric matrix (features × samples)
#' @param distance Correlation method to use. Options include:
#' - "cosine": Cosine similarity (scale-invariant)
#' - "pearson": Pearson correlation (linear relationship)
#' - "spearman": Spearman correlation (rank-based)
#' - "kendall": Kendall correlation (rank-based, more robust)
#'
#' @return A correlation function that can be applied to compute similarities
#' between the input matrices
#'
#' @details
#' This function returns an appropriate correlation function based on the specified
#' distance metric. For cosine similarity, it uses the optimized corCosine function;
#' for other methods, it returns the corresponding correlation function from R's
#' stats::cor function.
#'
#' Cosine similarity is particularly useful when the magnitude of vectors is less
#' important than their direction, such as in text analysis or recommendation systems.
#' Pearson correlation measures linear relationships, while Spearman and Kendall
#' correlations are rank-based and more robust to outliers.
#'
#' @examples
#' # Create example matrices
#' x <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
#' y <- matrix(rnorm(100 * 3), nrow = 100, ncol = 3)
#'
#' # Get cosine similarity 
#' cosine_similarity <- corFun(x, y, distance = "cosine")
#'
#' # Get Pearson correlation 
#' pearson_cor <- corFun(x, y, distance = "pearson")
#'
#' @seealso
#' [corCosine()] for cosine similarity, [stats::cor()] for standard correlation methods
#' @export
corFun <- function(x, y, distance) {
    if (distance == "cosine") {
        corCosine(x, y)
    } else {
        function(x, y) stats::cor(x, y, method = distance)
    }
}
