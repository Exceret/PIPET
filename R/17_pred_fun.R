#' @keywords internal
pred_fun <- function(
  n, # cell numbers
  distance,
  SC,
  mm = mm, # gene match
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
    cor_perm_max <- SigBridgeRUtils::rowMaxs3(
      x = DistToCor(
        x = vapply(
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
        )
      )
    )

    pred <- which.min(dist)
    cor_ranks <- rank(-c(cor[pred], cor_perm_max))
    pval <- cor_ranks[1] / length(cor_ranks)
  }

  return(c(
    prediction = pred, # prediction
    distance = dist, # distance
    pval = pval # p-value
  ))
}
