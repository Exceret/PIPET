# test_that("multiplication works", {
#   raw_corCosine <- function(x, y = NULL) {
#     x <- as.matrix(x)

#     if (is.null(y)) {
#       # 自相关情况,只需计算一次
#       x_norm <- sqrt(SigBridgeRUtils::colSums3(x^2))
#       result <- crossprod(x) / outer(x_norm, x_norm)
#       return(result)
#     }

#     y <- as.matrix(y)
#     x_norm <- sqrt(SigBridgeRUtils::colSums3(x^2))
#     y_norm <- sqrt(SigBridgeRUtils::colSums3(y^2))

#     crossprod(x, y) / outer(x_norm, y_norm)
#   }

#   mat <- matrix(runif(10000), 100)

#   microbenchmark::microbenchmark(
#     r = raw_corCosine(mat),
#     cpp = corCosine(mat),
#     rfast = Rfast::Dist(mat, "cosine")
#   )
#   # Unit: microseconds
#   #  expr      min        lq      mean   median       uq       max neval cld
#   #     r 1187.276 1211.8030 1343.5449 1235.001 1307.913 10347.917   100  a
#   #   cpp  171.662  204.5435  250.2805  233.799  288.937   510.698   100   b

#   r_res <- raw_corCosine(mat)
#   cpp_res <- corCosine(mat)
#   rfast_res <- Rfast::Dist(mat, "cosine")

#   tol <- max(abs(r_res - cpp_res))
#   tol2 <- max(abs(r_res - rfast_res))

#   expect_lt(tol, 1e-10)

#   mat2 <- matrix(rnbinom(10000, 10, 0.1), 100)

#   microbenchmark::microbenchmark(
#     r = raw_corCosine(mat2),
#     cpp = corCosine(mat2)
#   )
#   # Unit: microseconds
#   #  expr      min        lq      mean   median       uq       max neval cld
#   #     r 1187.276 1211.8030 1343.5449 1235.001 1307.913 10347.917   100  a
#   #   cpp  171.662  204.5435  250.2805  233.799  288.937   510.698   100   b

#   r_res2 <- raw_corCosine(mat2)
#   cpp_res2 <- corCosine(mat2)

#   tol2 <- max(abs(r_res2 - cpp_res2))

#   expect_lt(tol2, 1e-6)
# })
