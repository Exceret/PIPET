# test_that("multiplication works", {
#   set.seed(123)
#   v <- runif(100, max = 1)

#   raw_c <- function(x) {
#     if (rlang::is_installed("cheapr")) {
#       cheapr::sqrt_(1 / 2 * (1 - x))
#     } else {
#       sqrt(1 / 2 * (1 - x))
#     }
#   }

#   microbenchmark::microbenchmark(
#     r = raw_c(v),
#     cpp = CorToDist(v)
#   )
#   # Unit: microseconds
#   #  expr   min     lq    mean median     uq    max neval cld
#   #     r 2.434 2.4540 2.91696  2.465 2.5250 41.568   100   a
#   #   cpp 2.165 2.1845 2.25538  2.214 2.2395  3.407   100   a

#   v2 <- runif(200, max = 1)
#   microbenchmark::microbenchmark(
#     r = raw_c(v2),
#     cpp = CorToDist(v2)
#   )
#   # Unit: microseconds
#   #  expr   min    lq    mean median    uq    max neval cld
#   #     r 2.775 2.815 3.32823 2.8405 2.931 44.834   100   a
#   #   cpp 2.514 2.545 2.63397 2.5845 2.655  3.516   100   a

#   r <- raw_c(v)
#   cpp <- CorToDist(v)

#   expect_lt(max(abs(r - cpp)), 1e-6)

#   r2 <- raw_c(v2)
#   cpp2 <- CorToDist(v2)

#   expect_lt(max(abs(r2 - cpp2)), 1e-6)
# })

# test_that("multiplication works", {
#   set.seed(123)
#   v <- runif(100, max = 1)

#   raw_d <- \(x) 1 - 2 * x^2

#   microbenchmark::microbenchmark(
#     r = raw_d(v),
#     cpp = DistToCor(v)
#   )
#   # Unit: microseconds
#   #  expr   min    lq     mean median    uq      max neval cld
#   #     r 1.122 1.282 43.91213  1.393 1.668 4206.008   100   a
#   #   cpp 2.063 2.094  2.91018  2.194 2.490   49.433   100   a

#   v2 <- runif(200, max = 1)
#   microbenchmark::microbenchmark(
#     r = raw_d(v2),
#     cpp = DistToCor(v2)
#   )

#   # Unit: microseconds
#   #  expr   min     lq    mean median    uq     max neval cld
#   #     r 1.733 1.9235 7.50448  2.039 2.184 529.082   100   a
#   #   cpp 2.975 3.0255 3.30972  3.076 3.151  22.453   100   a

#   r <- raw_d(v)
#   cpp <- DistToCor(v)

#   expect_lt(max(abs(r - cpp)), 1e-6)

#   r2 <- raw_d(v2)
#   cpp2 <- DistToCor(v2)

#   expect_lt(max(abs(r2 - cpp2)), 1e-6)
# })
