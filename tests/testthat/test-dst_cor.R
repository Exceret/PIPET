test_that("multiplication works", {
  set.seed(123)
  v <- runif(100, max = 1)

  microbenchmark::microbenchmark(
    r = CorToDist(v),
    cpp = cor_to_dist_inplace(v)
  )
  # Unit: microseconds
  #  expr     min       lq      mean   median     uq      max neval cld
  #     r 515.297 523.1110 543.39359 531.3615 535.88 1864.327   100  a
  #   cpp   2.284   2.4545   4.05042   2.7550   5.48   19.026   100   b

  v2 <- runif(500, max = 1)
  microbenchmark::microbenchmark(
    r = CorToDist(v2),
    cpp = cor_to_dist_inplace(v2)
  )

  #     Unit: microseconds
  #  expr     min       lq      mean  median       uq      max neval cld
  #     r 517.801 522.7555 553.40572 533.470 538.0895 2337.433   100  a
  #   cpp   3.577   3.7870   5.80993   6.397   6.8980   32.751   100   b

  r <- CorToDist(v)
  cpp <- cor_to_dist_inplace(v)

  expect_lt(max(abs(r - cpp)), 1e-6)

  r2 <- CorToDist(v2)
  cpp2 <- cor_to_dist_inplace(v2)

  expect_lt(max(abs(r2 - cpp2)), 1e-6)
})

test_that("multiplication works", {
  set.seed(123)
  v <- runif(100, max = 1)

  microbenchmark::microbenchmark(
    r = DistToCor(v),
    cpp = dist_to_cor_inplace(v)
  )
  # Unit: microseconds
  #  expr   min    lq     mean median    uq      max neval cld
  #     r 1.122 1.282 43.91213  1.393 1.668 4206.008   100   a
  #   cpp 2.063 2.094  2.91018  2.194 2.490   49.433   100   a

  v2 <- runif(500, max = 1)
  microbenchmark::microbenchmark(
    r = DistToCor(v2),
    cpp = dist_to_cor_inplace(v2)
  )

  # Unit: microseconds
  #  expr   min     lq    mean median    uq     max neval cld
  #     r 1.733 1.9235 7.50448  2.039 2.184 529.082   100   a
  #   cpp 2.975 3.0255 3.30972  3.076 3.151  22.453   100   a

  r <- DistToCor(v)
  cpp <- dist_to_cor_inplace(v)

  expect_lt(max(abs(r - cpp)), 1e-6)

  r2 <- DistToCor(v2)
  cpp2 <- dist_to_cor_inplace(v2)

  expect_lt(max(abs(r2 - cpp2)), 1e-6)
})
