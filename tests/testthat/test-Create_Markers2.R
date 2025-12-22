test_that("Create_Markers2 function works correctly", {
  set.seed(123)

  # 生成模拟的 bulk RNA-seq 原始 counts 数据
  n_genes <- 1000
  n_samples <- 12

  # 注意：R 的 rnbinom(mu, size) 参数化：mean = mu，dispersion = 1/size
  bulk_counts <- matrix(
    rnbinom(n = n_genes * n_samples, mu = 1000, size = 1),
    nrow = n_genes,
    ncol = n_samples,
    dimnames = list(
      paste0("Gene", 1:n_genes),
      paste0("Sample", 1:n_samples)
    )
  )

  # 创建差异表达基因
  bulk_counts[1:50, 1:6] <- bulk_counts[1:50, 1:6] * 5 # Group 0 (samples 1–6) 上调
  bulk_counts[51:100, 7:12] <- bulk_counts[51:100, 7:12] * 5 # Group 1 (samples 7–12) 上调

  # 构建 colData
  colData <- data.frame(
    sample_id = colnames(bulk_counts),
    group = rep(c(0, 1), each = 6),
    batch = rep(c("Batch1", "Batch2"), times = 6),
    row.names = colnames(bulk_counts)
  )

  # ✅ ⭐ 关键步骤：log2 转换（+1 避免 log(0)）
  bulk_log2 <- log2(bulk_counts + 1)

  # 测试1: 基本功能 - 两分类比较
  test_that("two-class comparison works", {
    markers <- Create_Markers2(
      bulk_data = bulk_log2,
      colData = colData,
      class_col = "group",
      lg2FC = 1,
      p.adjust = 0.05,
      show_lg2FC = TRUE,
      verbose = TRUE
    )

    expect_s3_class(markers, "data.frame")
    expect_true(nrow(markers) > 0)
    expect_true(all(c("genes", "class") %in% colnames(markers)))
    expect_true("log2FoldChange" %in% colnames(markers))
    expect_true(all(markers$class %in% c("group_0", "group_1")))
  })

  # 测试2: 多分类比较
  test_that("multi-class comparison works", {
    colData$group2 <- sample(c(0, 1, 2, 3), nrow(colData), replace = TRUE)

    expect_warning({
      markers <- Create_Markers2(
        bulk_data = bulk_log2,
        colData = colData,
        class_col = "group2",
        lg2FC = 1,
        p.adjust = 0.05,
        show_lg2FC = TRUE,
        verbose = TRUE
      )
    })
  })
})
