# test_Create_Markers.R

test_that("Create_Markers function works correctly", {
    # 创建模拟数据
    set.seed(123)

    # 生成模拟的bulk RNA-seq数据
    n_genes <- 1000
    n_samples <- 12

    # 创建有差异表达的基因
    bulk_data <- matrix(
        rnbinom(n_genes * n_samples, mu = 1000, size = 1),
        nrow = n_genes,
        ncol = n_samples
    )
    rownames(bulk_data) <- paste0("Gene", 1:n_genes)
    colnames(bulk_data) <- paste0("Sample", 1:n_samples)

    # 在前100个基因中创建差异表达
    bulk_data[1:50, 1:6] <- bulk_data[1:50, 1:6] * 5 # GroupA上调
    bulk_data[51:100, 7:12] <- bulk_data[51:100, 7:12] * 5 # GroupB上调

    # 创建colData
    colData <- data.frame(
        sample_id = colnames(bulk_data),
        group = rep(c(0, 1), each = 6),
        batch = rep(c("Batch1", "Batch2"), 6)
    )
    rownames(colData) <- colnames(bulk_data)

    # 测试1: 基本功能 - 两分类比较
    test_that("two-class comparison works", {
        markers <- Create_Markers(
            bulk_data = bulk_data,
            colData = colData,
            class_col = "group",
            lg2FC = 1,
            p.adjust = 0.05,
            show_lg2FC = T,
            verbose = T
        )

        expect_s3_class(markers, "data.frame")
        expect_true(nrow(markers) > 0)
        expect_true(all(c("genes", "class") %in% colnames(markers)))
        expect_true("log2FoldChange" %in% colnames(markers))
        expect_true(all(markers$class %in% c("0", "1")))
    })

    # 测试2: 不显示log2FC
    test_that("hide log2FC works", {
        markers <- Create_Markers(
            bulk_data = bulk_data,
            colData = colData,
            class_col = "group",
            show_lg2FC = FALSE,
            verbose = FALSE
        )

        expect_false("log2FoldChange" %in% colnames(markers))
        expect_true(all(c("genes", "class") %in% colnames(markers)))
    })

    # 测试3: 多分类比较
    test_that("multi-class comparison works", {
        # 创建三分类数据
        colData_multi <- colData
        colData_multi$group <- rep(c(0, 1, 2), each = 4)

        # 调整数据以包含三组差异
        bulk_data_multi <- bulk_data
        bulk_data_multi[101:150, 9:12] <- bulk_data_multi[101:150, 9:12] * 5 # GroupC上调

        markers <- Create_Markers(
            bulk_data = bulk_data_multi,
            colData = colData_multi,
            class_col = "group",
            verbose = FALSE
        )

        expect_s3_class(markers, "data.frame")
        expect_true(nrow(markers) > 0)
        expect_true(all(markers$class %in% c("0", "1", "2")))
    })

    # 测试4: 错误处理 - class_col不存在
    test_that("error when class_col not in colData", {
        expect_error(
            Create_Markers(
                bulk_data = bulk_data,
                colData = colData,
                class_col = "nonexistent_column",
                verbose = FALSE
            )
        )
    })

    # 测试5: 错误处理 - 只有单一类别
    test_that("error when only one class", {
        colData_single <- colData
        colData_single$group <- "SingleGroup"

        expect_error(
            Create_Markers(
                bulk_data = bulk_data,
                colData = colData_single,
                class_col = "group",
                verbose = FALSE
            )
        )
    })

    # 测试6: 不同的lg2FC阈值
    test_that("different lg2FC thresholds work", {
        markers_strict <- Create_Markers(
            bulk_data = bulk_data,
            colData = colData,
            class_col = "group",
            lg2FC = 2, # 更严格的阈值
            p.adjust = 0.05,
            verbose = FALSE
        )

        markers_lenient <- Create_Markers(
            bulk_data = bulk_data,
            colData = colData,
            class_col = "group",
            lg2FC = 0.5, # 更宽松的阈值
            p.adjust = 0.05,
            verbose = FALSE
        )

        # 宽松阈值应该找到更多marker基因
        expect_true(nrow(markers_lenient) >= nrow(markers_strict))
    })

    # 测试7: 不同的p.adjust阈值
    test_that("different p.adjust thresholds work", {
        markers_strict <- Create_Markers(
            bulk_data = bulk_data,
            colData = colData,
            class_col = "group",
            lg2FC = 1,
            p.adjust = 0.01, # 更严格的阈值
            verbose = FALSE
        )

        markers_lenient <- Create_Markers(
            bulk_data = bulk_data,
            colData = colData,
            class_col = "group",
            lg2FC = 1,
            p.adjust = 0.1, # 更宽松的阈值
            verbose = FALSE
        )

        # 宽松阈值应该找到更多marker基因
        expect_true(nrow(markers_lenient) >= nrow(markers_strict))
    })

    # 测试8: 因子变量处理
    test_that("factor variable handling works", {
        colData_factor <- colData
        colData_factor$group <- factor(
            colData_factor$group,
            levels = c(0, 1)
        )

        expect_error(Create_Markers(
            bulk_data = bulk_data,
            colData = colData_factor,
            class_col = "group",
            verbose = FALSE
        ))
    })

    # 测试9: 输出结构验证
    test_that("output structure is correct", {
        markers <- Create_Markers(
            bulk_data = bulk_data,
            colData = colData,
            class_col = "group",
            verbose = FALSE
        )

        # 检查必要列
        expect_true(all(c("genes", "class") %in% colnames(markers)))

        # 检查基因名称不是空的
        expect_false(any(is.na(markers$genes)))
        expect_false(any(markers$genes == ""))

        # 检查类别不是空的
        expect_false(any(is.na(markers$class)))
        expect_false(any(markers$class == ""))

        # 如果显示log2FC，检查数值
        if ("log2FoldChange" %in% colnames(markers)) {
            expect_true(is.numeric(markers$log2FoldChange))
            expect_false(any(is.na(markers$log2FoldChange)))
        }
    })

    # 测试10: 边界情况 - 极小数据集
    test_that("small dataset works", {
        # 创建极小的测试数据集
        small_bulk <- matrix(c(100, 10, 10, 100, 50, 50, 50, 50), nrow = 2)
        rownames(small_bulk) <- c("Gene1", "Gene2")
        colnames(small_bulk) <- c("S1", "S2", "S3", "S4")

        small_colData <- data.frame(
            group = c(0, 0, 1, 1)
        )
        rownames(small_colData) <- colnames(small_bulk)

        markers <- Create_Markers(
            bulk_data = small_bulk,
            colData = small_colData,
            class_col = "group",
            p.adjust = 0.1, # 放宽阈值以适应小数据集
            verbose = FALSE
        )

        expect_s3_class(markers, "data.frame")
        # 小数据集可能找不到显著的marker，但函数不应该出错
    })
})

# 性能测试
test_that("performance with larger datasets", {
    skip_on_cran() # 在CRAN上跳过耗时测试

    # 创建更大的测试数据集
    set.seed(123)
    n_genes_large <- 5000
    n_samples_large <- 50

    bulk_data_large <- matrix(
        rnbinom(n_genes_large * n_samples_large, mu = 1000, size = 1),
        nrow = n_genes_large,
        ncol = n_samples_large
    )
    rownames(bulk_data_large) <- paste0("Gene", 1:n_genes_large)
    colnames(bulk_data_large) <- paste0("Sample", 1:n_samples_large)

    colData_large <- data.frame(
        group = rep(c(0, 1), each = n_samples_large / 2)
    )
    rownames(colData_large) <- colnames(bulk_data_large)

    # 测试运行时间
    expect_time_lt <- function(expr, time_limit = 30) {
        start_time <- Sys.time()
        force(expr)
        end_time <- Sys.time()
        expect_lt(
            as.numeric(difftime(end_time, start_time, units = "secs")),
            time_limit
        )
    }

    expect_time_lt(
        {
            markers <- Create_Markers(
                bulk_data = bulk_data_large,
                colData = colData_large,
                class_col = "group",
                verbose = FALSE
            )
        },
        time_limit = 60
    ) # 60秒时间限制
})

# 测试参数验证
test_that("parameter validation works", {
    set.seed(123)
    bulk_data <- matrix(rnbinom(100 * 10, mu = 1000, size = 1), nrow = 100)
    rownames(bulk_data) <- paste0("Gene", 1:100)
    colnames(bulk_data) <- paste0("Sample", 1:10)

    colData <- data.frame(
        group = rep(c(0, 1), each = 5)
    )
    rownames(colData) <- colnames(bulk_data)

    # 测试无效的lg2FC
    expect_error(
        Create_Markers(
            bulk_data = bulk_data,
            colData = colData,
            class_col = "group",
            lg2FC = -1, # 无效的负值
            verbose = FALSE
        )
    )

    # 测试无效的p.adjust
    expect_error(
        Create_Markers(
            bulk_data = bulk_data,
            colData = colData,
            class_col = "group",
            p.adjust = 1.5, # 无效的值
            verbose = FALSE
        )
    )
})
