library(SummarizedExperiment)
library(ggplot2)

res <- data.frame(
    lfc = c(2, -1.5, 0.1, 3, NA, -0.2),
    q = c(0.001, 0.01, 0.9, 0.04, 0.01, NA),
    row.names = paste0("f", 1:6))

test_that("plotVolcano checks its arguments", {
    expect_error(plotVolcano(res, "none", "q"), "'effect'")
    expect_error(plotVolcano(res, "lfc", "none"), "'p.col'")
    res$chr <- "a"
    expect_error(plotVolcano(res, "chr", "q"), "numeric")
    expect_error(plotVolcano(res, "lfc", "q", threshold = 2), "threshold")
    expect_error(plotVolcano(res, "lfc", "q", label.top = -1), "label.top")
    expect_error(plotVolcano(res[5:6, ], "lfc", "q"), "No feature")
    res$q[1] <- 2
    expect_error(plotVolcano(res, "lfc", "q"), "between 0 and 1")
})

test_that("plotVolcano draws effects against -log10 p", {
    p <- plotVolcano(res, "lfc", "q")
    # Rows with NA effect or p are dropped.
    expect_identical(p$data$.feature, paste0("f", 1:4))
    ld <- layer_data(p, 3L)
    expect_equal(ld$y, -log10(c(0.001, 0.01, 0.9, 0.04)))
    expect_identical(p$data$significant, c(TRUE, TRUE, FALSE, TRUE))
    expect_identical(ld$colour, c("#0072B2", "#0072B2", "grey80",
        "#0072B2"))
    expect_equal(layer_data(p, 2L)$yintercept, -log10(0.05))
})

test_that("plotVolcano labels the top significant features", {
    p <- plotVolcano(res, "lfc", "q", label.top = 2)
    expect_identical(p$data$.label, c("f1", "f2", "", ""))
    p <- plotVolcano(res, "lfc", "q", threshold = 0.005, label.top = 5)
    expect_identical(p$data$.label, c("f1", "", "", ""))
    p <- plotVolcano(res, "lfc", "q", label.top = 0)
    expect_true(all(p$data$.label == ""))
})

test_that("plotVolcano labels from a feature column", {
    long <- data.frame(feature = rownames(res), res, row.names = NULL)
    p <- plotVolcano(long, "lfc", "q", label.top = 2,
        feature.col = "feature")
    expect_identical(p$data$.label, c("f1", "f2", "", ""))
    expect_error(plotVolcano(long, "lfc", "q", feature.col = "none"),
        "'feature.col'")
})

test_that("plotVolcano reads rowData", {
    se <- SummarizedExperiment(
        assays = list(counts = matrix(1, 6L, 2L)),
        rowData = DataFrame(res))
    rownames(se) <- rownames(res)
    p <- plotVolcano(se, "lfc", "q")
    expect_identical(p$data, plotVolcano(res, "lfc", "q")$data)
    expect_error(plotVolcano(se, "none", "q"), "'effect'")
    rowData(se)$genus <- toupper(rownames(se))
    p <- plotVolcano(se, "lfc", "q", label.top = 1, feature.col = "genus")
    expect_identical(p$data$.label[1L], "F1")
})
