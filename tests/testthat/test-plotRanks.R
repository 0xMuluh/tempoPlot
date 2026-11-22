library(SummarizedExperiment)
library(ggplot2)

toy <- SummarizedExperiment(
    assays = list(counts = matrix(c(
        9, 1, 5,
        5, 9, 1,
        1, 5, 9,
        0, 0, 0), 4L, byrow = TRUE,
        dimnames = list(c("a", "b", "c", "d"), paste0("c", 1:3)))),
    colData = DataFrame(t = c(1, 2, 3)))

test_that("plotRanks checks its arguments", {
    expect_error(plotRanks(toy, "t", assay.type = "none"), "'assay.type'")
    expect_error(plotRanks(toy, "t", top = 0), "'top'")
    expect_error(plotRanks(toy, "t", coord = "flip"))
})

test_that("plotRanks ranks features within samples", {
    df <- plotRanks(toy, "t", top = 3)$data
    a <- df[df$.feature == "a", ]
    expect_equal(a$.rank, c(1, 3, 2))
    expect_setequal(as.character(df$.feature), c("a", "b", "c"))
    # Equal means keep row order; d is never in the top 3.
    expect_false("d" %in% df$.feature)
})

test_that("plotRanks ranks a single feature and draws polar clocks", {
    df <- plotRanks(toy["a", ], "t")$data
    expect_identical(levels(df$.feature), "a")
    expect_equal(df$.rank, c(1, 1, 1))
    p <- plotRanks(toy, "t", coord = "polar")
    expect_s3_class(p$coordinates, "CoordPolar")
})
