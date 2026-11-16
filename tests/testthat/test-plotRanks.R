library(SummarizedExperiment)
library(ggplot2)

toy <- SummarizedExperiment(
    assays = list(counts = matrix(c(
        9, 1, 5,
        5, 9, 1,
        1, 5, 9,
        0, 0, 0), 4L, byrow = TRUE,
        dimnames = list(c("a", "b", "c", "d"), paste0("c", 1:3)))),
    colData = DataFrame(t = c(1, 2, 3)),
    rowData = DataFrame(genus = c("g1", "g1", "g2", "g3")))

test_that("plotRanks checks its arguments", {
    expect_error(plotRanks(toy, "t", assay.type = "none"), "'assay.type'")
    expect_error(plotRanks(toy, "t", top = 0), "'top'")
    expect_error(plotRanks(toy, "t", rank = "none"), "'rank'")
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

test_that("plotRanks agglomerates by rank and draws polar clocks", {
    df <- plotRanks(toy, "t", top = 1, rank = "genus")$data
    # g1 = a + b = 14, 10, 6; g2 = c = 1, 5, 9.
    expect_identical(levels(df$.feature), "g1")
    expect_equal(df$.rank, c(1, 1, 2))
    p <- plotRanks(toy, "t", coord = "polar")
    expect_s3_class(p$coordinates, "CoordPolar")
})
