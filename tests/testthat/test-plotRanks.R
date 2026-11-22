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

test_that("plotRanks selects by mean rank and takes features", {
    # e has the largest mean count but ranks last in two samples.
    big <- rbind(assay(toy), e = c(100, 0, 0))
    big <- SummarizedExperiment(assays = list(counts = big),
        colData = colData(toy))
    df <- plotRanks(big, "t", top = 1)$data
    expect_identical(levels(df$.feature), "a")
    df <- plotRanks(big, "t", features = c("e", "b"))$data
    expect_identical(levels(df$.feature), c("e", "b"))
    expect_error(plotRanks(big, "t", features = "z"), "'features'")
    ld <- layer_data(plotRanks(toy, "t", top = 2), 1L)
    expect_identical(unique(ld$colour), .okabe_ito[1:2])
})

test_that("plotRanks breaks lines at skipped visits and at max.gap", {
    two <- toy[, rep(seq_len(ncol(toy)), 2L)]
    colnames(two) <- paste0("s", seq_len(ncol(two)))
    two$id <- rep(c("s1", "s2"), each = 3L)
    two$v <- factor(c("v1", "v2", "v4", "v1", "v2", "v3"),
        levels = paste0("v", 1:4))
    df <- plotRanks(two, "v", "id", top = 1)$data
    expect_equal(df$.segment[df$.panel == "s1"], c(1, 1, 2))
    expect_equal(df$.segment[df$.panel == "s2"], c(1, 1, 1))
    p <- plotRanks(two, "v", "id", top = 1)
    expect_identical(length(unique(layer_data(p, 1L)$group)), 3L)
    expect_true("v3" %in% layer_scales(p)$x$get_limits())
    # Numeric times break only with max.gap.
    two$t <- c(0, 1, 5, 0, 1, 2)
    expect_true(all(plotRanks(two, "t", "id", top = 1)$data$.segment == 1))
    df <- plotRanks(two, "t", "id", top = 1, max.gap = 2)$data
    expect_equal(df$.segment[df$.panel == "s1"], c(1, 1, 2))
})
