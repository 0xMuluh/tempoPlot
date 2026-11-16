library(SummarizedExperiment)
library(SingleCellExperiment)
library(ggplot2)

# s1 moves east then north; s2 moves west, then stays.
toy <- SingleCellExperiment(
    assays = list(counts = matrix(1, 1L, 5L)),
    colData = DataFrame(id = c("s1", "s1", "s1", "s2", "s2"),
        t = c(0, 1, 2, 0, 1)))
reducedDim(toy, "PC") <- cbind(PC1 = c(0, 1, 1, 0, -2),
    PC2 = c(0, 0, 3, 0, 0))

test_that("plotDirection checks its arguments", {
    expect_error(plotDirection(toy, "t", NULL, "PC"), "'group'")
    expect_error(plotDirection(toy, "t", "id", "none"), "'dimred'")
    expect_error(plotDirection(toy, "t", "id", "PC", bins = 1), "'bins'")
})

test_that("plotDirection measures angles counterclockwise from east", {
    df <- plotDirection(toy, "t", "id", "PC")$data
    expect_equal(sort(df$.angle), c(0, 90, 180))
    expect_equal(sort(df$.length), c(1, 2, 3))
    expect_s3_class(plotDirection(toy, "t", "id", "PC")$coordinates,
        "CoordPolar")
})

test_that("plotDirection drops zero moves and splits by group.by", {
    reducedDim(toy, "PC")[5, ] <- 0
    df <- plotDirection(toy, "t", "id", "PC")$data
    expect_identical(nrow(df), 2L)
    toy$arm <- c("a", "b", "b", "b", "b")
    p <- plotDirection(toy, "t", "id", "PC", group.by = "arm")
    expect_setequal(p$data$arm, c("a", "b"))
})
