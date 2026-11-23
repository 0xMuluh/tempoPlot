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
        "CoordRadial")
})

test_that("plotDirection drops zero moves and splits by facet.by", {
    reducedDim(toy, "PC")[5, ] <- 0
    df <- plotDirection(toy, "t", "id", "PC")$data
    expect_identical(nrow(df), 2L)
    # Each move takes the arm of its later sample.
    toy$arm <- c("a", "b", "c", "b", "b")
    p <- plotDirection(toy, "t", "id", "PC", facet.by = "arm")
    expect_identical(p$data$arm[order(p$data$.angle)], c("b", "c"))
    expect_s3_class(p$facet, "FacetWrap")
})

test_that("plotDirection orders samples and skips gaps", {
    shuffled <- toy[, c(3, 5, 1, 4, 2)]
    df <- plotDirection(shuffled, "t", "id", "PC")$data
    expect_equal(sort(df$.angle), c(0, 90, 180))
    # s1 moves north from t = 1 to t = 5, beyond max.gap.
    toy$t[3] <- 5
    df <- plotDirection(toy, "t", "id", "PC", max.gap = 2)$data
    expect_equal(sort(df$.angle), c(0, 180))
})

test_that("plotDirection selects components and passes on ...", {
    reducedDim(toy, "PC") <- cbind(reducedDim(toy, "PC"), PC3 = 0)
    df <- plotDirection(toy, "t", "id", "PC", ncomponents = c(2, 1))$data
    expect_equal(sort(df$.angle), c(0, 90, 270))
    expect_error(plotDirection(toy, "t", "id", "PC", ncomponents = 3),
        "'ncomponents'")
    expect_true(all(layer_data(plotDirection(toy, "t", "id",
        "PC"))$fill == "#0072B2"))
    p <- plotDirection(toy, "t", "id", "PC", fill = "grey40")
    expect_true(all(layer_data(p)$fill == "grey40"))
})

test_that("plotDirection skips moves across visits without a sample", {
    # s1 is seen at visits 1, 2 and 4: its northward move skips visit 3.
    toy$v <- factor(c(1, 2, 4, 1, 2), levels = 1:4)
    df <- plotDirection(toy, "v", "id", "PC")$data
    expect_equal(sort(df$.angle), c(0, 180))
    df <- plotDirection(toy, "v", "id", "PC", max.gap = 2)$data
    expect_equal(sort(df$.angle), c(0, 90, 180))
})
