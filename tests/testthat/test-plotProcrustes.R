library(SummarizedExperiment)
library(SingleCellExperiment)
library(ggplot2)

set.seed(1)
a <- matrix(rnorm(40), 20L)
theta <- pi / 5
rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2L)
toy <- SingleCellExperiment(assays = list(counts = matrix(1, 1L, 20L)),
    colData = DataFrame(g = rep(c("x", "y"), 10)))
reducedDim(toy, "A") <- a
reducedDim(toy, "B") <- 3 * a %*% rot + 5
reducedDim(toy, "C") <- matrix(rnorm(60), 20L)

test_that("plotProcrustes checks its arguments", {
    expect_error(plotProcrustes(toy, "none"), "'dimred'")
    expect_error(plotProcrustes(toy, "A", "C", ncomponents = 3),
        "'ncomponents'")
})

test_that("plotProcrustes undoes rotation, scaling and shift", {
    p <- plotProcrustes(toy, "A", "B", colour.by = "g")
    expect_equal(p$data$.x2, p$data$.x, tolerance = 1e-8)
    expect_match(p$labels$subtitle, "m2 = 0.000")
    expect_identical(p$labels$colour, "g")
})

test_that("plotProcrustes agrees with vegan", {
    skip_if_not_installed("vegan")
    p <- plotProcrustes(toy, "A", "C")
    v <- vegan::procrustes(a, reducedDim(toy, "C")[, 1:2],
        symmetric = TRUE)
    expect_match(p$labels$subtitle, sprintf("m2 = %.3f", v$ss))
})

test_that("plotProcrustes matches samples across experiments", {
    e1 <- toy[, 1:15]
    e2 <- toy[, 20:6]
    colnames(e1) <- paste0("s", 1:15)
    colnames(e2) <- paste0("s", 20:6)
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        list(one = e1, two = e2))
    p <- plotProcrustes(mae, "one", "two", dimred = "A")
    expect_identical(nrow(p$data), 10L)
    expect_equal(p$data$.x2, p$data$.x, tolerance = 1e-8)
    expect_error(plotProcrustes(mae, "one", "none", dimred = "A"),
        "'experiment2'")
})
