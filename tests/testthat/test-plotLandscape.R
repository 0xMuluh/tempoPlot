library(SummarizedExperiment)
library(ggplot2)

set.seed(1)
toy <- SummarizedExperiment(
    assays = list(counts = matrix(1, 1L, 60L)),
    colData = DataFrame(v = c(rnorm(30, -2, 0.3), rnorm(30, 2, 0.3)),
        arm = rep(c("a", "b"), 30)))

test_that("plotLandscape checks its arguments", {
    expect_error(plotLandscape(toy, "none"))
    expect_error(plotLandscape(toy, "v", scale = "log"))
    expect_error(plotLandscape(toy, "v", group.by = "none"), "'group.by'")
    expect_error(plotLandscape(toy[, 1], "v"), "two values")
})

test_that("plotLandscape draws the density and the potential", {
    p <- plotLandscape(toy, "v", bw = 0.3)
    d <- density(toy$v, bw = 0.3)
    expect_equal(p$data$.y, d$y)
    # Two states: the density peaks near -2 and 2 and dips at 0.
    mid <- p$data$.y[which.min(abs(p$data$.x))]
    expect_lt(mid, max(p$data$.y) / 10)
    q <- plotLandscape(toy, "v", scale = "potential", bw = 0.3)
    expect_equal(q$data$.y, -log(d$y[d$y > 0]))
    expect_identical(q$labels$y, "potential (-log density)")
})

test_that("plotLandscape passes bw and adjust to density, ... to the line", {
    p <- plotLandscape(toy, "v", bw = 0.3, adjust = 2, linetype = "dashed",
        linewidth = 2)
    expect_equal(p$data$.y, density(toy$v, bw = 0.3, adjust = 2)$y)
    ld <- layer_data(p, 1L)
    expect_true(all(ld$linetype == "dashed" & ld$linewidth == 2))
})

test_that("plotLandscape draws one curve per group", {
    p <- plotLandscape(toy, "v", group.by = "arm")
    expect_setequal(p$data$.by, c("a", "b"))
    expect_identical(p$labels$colour, "arm")
})
