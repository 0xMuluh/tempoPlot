library(SummarizedExperiment)
library(ggplot2)

toy <- SummarizedExperiment(
    assays = list(counts = matrix(c(
        1, 3, 0, 4,
        10, 30, 20, 40,
        100, 300, 200, 400), 3L, byrow = TRUE,
        dimnames = list(c("a", "b", "c"), paste0("c", 1:4)))),
    colData = DataFrame(arm = c("x", "x", "y", "y")),
    rowData = DataFrame(score = c(1, 2, 3), label = c("p", "q", "r")))

test_that("plotFeatureScatter checks its arguments", {
    expect_error(plotFeatureScatter(toy, "median", "variance"), "'x.stat'")
    expect_error(plotFeatureScatter(toy, "mean", "label"), "'y.stat'")
    expect_error(plotFeatureScatter(toy, "mean", "variance",
        fit = "neutral"), "occupancy")
    # Relative abundances would give a flat curve at zero.
    assay(toy, "rel") <- t(t(assay(toy)) / colSums(assay(toy)))
    expect_error(plotFeatureScatter(toy, "mean", "occupancy",
        assay.type = "rel", fit = "neutral"), "counts")
})

test_that("plotFeatureScatter computes statistics per group", {
    df <- plotFeatureScatter(toy, "mean", "occupancy")$data
    expect_equal(df$.x, c(2, 25, 250))
    expect_equal(df$.y, c(0.75, 1, 1))
    df <- plotFeatureScatter(toy, "score", "cv", group.by = "arm")$data
    x <- df[df$.by == "x", ]
    expect_equal(x$.x, c(1, 2, 3))
    expect_equal(x$.y, rep(sd(c(1, 3)) / 2, 3))
})

test_that("plotFeatureScatter fits Taylor's law", {
    se <- SummarizedExperiment(
        assays = list(counts = outer(c(1, 10, 100), c(1, 3, 2, 4))))
    p <- plotFeatureScatter(se, "mean", "variance", fit = "power")
    curve <- p$layers[[2L]]$data
    # Variance grows with the square of the mean: slope 2.
    expect_match(unique(curve$.label), "slope 2.00")
    expect_identical(p$labels$subtitle, "slope 2.00")
})

test_that("plotFeatureScatter fits each group in its colour", {
    se <- SummarizedExperiment(
        assays = list(counts = rbind(a = c(1, 3, 0, 0),
            b = c(10, 30, 20, 40), c = c(100, 300, 200, 400),
            d = c(5, 6, 7, 9))),
        colData = DataFrame(arm = c("x", "x", "y", "y")))
    p <- plotFeatureScatter(se, "mean", "variance", fit = "power",
        group.by = "arm")
    # Feature a is zero in arm y, so it is not drawn there.
    expect_false(any(p$data$.by == "y" & p$data$.feature == "a"))
    expect_true(any(p$data$.by == "x" & p$data$.feature == "a"))
    ld <- layer_data(p, 2L)
    expect_identical(nrow(ld), 200L)
    expect_identical(sort(unique(ld$colour)),
        sort(unique(layer_data(p, 1L)$colour)))
    expect_match(p$labels$subtitle, "^x: slope .*; y: slope ")
})

test_that("plotFeatureScatter fits the neutral model", {
    set.seed(1)
    p <- c(0.3, 0.1, 0.03, 0.01, 0.003, 0.001)
    counts <- rmultinom(50, 2000, p)
    rownames(counts) <- letters[1:6]
    se <- SummarizedExperiment(assays = list(counts = counts))
    q <- plotFeatureScatter(se, "mean", "occupancy", fit = "neutral")
    curve <- q$layers[[2L]]$data
    expect_true(all(diff(curve$.y) >= 0))
    expect_match(unique(curve$.label), "^m = ")
})
