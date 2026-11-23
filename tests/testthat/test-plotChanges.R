library(SummarizedExperiment)
library(ggplot2)

# s1 at t = 0, 1, 2: a doubles twice; b is 0 in the middle, so it has no
# change. s2 at t = 0, 1, 3: a falls to a quarter, then doubles over two
# time units; b stays at 1, then falls to 0.
toy <- SummarizedExperiment(
    assays = list(counts = matrix(c(
        1, 2, 4, 4, 1, 2,
        4, 0, 1, 1, 1, 0), 2L, byrow = TRUE,
        dimnames = list(c("a", "b"), paste0("c", 1:6)))),
    colData = DataFrame(id = rep(c("s1", "s2"), each = 3),
        t = c(0, 1, 2, 0, 1, 3), arm = rep(c("x", "y"), each = 3)))
changes <- function(...) {
    .change_table(toy, "t", "id", "counts", group.by = NULL, ...)
}

test_that("plotChanges checks its arguments", {
    expect_error(plotChanges(toy, "t", "id", assay.type = "none"),
        "'assay.type'")
    expect_error(plotChanges(toy, "t", "id", fit = "normal"))
    expect_error(plotChanges(toy, "t", "id", bins = 0), "'bins'")
    expect_error(plotChanges(toy, "t", "id", detection = 5),
        "No feature is above 'detection'")
})

test_that("plotChanges takes log ratios of consecutive samples", {
    df <- changes(detection = 0, max.gap = NULL)
    expect_equal(sort(df$.change),
        sort(log(c(2, 2, 1 / 4, 2, 1))))
    expect_equal(sort(df$.interval), c(1, 1, 1, 1, 2))
    # max.gap = 1 drops the interval of two time units of s2.
    df <- changes(detection = 0, max.gap = 1)
    expect_equal(sort(df$.change), sort(log(c(2, 2, 1 / 4, 1))))
    # Above detection 1, only a in s1 from 2 to 4 remains.
    df <- changes(detection = 1, max.gap = NULL)
    expect_equal(df$.change, log(2))
})

test_that("plotChanges draws densities and reports the intervals", {
    p <- plotChanges(toy, "t", "id", bins = 4)
    expect_identical(get_labs(p)$subtitle, "5 changes, t intervals 1 to 2")
    # Bin densities integrate to one.
    width <- diff(pretty(range(log(c(2, 1 / 4))), n = 4))[1L]
    expect_equal(sum(p$data$.y) * width, 1)
    expect_true(all(p$data$.y > 0))
    p <- plotChanges(toy, "t", "id", max.gap = 1)
    expect_identical(get_labs(p)$subtitle, "4 changes, t intervals 1")
})

test_that("plotChanges fits a Laplace distribution per group", {
    p <- plotChanges(toy, "t", "id", max.gap = 1, fit = "laplace")
    # Changes log 2, log 2, -log 4, 0: median log(2) / 2, scale the mean
    # absolute deviation from it.
    v <- log(c(2, 2, 1 / 4, 1))
    b <- mean(abs(v - log(2) / 2))
    # The curve starts at the smallest change, -log 4.
    curve <- layer_data(p, 2L)
    expect_equal(10^curve$y[1L],
        exp(-abs(-log(4) - log(2) / 2) / b) / (2 * b))
    expect_match(get_labs(p)$subtitle, sprintf("Laplace scale %.2f", b))
    p <- plotChanges(toy, "t", "id", group.by = "arm", fit = "laplace")
    expect_setequal(p$data$.by, c("x", "y"))
    # In x, both changes are log 2, so no Laplace distribution fits.
    expect_match(get_labs(p)$subtitle, "Laplace scale y [0-9.]+$")
    expect_identical(length(unique(layer_data(p, 1L)$colour)), 2L)
})
