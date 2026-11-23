library(SummarizedExperiment)
library(ggplot2)

# Subject s1 at times 0-2, s2 at times 0-1. Means over the samples at each
# time: a 2, 3, 3; b 5, 5, 5 (constant); c 1, 2, 4. a and c form cluster
# x, b cluster y.
toy <- SummarizedExperiment(
    assays = list(counts = matrix(c(
        1, 2, 3, 3, 4,
        5, 5, 5, 5, 5,
        0, 2, 4, 2, 2), 3L, byrow = TRUE,
        dimnames = list(c("a", "b", "c"), paste0("c", 1:5)))),
    colData = DataFrame(id = c("s1", "s1", "s1", "s2", "s2"),
        t = c(0, 1, 2, 0, 1)),
    rowData = DataFrame(cluster = c("x", "y", "x")))
values <- function(p, f) {
    df <- p$data[p$data$.feature == f, ]
    df$.value[order(df$.time)]
}
cluster_mean <- function(p, k) {
    df <- p$layers[[2L]]$data
    df <- df[startsWith(as.character(df$.cluster), k), ]
    df$.value[order(df$.time)]
}

test_that("plotProfiles checks its arguments", {
    expect_error(plotProfiles(toy, "t", "id", clust.col = "none"),
        "'clust.col' must be one column name of rowData")
    expect_error(plotProfiles(toy, "t", "id", clust.col = NULL),
        "'clust.col'")
    expect_error(plotProfiles(toy, "t", "id", scale = "yes"), "'scale'")
    expect_error(plotProfiles(toy, "none", "id"), "'time.col'")
    expect_error(plotProfiles(toy, "t", "id", events = data.frame(x = 1)),
        "'events'")
})

test_that("plotProfiles draws mean values over samples at each time", {
    p <- plotProfiles(toy, "t", "id", scale = FALSE)
    expect_equal(values(p, "a"), c(2, 3, 3))
    expect_equal(values(p, "c"), c(1, 2, 4))
    # The bold line of x is the mean of a and c.
    expect_equal(cluster_mean(p, "x"), c(1.5, 2.5, 3.5))
    expect_identical(p$labels$y, "Mean counts")
})

test_that("plotProfiles scales each feature's mean trajectory", {
    p <- plotProfiles(toy, "t", "id")
    expect_equal(values(p, "a"), as.vector(scale(c(2, 3, 3))))
    # A constant trajectory is drawn at 0.
    expect_equal(values(p, "b"), c(0, 0, 0))
    expect_equal(cluster_mean(p, "x"), (as.vector(scale(c(2, 3, 3))) +
        as.vector(scale(c(1, 2, 4)))) / 2)
    expect_identical(p$labels$y, "Scaled mean counts")
})

test_that("plotProfiles draws a panel per cluster with its size", {
    p <- plotProfiles(toy, "t", "id")
    expect_s3_class(p$facet, "FacetWrap")
    expect_identical(levels(p$data$.cluster), c("x (n = 2)", "y (n = 1)"))
    expect_identical(class(p$layers[[1L]]$geom)[1L], "GeomLine")
    expect_identical(p$layers[[2L]]$aes_params$colour, .okabe_ito[[1L]])
    # Features without a cluster are left out; top keeps the largest.
    rowData(toy)$cluster[2L] <- NA
    expect_identical(levels(plotProfiles(toy, "t", "id")$data$.cluster),
        "x (n = 2)")
    p <- plotProfiles(toy, "t", "id", top = 1)
    expect_identical(unique(as.character(p$data$.feature)), "a")
})

test_that("plotProfiles orders panels by factor levels or sorted values", {
    rowData(toy)$cluster <- factor(c("x", "y", "x"), levels = c("y", "x"))
    p <- plotProfiles(toy, "t", "id")
    expect_identical(levels(p$data$.cluster), c("y (n = 1)", "x (n = 2)"))
    rowData(toy)$cluster <- c(10, 2, 10)
    p <- plotProfiles(toy, "t", "id")
    expect_identical(levels(p$data$.cluster), c("2 (n = 1)", "10 (n = 2)"))
})

test_that("plotProfiles breaks lines at gaps and keeps the time class", {
    toy$v <- factor(c(0, 1, 3, 0, 1), levels = 0:3)
    p <- plotProfiles(toy, "v", "id")
    expect_equal(p$data$.segment[p$data$.feature == "a"], c(1, 1, 2))
    expect_equal(p$layers[[2L]]$data$.segment[
        startsWith(as.character(p$layers[[2L]]$data$.cluster), "x")],
        c(1, 1, 2))
    p <- plotProfiles(toy, "t", "id", max.gap = 0.5)
    expect_equal(p$data$.segment[p$data$.feature == "a"], 1:3)
    toy$date <- as.Date("2020-01-01") + toy$t
    expect_s3_class(plotProfiles(toy, "date", "id")$data$.time, "Date")
})

test_that("plotProfiles marks events in every panel", {
    events <- data.frame(start = c(0, 1), end = c(0.5, NA),
        label = c("batch", "added"))
    p <- plotProfiles(toy, "t", "id", events = events)
    geoms <- unname(vapply(p$layers, function(l) class(l$geom)[1L],
        character(1)))
    expect_identical(geoms,
        c("GeomRect", "GeomText", "GeomLine", "GeomLine", "GeomVline"))
    expect_equal(p$layers[[5L]]$data$.start, 1)
})
