library(SummarizedExperiment)
library(ggplot2)

# a: two intervals, one of them a single time; b: one long interval that
# starts first; c: a single time.
iv <- data.frame(feature = c("a", "a", "b", "c"),
    start = c(5, 20, 0, 10), end = c(10, 20, 30, NA),
    dominant = c("x", "y", "x", "y"), method = c("m1", "m1", "m2", "m2"))

test_that("plotIntervals checks its arguments", {
    expect_error(plotIntervals(iv, feature.col = "none"), "'feature.col'")
    expect_error(plotIntervals(iv, colour.by = "none"), "'colour.by'")
    expect_error(plotIntervals(iv[, -2]), "'start' and 'end'")
    bad <- iv
    bad$end[1] <- 1
    expect_error(plotIntervals(bad), "before 'start'")
    bad <- iv
    bad$start <- as.Date("2020-01-01") + bad$start
    expect_error(plotIntervals(bad), "type of 'start'")
    expect_error(plotIntervals(SummarizedExperiment(), name = "none"),
        "'name'")
})

test_that("plotIntervals draws bars and points, earliest feature on top", {
    p <- plotIntervals(iv, colour.by = "dominant")
    # Top to bottom: b (starts at 0), a (5), c (10).
    expect_identical(levels(p$data$.feature), c("c", "a", "b"))
    bars <- layer_data(p, 1L)
    expect_equal(sort(bars$xend - bars$x), c(5, 30))
    points <- layer_data(p, 2L)
    expect_equal(sort(points$x), c(10, 20))
    expect_identical(p$labels$colour, "dominant")
    # The order holds although c has no bar.
    b <- ggplot_build(p)
    expect_identical(b$layout$panel_params[[1L]]$y$get_labels(),
        c("c", "a", "b"))
})

test_that("plotIntervals keeps the features with the longest intervals", {
    p <- plotIntervals(iv, top = 1)
    expect_identical(unique(as.character(p$data$.feature)), "b")
    p <- plotIntervals(iv, features = c("a", "c"))
    expect_setequal(as.character(p$data$.feature), c("a", "c"))
})

test_that("plotIntervals takes visits, dates, panels and metadata", {
    v <- iv
    lv <- paste0("v", 0:30)
    v$start <- factor(paste0("v", v$start), levels = lv)
    v$end <- factor(ifelse(is.na(v$end), NA, paste0("v", v$end)),
        levels = lv)
    p <- plotIntervals(v)
    expect_equal(sort(layer_data(p, 1L)$x), c(1, 6))
    d <- iv
    d$start <- as.Date("2020-01-01") + d$start
    d$end <- as.Date("2020-01-01") + d$end
    expect_s3_class(plotIntervals(d)$data$.start, "Date")
    p <- plotIntervals(iv, facet.by = "method")
    expect_s3_class(p$facet, "FacetWrap")
    se <- SummarizedExperiment()
    S4Vectors::metadata(se)$iv <- iv
    expect_identical(plotIntervals(se, name = "iv", top = 2)$data,
        plotIntervals(iv, top = 2)$data)
})
