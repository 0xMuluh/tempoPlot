library(SummarizedExperiment)
library(ggplot2)

data("hitchip1006", package = "miaTime")

# Subjects with 3 or more samples. All start at time 0, so age at each
# sample gives start times that differ between subjects.
x <- hitchip1006[, hitchip1006$subject %in%
    names(which(table(hitchip1006$subject) >= 3))]
x$age.time <- x$age + x$time
n_subjects <- length(unique(x$subject))

test_that("plotSwimmer checks its arguments", {
    expect_error(plotSwimmer(x, "time", group = NULL), "'group'")
    expect_error(
        plotSwimmer(x, "time", "subject", colour.by = "missing"),
        "'colour.by'"
    )
    expect_error(plotSwimmer(x, "time", "subject", order.by = "last"))
    x_chr <- x
    x_chr$time <- as.character(x_chr$time)
    expect_error(plotSwimmer(x_chr, "time", "subject"), "factor")
})

test_that("plotSwimmer draws a point per sample and a lane per subject", {
    p <- plotSwimmer(x, "time", "subject")
    expect_s3_class(p, "ggplot")
    expect_identical(nrow(layer_data(p, 2L)), ncol(x))
    expect_length(unique(layer_data(p, 1L)$group), n_subjects)

    colData(x)[["sampling time"]] <- x$time
    p <- plotSwimmer(x, "sampling time", "subject", size = 4)
    expect_true(all(layer_data(p, 2L)$size == 4))
})

test_that("plotSwimmer orders subjects by first time, first at the top", {
    ld <- layer_data(plotSwimmer(x, "age.time", "subject"), 2L)
    # tapply() returns rows bottom to top.
    start <- rev(tapply(ld$x, ld$y, min))
    expect_gt(length(unique(start)), 1L)
    expect_true(all(diff(start) >= 0))
})

test_that("plotSwimmer orders subjects by number of samples", {
    ld <- layer_data(plotSwimmer(x, "time", "subject", order.by = "n"), 2L)
    n <- rev(as.vector(table(ld$y)))
    expect_gt(length(unique(n)), 1L)
    expect_true(all(diff(n) <= 0))
})

test_that("plotSwimmer breaks lanes at gaps larger than max.gap", {
    gaps <- unlist(tapply(x$time, x$subject, function(t) diff(sort(t))))
    p <- plotSwimmer(x, "time", "subject", max.gap = 5)
    expect_length(
        unique(layer_data(p, 1L)$group),
        n_subjects + sum(gaps > 5)
    )
})

test_that("plotSwimmer colours points by colour.by", {
    p <- plotSwimmer(x, "time", "subject", colour.by = "bmi_group")
    expect_length(
        unique(layer_data(p, 2L)$colour),
        length(unique(x$bmi_group))
    )
    p <- plotSwimmer(x, "time", "subject")
    expect_length(unique(layer_data(p, 2L)$colour), 1L)
})

test_that("plotSwimmer keeps subjects with a single sample", {
    expect_warning(
        p <- plotSwimmer(hitchip1006, "time", "subject"),
        "single sample"
    )
    ld <- layer_data(p, 2L)
    expect_identical(nrow(ld), ncol(hitchip1006))
    expect_length(unique(ld$y), length(unique(hitchip1006$subject)))
})

test_that("plotSwimmer draws visits in level order with unused levels", {
    x$visit <- factor(x$time, levels = c(-1, sort(unique(x$time))))
    ld <- layer_data(plotSwimmer(x, "visit", "subject"), 2L)
    expect_identical(
        sort(unique(as.integer(ld$x))),
        sort(unique(as.integer(x$visit)))
    )
})

test_that("plotSwimmer draws Date times on a continuous axis", {
    x$date <- as.Date("2000-01-01") + round(x$time * 30)
    p <- plotSwimmer(x, "date", "subject")
    expect_s3_class(layer_scales(p)$x, "ScaleContinuousDate")
    expect_equal(sort(layer_data(p, 2L)$x), sort(as.numeric(x$date)))
})
