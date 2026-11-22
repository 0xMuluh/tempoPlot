library(SummarizedExperiment)
library(ggplot2)

data("hitchip1006", package = "miaTime")
x <- hitchip1006[, hitchip1006$subject %in%
    names(which(table(hitchip1006$subject) >= 2))]
gaps <- unlist(tapply(x$time, x$subject, function(t) diff(sort(t))))

test_that("plotGaps checks its arguments", {
    expect_error(plotGaps(x, "time", group = NULL), "'group'")
    expect_error(plotGaps(x, "time", "subject", group.by = "none"),
        "'group.by'")
    expect_error(plotGaps(x, "time", "subject", layout = "violin"))
    one <- hitchip1006[, !duplicated(hitchip1006$subject)]
    expect_error(suppressWarnings(plotGaps(one, "time", "subject")),
        "more than one sample")
})

test_that("plotGaps draws the intervals within subjects", {
    p <- plotGaps(x, "time", "subject")
    expect_equal(sort(p$data$interval), sort(unname(gaps)))
    expect_identical(class(p$layers[[1L]]$geom)[1L], "GeomBar")
    p <- plotGaps(x, "time", "subject", layout = "density")
    expect_identical(class(p$layers[[1L]]$geom)[1L], "GeomDensity")
})

test_that("plotGaps splits intervals by group.by of the later sample", {
    p <- plotGaps(x, "time", "subject", group.by = "sex", layout = "box")
    expect_identical(class(p$layers[[1L]]$geom)[1L], "GeomBoxplot")
    expect_setequal(as.character(p$data$sex), unique(x$sex))
    expect_length(unique(layer_data(p)$x), length(unique(x$sex)))
    p <- plotGaps(x, "time", "subject", group.by = "sex")
    expect_setequal(layer_data(p)$fill,
        c(.okabe_ito[1:2], if (anyNA(x$sex)) "grey50"))
    expect_identical(p$theme$legend.position, "bottom")
    p <- plotGaps(x, "time", "subject", group.by = "sex",
        layout = "density")
    expect_setequal(layer_data(p)$colour,
        c(.okabe_ito[1:2], if (anyNA(x$sex)) "grey50"))
})

test_that("plotGaps measures intervals between consecutive samples", {
    toy <- SummarizedExperiment(
        assays = list(counts = matrix(1, 1L, 5L)),
        colData = DataFrame(id = c("b", "a", "a", "b", "a"),
            t = c(10, 5, 0, 0, 1), arm = c("x", "y", "x", "x", "z")))
    df <- plotGaps(toy, "t", "id", group.by = "arm")$data
    # a: 0 -> 1 -> 5, b: 0 -> 10; arm of the later sample.
    expect_identical(df$interval, c(1, 4, 10))
    expect_identical(df$arm, c("z", "y", "x"))
})

test_that("plotGaps counts visits in level steps and Dates in days", {
    x$visit <- factor(x$time)
    p <- plotGaps(x, "visit", "subject")
    expect_true(all(p$data$interval == round(p$data$interval)))
    x$date <- as.Date("2000-01-01") + round(x$time * 365)
    p <- plotGaps(x, "date", "subject")
    expect_equal(
        sort(p$data$interval),
        sort(unname(unlist(tapply(x$date, x$subject,
            function(t) diff(as.numeric(sort(t)))))))
    )
})
