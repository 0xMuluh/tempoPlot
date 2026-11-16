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
    expect_length(unique(layer_data(p)$fill), length(unique(x$sex)))
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
