library(SummarizedExperiment)
library(ggplot2)

effect <- data.frame(t = c(0, 1, 2, 0, 1, 2), estimate = c(0, 1, 2, 0, -1, -2),
    lower = c(-1, 0, 1, -1, -2, -3), upper = c(1, 2, 3, 1, 0, -1),
    contrast = rep(c("a", "b"), each = 3L))

test_that("plotTimeEffect checks its arguments", {
    expect_error(plotTimeEffect(effect, "time"), "'time.col'")
    expect_error(plotTimeEffect(effect, "t", group.by = "arm"), "'group.by'")
    expect_error(plotTimeEffect(effect, "contrast"), "numeric, Date")
    hour <- transform(effect, t = as.POSIXct("2020-01-01", tz = "UTC") +
        t * 3600)
    expect_s3_class(layer_scales(plotTimeEffect(hour, "t"))$x,
        "ScaleContinuousDatetime")
    expect_error(plotTimeEffect(effect[, -2L], "t"), "estimate")
    expect_error(plotTimeEffect(effect[, -3L], "t"), "together")
    tse <- SummarizedExperiment(assays = list(counts = matrix(1, 1L, 1L)))
    expect_error(plotTimeEffect(tse, "effect", "t"), "'name'")
})

test_that("plotTimeEffect draws estimates and intervals against zero", {
    p <- plotTimeEffect(effect, "t", group.by = "contrast")
    geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1))
    expect_identical(unname(geoms), c("GeomHline", "GeomRibbon", "GeomLine"))
    expect_equal(layer_data(p, 1L)$yintercept, 0)
    expect_equal(sort(layer_data(p, 2L)$ymin), sort(effect$lower))
    line <- layer_data(p, 3L)
    expect_equal(sort(line$y), sort(effect$estimate))
    expect_identical(unique(line$colour), .okabe_ito[1:2])
    expect_identical(unique(layer_data(p, 2L)$fill), .okabe_ito[1:2])
    expect_identical(p$theme$legend.position, "bottom")
    expect_identical(p$labels$colour, "contrast")

    # Without an interval or groups: one line; ... goes to the line.
    p <- plotTimeEffect(effect[effect$contrast == "a", 1:2], "t",
        linewidth = 3)
    expect_length(p$layers, 2L)
    expect_true(all(layer_data(p, 2L)$linewidth == 3))
})

test_that("plotTimeEffect reads results from metadata(x)", {
    tse <- SummarizedExperiment(assays = list(counts = matrix(1, 1L, 1L)))
    S4Vectors::metadata(tse)$effect <- effect
    p <- plotTimeEffect(tse, "effect", "t", group.by = "contrast")
    expect_equal(sort(layer_data(p, 3L)$y), sort(effect$estimate))
})
