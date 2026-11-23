library(SummarizedExperiment)
library(ggplot2)

effect <- data.frame(t = c(0, 1, 2, 0, 1, 2), estimate = c(0, 1, 2, 0, -1, -2),
    lower = c(-1, 0, 1, -1, -2, -3), upper = c(1, 2, 3, 1, 0, -1),
    contrast = rep(c("a", "b"), each = 3L))
# Two permutations of the effect of contrast a.
null <- data.frame(permutation = rep(1:2, each = 3L), t = c(0, 1, 2, 0, 1, 2),
    estimate = c(0.5, -0.5, 0, -0.5, 0.5, 0))

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
    expect_error(plotTimeEffect(effect, "t", p.col = "p"), "'p.col'")
    expect_error(plotTimeEffect(effect, "t", size.by = "n"), "'size.by'")
    expect_error(plotTimeEffect(effect, "t", null = effect), "'null'")
    expect_error(plotTimeEffect(effect, "t",
        null = transform(null, t = factor(t))), "type of 'x'")
    expect_error(plotTimeEffect(effect, "t", p.col = "estimate",
        threshold = 0), "'threshold'")
    expect_error(plotTimeEffect(effect, "t", p.col = "estimate",
        null = null), "'p.col'")
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

test_that("plotTimeEffect draws permutations behind the effect", {
    a <- effect[effect$contrast == "a", ]
    p <- plotTimeEffect(a, "t", null = null)
    geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1))
    expect_identical(unname(geoms),
        c("GeomHline", "GeomLine", "GeomRibbon", "GeomLine"))
    curves <- layer_data(p, 2L)
    expect_length(unique(curves$group), 2L)
    expect_equal(curves$y, null$estimate)
    expect_true(all(curves$colour == "grey50" & curves$alpha == 0.2))

    # With the group.by column, each null takes the colour of its group.
    grouped <- rbind(transform(null, contrast = "a"),
        transform(null, contrast = "b"))
    p <- plotTimeEffect(effect, "t", group.by = "contrast", null = grouped)
    curves <- layer_data(p, 2L)
    expect_length(unique(curves$group), 4L)
    expect_setequal(curves$colour, .okabe_ito[1:2])
})

test_that("plotTimeEffect draws p-values on a log scale", {
    effect$p <- c(0.01, 0.1, 0.001, 0.5, 0.05, 1)
    effect$n <- c(10, 20, 30, 10, 20, 30)
    p <- plotTimeEffect(effect, "t", group.by = "contrast", p.col = "p",
        threshold = 0.01, size.by = "n")
    geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1))
    # No interval: lower and upper belong to the estimate.
    expect_identical(unname(geoms), c("GeomHline", "GeomLine", "GeomPoint"))
    expect_equal(layer_data(p, 1L)$yintercept, -2)
    expect_identical(layer_data(p, 1L)$linetype, "dashed")
    expect_equal(sort(layer_data(p, 2L)$y), sort(log10(effect$p)))
    points <- layer_data(p, 3L)
    expect_identical(order(points$size), order(effect$n))
    expect_identical(p$labels$y, "p")
})
