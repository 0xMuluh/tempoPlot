library(SummarizedExperiment)
library(ggplot2)

# Two features, two subjects; s1 at times 0, 1, 3 and s2 at 0, 1.
toy <- SummarizedExperiment(
    assays = list(counts = matrix(c(
        1, 2, 4, 0, 3,
        5, 0, 0, 1, 1), 2L, byrow = TRUE,
        dimnames = list(c("a", "b"), paste0("c", 1:5)))),
    colData = DataFrame(id = c("s1", "s1", "s1", "s2", "s2"),
        t = c(0, 1, 3, 0, 1), arm = c("x", "x", "x", "y", "y")))
cells <- function(p) p$data[order(p$data$.panel, p$data$.feature,
    p$data$.time), ]

test_that("plotTiles checks its arguments", {
    expect_error(plotTiles(toy, "t", "id", assay.type = "none"),
        "'assay.type'")
    expect_error(plotTiles(toy, "t", "id", features = "z"), "'features'")
    expect_error(plotTiles(toy, "t", "id", top = 0), "'top'")
    expect_error(plotTiles(toy, "t", "id", value = "effect"))
    expect_error(plotTiles(toy, "t", "id", group.by = "none"), "'group.by'")
})

test_that("plotTiles draws a tile per feature and sample", {
    p <- plotTiles(toy, "t", "id")
    df <- cells(p)
    expect_identical(nrow(df), 10L)
    expect_equal(df$.value[df$.panel == "s1" & df$.feature == "a"],
        c(1, 2, 4))
    # Edges of s1 at 0, 1, 3: -0.5, 0.5, 2, 4.
    s1 <- df[df$.panel == "s1" & df$.feature == "a", ]
    expect_equal(s1$.width, c(1, 1.5, 2))
    expect_equal(s1$.mid, c(0, 1.25, 3))
    expect_s3_class(p$facet, "FacetWrap")
    # Feature a has the larger mean, so it is on top.
    expect_identical(levels(df$.feature), c("b", "a"))
})

test_that("plotTiles selects features", {
    expect_identical(levels(plotTiles(toy, "t", "id", top = 1)$data$.feature),
        "a")
    p <- plotTiles(toy, "t", "id", features = "b")
    expect_identical(unique(as.character(p$data$.feature)), "b")
})

test_that("plotTiles computes change and occupancy", {
    df <- cells(plotTiles(toy, "t", "id", value = "change"))
    expect_equal(df$.value[df$.panel == "s1" & df$.feature == "a"],
        c(0, 1, 3))
    expect_equal(df$.value[df$.panel == "s2" & df$.feature == "b"],
        c(0, 0))
    df <- cells(plotTiles(toy, "t", "id", value = "occupancy"))
    expect_equal(df$.value[df$.panel == "s1" & df$.feature == "b"],
        c(1, 0, 0))
    # a in s1 is 1, 2, 4: present above 1 at the last two times.
    df <- cells(plotTiles(toy, "t", "id", value = "occupancy",
        detection = 1))
    expect_equal(df$.value[df$.panel == "s1" & df$.feature == "a"],
        c(0, 1, 1))
    expect_error(plotTiles(toy, "t", "id", value = "occupancy",
        detection = "a"), "'detection'")
})

test_that("plotTiles averages subjects with level = 'group'", {
    p <- plotTiles(toy, "t", "id", value = "occupancy", level = "group")
    df <- cells(p)
    # At time 0, a is present in s1 only; b in both.
    expect_equal(df$.value[df$.feature == "a" & df$.time == 0], 0.5)
    expect_equal(df$.value[df$.feature == "b" & df$.time == 0], 1)
    expect_identical(p$labels$fill, "occupancy")
    p <- plotTiles(toy, "t", "id", level = "group", group.by = "arm")
    expect_setequal(p$data$.panel, c("x", "y"))
})

test_that("plotTiles orders by clustering and handles visits and Dates", {
    p <- plotTiles(toy, "t", "id", order.by = "cluster")
    expect_setequal(levels(p$data$.feature), c("a", "b"))
    toy$v <- factor(toy$t, levels = c(0, 1, 2, 3))
    df <- plotTiles(toy, "v", "id")$data
    expect_true(all(df$.width == 1))
    toy$date <- as.Date("2020-01-01") + toy$t
    expect_s3_class(plotTiles(toy, "date", "id")$data$.mid, "Date")
    toy$hour <- as.POSIXct("2020-01-01", tz = "UTC") + toy$t * 3600
    expect_s3_class(plotTiles(toy, "hour", "id")$data$.mid, "POSIXct")
    # A tile reaches at most max.gap / 2, so the gap from 1 to 3 is blank.
    s1 <- cells(plotTiles(toy, "t", "id", max.gap = 1.5))
    s1 <- s1[s1$.panel == "s1" & s1$.feature == "a", ]
    expect_equal(s1$.lower, c(-0.5, 0.5, 2.25))
    expect_equal(s1$.upper, c(0.5, 1.75, 3.75))
    expect_error(plotTiles(toy, "t", "id", max.gap = 0), "'max.gap'")
    one <- toy[, toy$id == "s1"]
    p <- plotTiles(one, "t")
    expect_identical(nrow(p$data), 6L)
})

test_that("plotTiles draws effects from a results table", {
    res <- data.frame(
        feature = rep(c("a", "b", "c"), each = 2),
        visit = factor(rep(c("v1", "v2"), 3)),
        lfc = c(1, -2, 0.1, 0.2, 3, 4),
        q = c(0.01, 0.2, 0.5, 0.5, 0.04, 0.001))
    p <- plotTiles(res, "visit", effect = "lfc", p.col = "q")
    df <- p$data
    expect_identical(nrow(df), 6L)
    expect_equal(df$.value[df$.feature == "a"], c(1, -2))
    # Visits are evenly spaced.
    expect_true(all(df$.width == 1))
    dots <- p$layers[[2L]]$data
    expect_setequal(paste(dots$.feature, dots$.time),
        c("a 1", "c 1", "c 2"))
    expect_identical(p$labels$fill, "lfc")
    # The largest mean absolute effects: c, then a.
    p <- plotTiles(res, "visit", effect = "lfc", top = 2)
    expect_setequal(as.character(p$data$.feature), c("a", "c"))
    expect_length(p$layers, 1L)
})

test_that("plotTiles checks results tables", {
    res <- data.frame(feature = c("a", "a"), t = c(1, 1), lfc = 1:2)
    expect_error(plotTiles(res, "t", effect = "none"), "'effect'")
    expect_error(plotTiles(res, "t", effect = "lfc"), "one row per")
    res$t <- c(1, 2)
    expect_error(plotTiles(res, "t", effect = "lfc", threshold = 2),
        "'threshold'")
    expect_s3_class(plotTiles(res, "t", effect = "lfc"), "ggplot")
    res$t <- c("week2", "week10")
    expect_error(plotTiles(res, "t", effect = "lfc"), "factor")
})

test_that("plotTiles keeps the time type of results tables", {
    res <- data.frame(feature = "a", lfc = 1:3,
        day = as.Date("2020-01-01") + c(0, 1, 30))
    df <- plotTiles(res, "day", effect = "lfc", max.gap = 2)$data
    df <- df[order(df$.time), ]
    expect_s3_class(df$.mid, "Date")
    # Dates keep their spacing; the 29-day gap stays blank.
    expect_equal(df$.width, c(1, 1.5, 2))
    # POSIXct times are in hours.
    res$hour <- as.POSIXct("2020-01-01", tz = "UTC") + c(0, 1, 30) * 3600
    df <- plotTiles(res, "hour", effect = "lfc", max.gap = 2)$data
    df <- df[order(df$.time), ]
    expect_s3_class(df$.mid, "POSIXct")
    expect_equal(df$.width, c(1, 1.5, 2))
    res$visit <- factor(c("w2", "w10", "w20"), levels = c("w2", "w10", "w20"))
    df <- plotTiles(res, "visit", effect = "lfc")$data
    expect_equal(df$.time[order(df$.value)], 1:3)
})

test_that("plotTiles fills changes blue below zero, red above, NA grey", {
    # Change of b in s1 is 0, -5, -5 (the lowest); of a in s1 0, 1, 3.
    p <- plotTiles(toy, "t", "id", value = "change")
    fill <- ggplot_build(p)$plot$scales$get_scales("fill")
    expect_identical(col2rgb(fill$map(c(-5, 0))),
        col2rgb(c("#2166AC", "grey90")))
    rgb <- col2rgb(fill$map(3))
    expect_gt(rgb["red", 1L], rgb["blue", 1L])
    assay(toy)[1L, 2L] <- NA
    p <- plotTiles(toy, "t", "id")
    ld <- layer_data(p, 1L)
    expect_identical(ld$fill[is.na(p$data$.value)], "grey50")
    expect_identical(ld$fill[which.max(p$data$.value)], "#FDE725")
    expect_identical(p$theme$legend.position, "bottom")
})

test_that("plotTiles draws one panel per subject of a results table", {
    res <- data.frame(feature = rep(c("a", "b"), each = 4),
        id = rep(rep(c("s1", "s2"), each = 2), 2), t = rep(c(1, 2), 4),
        lfc = 1:8)
    p <- plotTiles(res, "t", effect = "lfc", group = "id")
    expect_setequal(p$data$.panel, c("s1", "s2"))
    expect_s3_class(p$facet, "FacetWrap")
    d <- p$data[p$data$.panel == "s2" & p$data$.feature == "b", ]
    expect_equal(d$.value[order(d$.time)], 7:8)
    expect_error(plotTiles(res, "t", effect = "lfc"), "one row per")
    expect_error(plotTiles(res, "t", effect = "lfc", group = "none"),
        "'group'")
})

test_that("plotTiles draws categorical effects", {
    # c has three cells with a strain, a two, b one.
    res <- data.frame(feature = rep(c("a", "b", "c"), each = 3),
        t = rep(1:3, 3),
        strain = c("s1", "s1", NA, NA, "s2", NA, "s1", "s2", "s2"),
        q = c(0.01, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.01))
    p <- plotTiles(res, "t", effect = "strain", p.col = "q")
    expect_identical(levels(p$data$.feature), c("b", "a", "c"))
    ld <- layer_data(p, 1L)
    expect_identical(unique(ld$fill[p$data$.value %in% "s1"]), "#0072B2")
    expect_identical(unique(ld$fill[p$data$.value %in% "s2"]), "#E69F00")
    expect_true(all(ld$fill[is.na(p$data$.value)] == "grey50"))
    expect_setequal(paste(p$layers[[2L]]$data$.feature,
        p$layers[[2L]]$data$.time), c("a 1", "c 3"))
    expect_identical(p$labels$fill, "strain")
    p <- plotTiles(res, "t", effect = "strain", top = 2)
    expect_setequal(as.character(p$data$.feature), c("a", "c"))
    # Factors keep their levels; logicals are categories too.
    res$strain <- factor(res$strain, levels = c("s2", "s1"))
    p <- plotTiles(res, "t", effect = "strain")
    expect_identical(unique(layer_data(p, 1L)$fill[p$data$.value %in%
        "s2"]), "#0072B2")
    res$shared <- res$strain == "s1"
    expect_s3_class(plotTiles(res, "t", effect = "shared"), "ggplot")
    res$day <- as.Date("2020-01-01") + res$t
    expect_error(plotTiles(res, "t", effect = "day"), "'effect'")
    res$q <- as.character(res$q)
    expect_error(plotTiles(res, "t", effect = "shared", p.col = "q"),
        "'p.col'")
})

test_that("plotTiles clusters categorical profiles by matching cells", {
    # a and c share their strains at every time; b differs everywhere.
    res <- data.frame(feature = rep(c("a", "b", "c"), each = 3),
        t = rep(1:3, 3),
        strain = c("x", "x", "y", "y", "y", "x", "x", "x", "y"))
    lv <- levels(plotTiles(res, "t", effect = "strain",
        order.by = "cluster")$data$.feature)
    expect_identical(abs(diff(match(c("a", "c"), lv))), 1L)
    expect_false(lv[2L] == "b")
})
