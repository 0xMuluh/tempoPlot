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
    one <- toy[, toy$id == "s1"]
    p <- plotTiles(one, "t")
    expect_identical(nrow(p$data), 6L)
})
