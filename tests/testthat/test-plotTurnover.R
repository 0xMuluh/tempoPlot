library(SummarizedExperiment)
library(ggplot2)

# s1: {a, b} -> {b, c} -> {b, c}; s2: {a} -> {}.
toy <- SummarizedExperiment(
    assays = list(counts = matrix(c(
        1, 0, 0, 2, 0,
        1, 1, 1, 0, 0,
        0, 3, 1, 0, 0), 3L, byrow = TRUE,
        dimnames = list(c("a", "b", "c"), paste0("c", 1:5)))),
    colData = DataFrame(id = c("s1", "s1", "s1", "s2", "s2"),
        t = c(0, 1, 2, 0, 3), arm = c("x", "x", "x", "y", "y")))
value <- function(df, panel, time, part) {
    df$.value[df$.panel == panel & df$.time == time & df$.part == part]
}

test_that("plotTurnover checks its arguments", {
    expect_error(plotTurnover(toy, "t", "id", assay.type = "none"),
        "'assay.type'")
    expect_error(plotTurnover(toy, "t", "id", type = "total"))
    expect_error(suppressWarnings(plotTurnover(toy[, c(1, 4)], "t", "id")),
        "two consecutive samples")
    expect_error(plotTurnover(toy, "t", "id", detection = NA), "'detection'")
})

test_that("plotTurnover counts gains and losses", {
    p <- plotTurnover(toy, "t", "id")
    # Losses are drawn down but labelled as counts.
    expect_false(any(grepl("-", layer_scales(p)$y$get_labels())))
    df <- p$data
    expect_identical(nrow(df), 6L)
    expect_equal(value(df, "s1", 1, "gained"), 1)
    expect_equal(value(df, "s1", 1, "lost"), -1)
    expect_equal(value(df, "s1", 2, "gained"), 0)
    expect_equal(value(df, "s2", 3, "lost"), -1)
})

test_that("plotTurnover computes turnover components", {
    df <- plotTurnover(toy, "t", "id", type = "components")$data
    # s1 0 -> 1: 3 features in either sample, 1 appeared, 1 disappeared.
    expect_equal(value(df, "s1", 1, "appearance"), 1 / 3)
    expect_equal(value(df, "s1", 1, "disappearance"), 1 / 3)
    expect_equal(value(df, "s2", 3, "disappearance"), 1)
})

test_that("plotTurnover averages subjects by group", {
    toy$t2 <- c(0, 1, 2, 0, 1)
    df <- plotTurnover(toy, "t2", "id", level = "group")$data
    # At time 1, s1 lost 1 and s2 lost 1.
    expect_equal(value(df, "", 1, "lost"), -1)
    expect_equal(value(df, "", 1, "gained"), 0.5)
    p <- plotTurnover(toy, "t2", "id", level = "group", group.by = "arm")
    expect_setequal(p$data$.panel, c("x", "y"))
    expect_s3_class(p$facet, "FacetWrap")
})

test_that("plotTurnover spans intervals and stacks components", {
    df <- plotTurnover(toy, "t", "id")$data
    s2 <- df[df$.panel == "s2", ]
    expect_equal(unique(s2$.from), 0)
    expect_equal(unique(s2$.time), 3)
    df <- plotTurnover(toy, "t", "id", type = "components")$data
    s1 <- df[df$.panel == "s1" & df$.time == 1, ]
    expect_equal(s1$.ymin, c(0, 1 / 3))
    expect_equal(s1$.ymax, c(1 / 3, 2 / 3))
    # Date and POSIXct intervals keep their class.
    toy$date <- as.Date("2020-01-01") + toy$t
    expect_s3_class(plotTurnover(toy, "date", "id")$data$.from, "Date")
    toy$hour <- as.POSIXct("2020-01-01", tz = "UTC") + toy$t * 3600
    df <- plotTurnover(toy, "hour", "id")$data
    expect_s3_class(df$.time, "POSIXct")
    expect_equal(unique(df$.time[df$.panel == "s2"]), toy$hour[5L])
})

test_that("plotTurnover counts presence above detection", {
    # Above 1: s1 {} -> {c} -> {}; s2 {a} -> {}.
    df <- plotTurnover(toy, "t", "id", detection = 1)$data
    expect_equal(value(df, "s1", 1, "gained"), 1)
    expect_equal(value(df, "s1", 2, "lost"), -1)
    expect_equal(value(df, "s2", 3, "lost"), -1)
})

test_that("plotTurnover skips intervals across gaps", {
    # s2 misses visit 2: its interval from visit 1 to 3 is skipped.
    toy$v <- factor(c(1, 2, 3, 1, 3), levels = 1:3)
    df <- plotTurnover(toy, "v", "id")$data
    expect_setequal(df$.panel, "s1")
    df <- plotTurnover(toy, "v", "id", max.gap = 2)$data
    expect_equal(value(df, "s2", 3 - 0.05, "lost"), -1)
    # Numeric times are skipped only beyond max.gap.
    df <- plotTurnover(toy, "t", "id", max.gap = 2)$data
    expect_setequal(df$.panel, "s1")
    expect_error(plotTurnover(toy[, 4:5], "v", "id"),
        "two consecutive samples")
})

test_that("plotTurnover fills parts with the categorical palette", {
    p <- plotTurnover(toy, "t", "id")
    ld <- layer_data(p, 1L)
    expect_identical(ld$fill[p$data$.part == "gained"][1L], .okabe_ito[1L])
    expect_identical(ld$fill[p$data$.part == "lost"][1L], .okabe_ito[2L])
    expect_identical(p$theme$legend.position, "bottom")
})
