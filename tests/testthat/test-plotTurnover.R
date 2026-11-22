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
        "more than one sample")
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
})
