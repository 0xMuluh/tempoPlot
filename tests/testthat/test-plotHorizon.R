library(SummarizedExperiment)
library(ggplot2)

# Feature a in s1 deviates 0, +3, -6 from its median 0 at times 0-2.
toy <- SummarizedExperiment(
    assays = list(counts = matrix(c(
        0, 3, -6, 1, 1,
        2, 2, 2, 5, 9), 2L, byrow = TRUE,
        dimnames = list(c("a", "b"), paste0("c", 1:5)))),
    colData = DataFrame(id = c("s1", "s1", "s1", "s2", "s2"),
        t = c(0, 1, 2, 0, 1)))
band <- function(df, f, s, b) {
    df <- df[df$.feature == f & df$.panel == s & df$.band == b, ]
    df$.y[order(df$.time)]
}

test_that("plotHorizon checks its arguments", {
    expect_error(plotHorizon(toy, "t", "id", bands = 0), "'bands'")
    expect_error(plotHorizon(toy, "t", "id", centre = "max"))
    expect_error(plotHorizon(toy, "t", "id", assay.type = "none"),
        "'assay.type'")
})

test_that("plotHorizon folds deviations into bands", {
    df <- plotHorizon(toy, "t", "id", bands = 3)$data
    # Largest |deviation| of a is 6, so a band is 2 high.
    expect_equal(band(df, "a", "s1", "1"), c(0, 1, 0))
    expect_equal(band(df, "a", "s1", "2"), c(0, 0.5, 0))
    expect_equal(band(df, "a", "s1", "-3"), c(0, 0, 1))
    # A constant series has no deviation.
    expect_true(all(band(df, "b", "s1", "1") == 0))
    expect_identical(nlevels(df$.band), 6L)
})

test_that("plotHorizon centres on baseline and mean", {
    df <- plotHorizon(toy, "t", "id", bands = 1, centre = "baseline")$data
    # b in s2: 5, 9 from baseline 5 gives 0, +4 (the largest, one band).
    expect_equal(band(df, "b", "s2", "1"), c(0, 1))
    df <- plotHorizon(toy, "t", "id", bands = 1, centre = "mean")$data
    expect_equal(band(df, "b", "s2", "1"), c(0, 1))
    expect_equal(band(df, "b", "s2", "-1"), c(1, 0))
})

test_that("plotHorizon draws a strip per feature and panel", {
    p <- plotHorizon(toy, "t", "id")
    expect_s3_class(p$facet, "FacetGrid")
    # b has the larger mean, so it comes first (top).
    expect_identical(levels(p$data$.feature), c("b", "a"))
    expect_identical(class(p$layers[[1L]]$geom)[1L], "GeomRibbon")
    p <- plotHorizon(toy, "t", "id", level = "group")
    expect_setequal(p$data$.panel, "")
    toy$date <- as.Date("2020-01-01") + toy$t
    expect_s3_class(plotHorizon(toy, "date", "id")$data$.time, "Date")
})
