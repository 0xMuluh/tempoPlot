library(SummarizedExperiment)
library(ggplot2)

data("hitchip1006", package = "miaTime")
tse <- mia::transformAssay(hitchip1006, method = "hellinger")
lags <- c(1, 2, 3, 6, 9)
tse <- miaTime::addVariogram(tse, time.col = "time", group = "subject",
    assay.type = "hellinger", method = "euclidean", squared = TRUE,
    time.lag = lags)
tse <- miaTime::addVariogram(tse, time.col = "time", group = "subject",
    assay.type = "hellinger", method = "euclidean", squared = TRUE,
    time.lag = c(2, 6, 8), split.by = "sex", name = "by_sex")
v <- as.data.frame(metadata(tse)$variogram)
vs <- as.data.frame(metadata(tse)$by_sex)

test_that("plotLagProfile checks its arguments", {
    expect_error(plotLagProfile(hitchip1006), "addVariogram")
    expect_error(plotLagProfile(tse, name = "missing"), "addVariogram")
    bad <- tse
    metadata(bad)$variogram <- v[, c("lag", "estimate")]
    expect_error(plotLagProfile(bad), "lower and upper")
    expect_error(plotLagProfile(tse, colour.by = "sex"), "'colour.by'")
})

test_that("plotLagProfile draws supported estimates with intervals", {
    p <- plotLagProfile(tse, size = 2)
    ld <- layer_data(p, 2L)
    ok <- !is.na(v$estimate)
    expect_true(any(!ok))
    expect_equal(ld$x, v$lag[ok])
    expect_equal(ld$y, v$estimate[ok])
    expect_equal(ld$ymin, v$lower[ok])
    expect_equal(ld$ymax, v$upper[ok])
    expect_true(all(ld$size == 2))
})

test_that("plotLagProfile colours levels and leaves out ratios", {
    p <- plotLagProfile(tse, name = "by_sex", colour.by = "sex")
    ld <- layer_data(p, 2L)
    est <- vs[vs$type == "dissimilarity" & !is.na(vs$estimate), ]
    expect_true(any(vs$type == "ratio"))
    expect_identical(nrow(ld), nrow(est))
    expect_equal(sort(ld$y), sort(est$estimate))
    expect_length(unique(ld$colour), 2L)
    # Levels at the same lag are shifted apart.
    at2 <- ld[abs(ld$x - 2) < 0.5, ]
    expect_length(unique(at2$x), nrow(at2))
})
