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

test_that("plotLagProfile draws one line per level and leaves out ratios", {
    p <- plotLagProfile(tse, name = "by_sex")
    expect_length(unique(layer_data(p, 1L)$group), 2L)
    ld <- layer_data(p, 2L)
    est <- vs[vs$type == "dissimilarity" & !is.na(vs$estimate), ]
    expect_true(any(vs$type == "ratio"))
    expect_identical(nrow(ld), nrow(est))
    expect_equal(sort(ld$y), sort(est$estimate))
    expect_setequal(ld$colour, c("#0072B2", "#E69F00"))
    # Levels at the same lag are shifted apart.
    at2 <- ld[abs(ld$x - 2) < 0.5, ]
    expect_length(unique(at2$x), nrow(at2))
})

test_that("plotLagProfile draws ratios between split.by levels", {
    vg <- data.frame(lag = rep(c(2, 6), each = 3),
        sex = rep(c("female", "male", "male / female"), 2),
        type = rep(c("dissimilarity", "dissimilarity", "ratio"), 2),
        estimate = c(0.2, 0.3, 1.5, 0.4, 0.2, 0.5),
        lower = c(0.1, 0.2, 1.1, 0.3, 0.1, 0.3),
        upper = c(0.3, 0.4, 2.0, 0.5, 0.3, 0.9))
    toy <- SummarizedExperiment(assays = list(counts = matrix(1, 1L, 1L)))
    S4Vectors::metadata(toy)$v <- vg
    p <- plotLagProfile(toy, "v", type = "ratio")
    expect_equal(p$data$estimate, c(1.5, 0.5))
    expect_identical(p$labels$y, "Ratio of D(L)")
    expect_s3_class(p$scales$get_scales("y")$trans, "transform")
    p <- plotLagProfile(toy, "v")
    expect_identical(nrow(p$data), 4L)
    # Two comparisons with the first level are two lines.
    three <- rbind(vg, transform(vg[vg$type == "ratio", ],
        sex = "other / female"))
    S4Vectors::metadata(toy)$v3 <- three
    p <- plotLagProfile(toy, "v3", type = "ratio")
    expect_length(unique(layer_data(p, 1L)$group), 2L)
    expect_identical(p$labels$colour, "sex")
    S4Vectors::metadata(toy)$w <- vg[vg$type == "dissimilarity", -3L]
    expect_error(plotLagProfile(toy, "w", type = "ratio"), "split.by")
})

test_that("plotLagProfile reads a variogram table directly", {
    expect_equal(plotLagProfile(v)$data, plotLagProfile(tse)$data)
    expect_error(plotLagProfile(v, facet.by = "subject"),
        "'facet.by' must be NULL or one column name of 'x'")
    expect_error(plotLagProfile(v[, c("lag", "estimate")]), "acf")
})

# Autocorrelations of features f1 and f2 in subjects a and b at lags of 1
# to 3 steps of 2 days. Subject b has 4, 1 and 0 pairs of samples: no
# estimate at lag 3.
acfs <- data.frame(
    FeatureID = rep(rep(c("f1", "f2"), each = 3L), 2L),
    subject = rep(c("a", "b"), each = 6L),
    lag = rep(1:3, 4L), time.lag = rep(c(2, 4, 6), 4L),
    acf = c(0.8, 0.5, 0.2, 0.1, -0.1, 0, 0.6, 0.3, NA, 0.2, 0, NA),
    n.pairs = c(16, 9, 4, 16, 9, 4, 4, 1, 0, 4, 1, 0))

test_that("plotLagProfile draws a line per feature and subject", {
    p <- plotLagProfile(acfs, group = "subject", linewidth = 2)
    expect_true(all(layer_data(p, 1L)$yintercept == 0))
    ld <- layer_data(p, 3L)
    ok <- !is.na(acfs$acf)
    # x is the lag in days; lags without an estimate are dropped.
    expect_equal(ld$x, acfs$time.lag[ok])
    expect_equal(ld$y, acfs$acf[ok])
    expect_length(unique(ld$group), 4L)
    expect_true(all(ld$colour == "#0072B2"))
    expect_true(all(ld$linewidth == 2))
    # One panel per subject.
    expect_equal(as.integer(ld$PANEL), rep(1:2, c(6L, 4L)))
    # One band per subject, +-1.96 / sqrt(n.pairs), lags with pairs only.
    band <- layer_data(p, 2L)
    expect_equal(band$ymax, 1.96 / c(4, 3, 2, 2, 1))
    expect_equal(band$ymin, -band$ymax)
    expect_equal(as.integer(band$PANEL), c(1L, 1L, 1L, 2L, 2L))
})

test_that("plotLagProfile colours, facets and reads autocorrelations", {
    p <- plotLagProfile(acfs, group = "subject", colour.by = "FeatureID",
        facet.by = NULL)
    ld <- layer_data(p, 3L)
    expect_setequal(ld$colour, c("#0072B2", "#E69F00"))
    expect_identical(p$labels$colour, "FeatureID")
    expect_length(unique(ld$PANEL), 1L)
    expect_length(unique(ld$group), 4L)
    # Without time.lag, x is the lag in steps; without n.pairs, no band.
    one <- acfs[acfs$FeatureID == "f1" & acfs$subject == "a",
        c("lag", "acf")]
    p <- plotLagProfile(one, feature.col = NULL)
    expect_length(p$layers, 2L)
    expect_equal(layer_data(p, 2L)$x, 1:3)
    expect_error(plotLagProfile(one), "'feature.col' must be NULL")
    expect_error(plotLagProfile(acfs, group = "id"), "'group' must be NULL")
    toy <- SummarizedExperiment(assays = list(counts = matrix(1, 1L, 1L)))
    S4Vectors::metadata(toy)$autocorrelation <- S4Vectors::DataFrame(acfs)
    expect_equal(
        plotLagProfile(toy, "autocorrelation", group = "subject")$data,
        plotLagProfile(acfs, group = "subject")$data)
})

test_that("plotLagProfile draws the fitted variogram of each level", {
    model <- data.frame(sex = rep(c("female", "male"), each = 3),
        type = "parameter", form = "exponential",
        parameter = c("nugget", "partial.sill", "timescale"),
        estimate = c(0.01, 0.2, 1, 0, 0.3, 2), min.lag = 0.2, max.lag = 8)
    p <- plotLagProfile(tse, name = "by_sex", model = model)
    ld <- layer_data(p, 1L)
    expect_length(unique(ld$group), 2L)
    expect_equal(range(ld$x), c(0.2, 8))
    f <- ld[ld$colour == "#0072B2", ]
    expect_equal(f$y, 0.01 + 0.2 * (1 - exp(-f$x / 1)))
    # Curves and estimates share the colour of their level.
    expect_setequal(ld$colour, layer_data(p, 2L)$colour)
    # Stored by addVariogramModel()
    metadata(tse)$model <- model
    expect_equal(layer_data(plotLagProfile(tse, name = "by_sex",
        model = "model"), 1L), ld)
})

test_that("plotLagProfile evaluates each form of getVariogramModel", {
    h <- c(0.5, 2, 6)
    est <- c(nugget = 0.1, partial.sill = 0.2, timescale = 2,
        partial.sill.1 = 0.1, partial.sill.2 = 0.3, timescale.1 = 1,
        timescale.2 = 5, scale = 0.05, exponent = 0.5, slope = 0.02,
        seasonal = 0.04)
    f <- tempoPlot:::.variogram_form
    expect_equal(f(h, "exponential", est), 0.1 + 0.2 * (1 - exp(-h / 2)))
    expect_equal(f(h, "two.timescale", est), 0.1 +
        0.1 * (1 - exp(-h)) + 0.3 * (1 - exp(-h / 5)))
    expect_equal(f(h, "power", est), 0.1 + 0.05 * sqrt(h))
    expect_equal(f(h, "linear", est), 0.1 + 0.02 * h)
    expect_equal(f(h, "linear.seasonal", est, period = 4),
        0.1 + 0.02 * h + 0.04 * (1 - cos(2 * pi * h / 4)))
})

test_that("plotLagProfile checks the model", {
    model <- data.frame(type = "parameter", form = "linear",
        parameter = c("nugget", "slope"), estimate = c(0.1, 0.02))
    # Without min.lag, the curve spans the lags of the variogram.
    ld <- layer_data(plotLagProfile(v, model = model), 1L)
    expect_equal(range(ld$x), range(v$lag[!is.na(v$estimate)]))
    expect_error(plotLagProfile(vs, model = model), "same split.by")
    expect_error(plotLagProfile(vs, type = "ratio", model = model),
        "dissimilarity")
    expect_error(plotLagProfile(v, model = model[, -2]), "getVariogramModel")
    expect_error(plotLagProfile(tse, model = "missing"), "addVariogramModel")
    model$form <- "linear.seasonal"
    expect_error(plotLagProfile(v, model = model), "period")
})

test_that("plotLagProfile draws a getVariogramModel fit", {
    skip_if_not("getVariogramModel" %in% getNamespaceExports("miaTime"))
    set.seed(1)
    fit <- miaTime::getVariogramModel(tse, time.col = "time",
        group = "subject", assay.type = "hellinger", method = "euclidean",
        squared = TRUE, form = "exponential")
    ld <- layer_data(plotLagProfile(tse, model = fit), 1L)
    est <- stats::setNames(fit$estimate, fit$parameter)
    expect_equal(ld$y, est[["nugget"]] + est[["partial.sill"]] *
        (1 - exp(-ld$x / est[["timescale"]])))
})
