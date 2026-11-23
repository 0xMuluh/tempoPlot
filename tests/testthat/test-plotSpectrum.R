library(SummarizedExperiment)
library(ggplot2)

# A sine with period 10 sampled unevenly.
set.seed(1)
t <- sort(runif(120, 0, 100))
toy <- SummarizedExperiment(
    assays = list(counts = matrix(1, 1L, 120L)),
    colData = DataFrame(t = t, v = sin(2 * pi * t / 10),
        id = rep(c("a", "b"), 60), f = factor(seq_len(120))))

test_that("plotSpectrum checks its arguments", {
    expect_error(plotSpectrum(toy, "v", "f"), "numeric, Date or POSIXct")
    expect_error(plotSpectrum(toy[, 1:3], "v", "t"), "four or more")
    se <- toy
    rownames(se) <- "taxon"
    expect_error(plotSpectrum(se, "taxon", "t", assay.type = "clr"),
        "'assay.type' must be one of assayNames")
})

test_that("plotSpectrum finds the period of a sine", {
    df <- plotSpectrum(toy, "v", "t")$data
    expect_equal(df$frequency[which.max(df$power)], 0.1, tolerance = 0.05)
    p <- plotSpectrum(toy, "v", "t", group = "id")
    expect_setequal(p$data$.subject, c("a", "b"))
    expect_identical(p$labels$colour, "id")
    expect_setequal(layer_data(p)$colour, c("#0072B2", "#E69F00"))
    for (s in c("a", "b")) {
        d <- p$data[p$data$.subject == s, ]
        expect_equal(d$frequency[which.max(d$power)], 0.1,
            tolerance = 0.05)
    }
})

test_that("plotSpectrum gives Date frequencies per day", {
    days <- toy
    days$t <- as.Date("2020-01-01") + round(days$t)
    days <- days[, !duplicated(days$t)]
    df <- plotSpectrum(days, "v", "t")$data
    expect_equal(df$frequency[which.max(df$power)], 0.1, tolerance = 0.05)
})

test_that("plotSpectrum gives POSIXct frequencies per hour", {
    hours <- toy
    hours$t <- as.POSIXct("2020-01-01", tz = "UTC") + hours$t * 3600
    df <- plotSpectrum(hours, "v", "t")$data
    expect_equal(df$frequency[which.max(df$power)], 0.1, tolerance = 0.05)
})

test_that("plotSpectrum reads features from counts by default", {
    se <- toy
    assays(se) <- list(other = matrix(0, 1L, 120L),
        counts = matrix(se$v, 1L))
    rownames(se) <- "g"
    df <- plotSpectrum(se, "g", "t")$data
    expect_equal(df$frequency[which.max(df$power)], 0.1, tolerance = 0.05)
})

test_that("plotSpectrum fits a straight line on log-log axes", {
    p <- plotSpectrum(toy, "v", "t", group = "id", fit = "power")
    expect_s3_class(p$scales$get_scales("x")$trans, "transform")
    expect_s3_class(p$scales$get_scales("y")$trans, "transform")
    ld <- layer_data(p, 2L)
    expect_true(all(ld$linetype == "dashed"))
    for (s in c("a", "b")) {
        d <- p$data[p$data$.subject == s, ]
        b <- coef(lm(log10(power) ~ log10(frequency), d))[[2L]]
        fitted <- ld[ld$group == match(s, c("a", "b")), ]
        expect_equal(coef(lm(y ~ x, fitted))[[2L]], b)
        expect_match(p$labels$subtitle, sprintf("%s: slope %.2f", s, b),
            fixed = TRUE)
    }
    p <- plotSpectrum(toy, "v", "t", fit = "power")
    expect_match(p$labels$subtitle, "^slope -?[0-9.]+$")
    expect_true(all(layer_data(p, 2L)$colour == "black"))
})
