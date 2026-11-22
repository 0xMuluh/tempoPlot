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
    expect_error(plotSpectrum(toy, "v", "f"), "numeric or Date")
    expect_error(plotSpectrum(toy[, 1:3], "v", "t"), "four or more")
})

test_that("plotSpectrum finds the period of a sine", {
    df <- plotSpectrum(toy, "v", "t")$data
    expect_equal(df$frequency[which.max(df$power)], 0.1, tolerance = 0.05)
    p <- plotSpectrum(toy, "v", "t", group = "id")
    expect_setequal(p$data$.subject, c("a", "b"))
    expect_identical(p$labels$colour, "id")
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
