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
    expect_error(plotSpectrum(toy, "v", "t", method = "wavelet"))
    expect_error(plotSpectrum(toy[, 1:3], "v", "t"), "four or more")
})

test_that("plotSpectrum finds the period of a sine", {
    df <- plotSpectrum(toy, "v", "t")$data
    expect_equal(df$frequency[which.max(df$power)], 0.1, tolerance = 0.05)
    p <- plotSpectrum(toy, "v", "t", group = "id")
    expect_setequal(p$data$.subject, c("a", "b"))
    for (s in c("a", "b")) {
        d <- p$data[p$data$.subject == s, ]
        expect_equal(d$frequency[which.max(d$power)], 0.1,
            tolerance = 0.05)
    }
})

test_that("plotSpectrum uses the median interval for fft", {
    even <- toy
    even$t <- seq(0, by = 0.5, length.out = 120)
    even$v <- sin(2 * pi * even$t / 10)
    df <- plotSpectrum(even, "v", "t", method = "fft")$data
    expect_equal(df$frequency[which.max(df$power)], 0.1, tolerance = 0.05)
})
