library(SummarizedExperiment)
library(ggplot2)

# A series of n days whose periodogram falls as frequency^-beta: cosines
# at the Fourier frequencies with amplitude f^(-beta / 2), random phases.
# Long series keep the slope estimates within 0.3 of -beta.
powerlaw <- function(beta, n = 400L) {
    f <- seq_len(n / 2L) / n
    phase <- runif(length(f), 0, 2 * pi)
    colSums(f^(-beta / 2) * cos(2 * pi * outer(f, seq_len(n)) + phase))
}
# Slopes 0, -1, -2 and -3 in subjects a and b, and a constant feature.
set.seed(1)
mat <- rbind(
    white = c(powerlaw(0), powerlaw(0)),
    pink = c(powerlaw(1), powerlaw(1)),
    brown = c(powerlaw(2), powerlaw(2)),
    black = c(powerlaw(3), powerlaw(3)),
    constant = 1)
toy <- SummarizedExperiment(assays = list(counts = mat),
    colData = DataFrame(day = rep(1:400, 2L), id = rep(c("a", "b"),
        each = 400L), visit = factor(rep(1:400, 2L))))

test_that("plotNoise classifies series of known spectral slope", {
    p <- plotNoise(toy, "day", "id")
    df <- p$data
    # The constant feature has no colour.
    expect_setequal(df$.feature, c("white", "pink", "brown", "black"))
    beta <- c(white = 0, pink = 1, brown = 2, black = 3)
    expect_true(all(abs(df$slope + beta[df$.feature]) < 0.3))
    expect_identical(as.character(df$noise), df$.feature)
    # One feature of four in each colour per subject.
    ld <- layer_data(p)
    expect_equal(ld$ymax - ld$ymin, rep(0.25, 8L))
    expect_setequal(ld$fill, c("white", "pink", "sienna4", "black"))
    expect_identical(levels(df$noise), c("white", "pink", "brown", "black"))
    expect_identical(p$labels$x, "id")
})

test_that("plotNoise tells white noise from a random walk", {
    set.seed(2)
    walk <- SummarizedExperiment(colData = colData(toy),
        assays = list(counts = rbind(noise = rnorm(800),
            walk = c(cumsum(rnorm(400)), cumsum(rnorm(400))))))
    df <- plotNoise(walk, "day", "id")$data
    white <- df$.feature == "noise"
    expect_true(all(df$noise[white] == "white"))
    expect_true(all(df$slope[!white] < -1))
})

test_that("plotNoise checks its input", {
    expect_error(plotNoise(toy, "visit", "id"), "numeric, Date or POSIXct")
    expect_error(plotNoise(toy["constant", ], "day", "id"),
        "No feature varies")
    expect_error(plotNoise(toy, "day", "id", assay.type = "clr"),
        "'assay.type' must be one of assayNames")
})

test_that(".noise_colour uses the thresholds of Faust et al.", {
    expect_identical(as.character(.noise_colour(
        c(-3, -2.25, -2, -1.75, -1, -0.5, 0.3))),
        c("black", "brown", "brown", "pink", "pink", "white", "white"))
})
