library(SummarizedExperiment)
library(ggplot2)

# Subjects a and b at x = 1, 2, 3. Centred on the subject means, x is
# (-1, 0, 1) in both, y is (-1, 1, 0) in a and (-1, -1, 2) in b, so the
# common slope is 4 / 4 = 1, r_rm = 4 / sqrt(4 * 8) = 1 / sqrt(2) and
# df = 6 samples - 2 subjects - 1 = 3. Sample s7 has no subject.
toy <- SummarizedExperiment(
    assays = list(counts = matrix(c(1, 3, 2, 10, 10, 13, 5), 1L,
        dimnames = list("f", paste0("s", 1:7)))),
    colData = DataFrame(id = c(1, 1, 1, 2, 2, 2, NA),
        u = c(1, 2, 3, 1, 2, 3, 4), flat = c(1, 1, 1, 2, 2, 2, 3)))

test_that("plotRmcorr checks its arguments", {
    expect_error(plotRmcorr(toy, "u", "f", group = "none"), "'group'")
    expect_error(plotRmcorr(toy, "none", "f", group = "id"), "'x.var'")
    expect_error(plotRmcorr(toy, "u", "f", group = "id",
        assay.type = "none"), "'assay.type'")
    expect_error(plotRmcorr(toy, "flat", "f", group = "id"),
        "vary within subjects")
    expect_error(plotRmcorr(toy[, c(1, 4)], "u", "f", group = "id"),
        "vary within subjects")
})

test_that("plotRmcorr computes the repeated-measures correlation", {
    p <- plotRmcorr(toy, "u", "f", group = "id")
    r <- 1 / sqrt(2)
    t <- r * sqrt(3 / (1 - r^2))
    expect_identical(get_labs(p)$subtitle, sprintf(
        "r_rm = %.2f, df = %d, p = %.2g", r, 3L, 2 * pt(-t, 3)))
    # The same from the analysis of covariance of Bakdash and Marusich.
    df <- data.frame(id = factor(toy$id), x = toy$u,
        y = assay(toy)[1L, ])[1:6, ]
    fit <- lm(y ~ id + x, data = df)
    a <- anova(fit)
    rmc <- .rmcorr(df$x, df$y, df$id)
    expect_equal(rmc$r, sign(coef(fit)[["x"]]) *
        sqrt(a["x", "Sum Sq"] / sum(a[c("x", "Residuals"), "Sum Sq"])))
    expect_equal(rmc$df, a["Residuals", "Df"])
    expect_equal(rmc$p, a["x", "Pr(>F)"])
    expect_equal(rmc$slope, coef(fit)[["x"]])
})

test_that("plotRmcorr draws parallel lines through the subject means", {
    p <- plotRmcorr(toy, "u", "f", group = "id")
    # The sample without a subject is dropped; numeric ids are discrete.
    expect_identical(nrow(p$data), 6L)
    expect_s3_class(p$data$.subject, "factor")
    # a: y = 2 + (x - 2); b: y = 11 + (x - 2).
    expect_equal(p$data$.fit, c(1, 2, 3, 10, 11, 12))
    lines <- layer_data(p, 1L)
    expect_identical(length(unique(lines$colour)), 2L)
    expect_identical(get_labs(p)$x, "u")
    expect_identical(get_labs(p)$y, "f")
})

test_that("plotRmcorr drops the legend beyond six subjects", {
    expect_identical(plotRmcorr(toy, "u", "f", group = "id")$scales$
        get_scales("colour")$guide, "legend")
    many <- toy[, rep(1:6, 4)]
    colnames(many) <- NULL
    many$id <- rep(1:8, each = 3)
    many$u <- rep(1:3, 8)
    p <- plotRmcorr(many, "u", "f", group = "id")
    expect_identical(p$scales$get_scales("colour")$guide, "none")
})
