library(SummarizedExperiment)
library(ggplot2)

# s1: (4, 0, 2) -> (2, 3, 2); s2 only at t = 1.
toy <- SummarizedExperiment(
    assays = list(counts = matrix(c(
        4, 2, 1,
        0, 3, 1,
        2, 2, 1), 3L, byrow = TRUE,
        dimnames = list(c("a", "b", "c"), paste0("c", 1:3)))),
    colData = DataFrame(id = c("s1", "s1", "s2"), t = c(1, 2, 1),
        sex = c("f", "m", "f")))

test_that("plotTBI checks its arguments", {
    expect_error(suppressWarnings(plotTBI(toy, "t", NULL, c(1, 2))),
        "'group'")
    expect_error(suppressWarnings(plotTBI(toy, "t", "id", c(1, 9))),
        "'time.point'")
    expect_error(suppressWarnings(plotTBI(toy, "t", "id", c(1, 1))),
        "two different values")
    expect_error(suppressWarnings(plotTBI(toy, "t", "id", c(1, 2),
        method = "jaccard")))
})

test_that("plotTBI splits Bray-Curtis into losses and gains", {
    p <- suppressWarnings(plotTBI(toy, "t", "id", c(1, 2), method = "bray",
        colour.by = "sex"))
    df <- p$data
    expect_identical(df$.subject, "s1")
    # A = 2 + 0 + 2 = 4, B = 2, C = 3; denominator 13.
    expect_equal(df$losses, 2 / 13)
    expect_equal(df$gains, 3 / 13)
    expect_equal(df$losses + df$gains,
        sum(abs(c(4, 0, 2) - c(2, 3, 2))) / sum(c(4, 0, 2, 2, 3, 2)))
    expect_identical(df$sex, "f")
    expect_identical(layer_data(p, 2L)$colour, "#0072B2")
})

test_that("plotTBI uses presence and absence by default", {
    df <- suppressWarnings(plotTBI(toy, "t", "id", c(1, 2)))$data
    # Presence (1, 0, 1) -> (1, 1, 1): A = 2, B = 0, C = 1.
    expect_equal(df$losses, 0)
    expect_equal(df$gains, 1 / 5)
})

test_that("plotTBI bray losses equal gains on relative abundances", {
    # B - C is the difference of the sample totals.
    rel <- toy
    assay(rel) <- t(t(assay(rel)) / colSums(assay(rel)))
    df <- suppressWarnings(plotTBI(rel, "t", "id", c(1, 2),
        method = "bray"))$data
    expect_equal(df$losses, df$gains)
})

test_that("plotTBI reads counts by default", {
    se <- toy
    assays(se) <- list(other = assay(se) * 0, counts = assay(se))
    df <- suppressWarnings(plotTBI(se, "t", "id", c(1, 2)))$data
    expect_equal(df$gains, 1 / 5)
    assays(se) <- list(other = assay(se))
    expect_error(suppressWarnings(plotTBI(se, "t", "id", c(1, 2))),
        "'assay.type'")
})
