library(SummarizedExperiment)
library(ggplot2)

# Partner P has a and b. Subject f has b and c at its event (time 0), then
# a, b, c, d and later only e; subject g has a, then a and e.
counts <- matrix(c(
    5, 0, 1, 0, 1, 1,
    5, 1, 2, 0, 0, 0,
    0, 1, 3, 0, 0, 0,
    0, 0, 4, 0, 0, 0,
    0, 0, 0, 0.5, 0, 1), 5L, byrow = TRUE,
    dimnames = list(letters[1:5], c("P", "f0", "f1", "f2", "g0", "g1")))
toy <- SummarizedExperiment(assays = list(counts = counts),
    colData = DataFrame(id = c("P", "f", "f", "f", "g", "g"),
        t = c(NA, 0, 1, 2, 0, 1), ev = c(NA, 0, 0, 0, 0, 0),
        partner = c(NA, "P", "P", "P", "P", "P"),
        arm = c(NA, "x", "x", "x", "x", "x")))
shares <- function(df) {
    split(df$.value, paste(df$id, df$t))
}

test_that("plotOrigin checks its arguments", {
    expect_error(plotOrigin(toy, "t", "id", "none"), "'partner.col'")
    expect_error(plotOrigin(toy, "t", "id", "partner", value = "count"))
    expect_error(plotOrigin(toy, "t", "id", "partner", assay.type = "x"),
        "'assay.type'")
    expect_error(plotOrigin(toy, "t", "id", "partner", align.to = "none"),
        "'align.to'")
    visits <- toy
    visits$t <- factor(c(NA, 1, 2, 3, 1, 2))
    visits$ev <- factor(c(NA, 1, 1, 1, 1, 1))
    expect_error(plotOrigin(visits, "t", "id", "partner", align.to = "ev"),
        "numeric, Date or POSIXct")
})

test_that("plotOrigin splits features by partner and reference", {
    p <- plotOrigin(toy, "t", "id", "partner", align.to = "ev")
    df <- p$data
    expect_identical(levels(df$.feature), c("partner", "both", "own", "new"))
    # Reference samples (time 0) are not drawn; f1 has one feature of each
    # origin, f2 only the new e, g1 a from both and the new e.
    expect_equal(shares(df), list(`f 1` = rep(0.25, 4),
        `f 2` = c(0, 0, 0, 1), `g 1` = c(0, 0.5, 0, 0.5)))
    df <- plotOrigin(toy, "t", "id", "partner", value = "abundance",
        align.to = "ev")$data
    expect_equal(shares(df)$`f 1`, c(0.1, 0.2, 0.3, 0.4))
    expect_identical(p$labels$x, "t since ev")
    expect_identical(p$labels$fill, "origin")
    # "new" is grey.
    ld <- layer_data(p, 1L)
    expect_identical(unique(ld$fill[ld$group == 4L]), "grey50")
    expect_identical(unique(ld$fill[ld$group == 1L]), "#0072B2")
})

test_that("plotOrigin without a reference splits partner and new", {
    df <- plotOrigin(toy, "t", "id", "partner")$data
    expect_identical(levels(df$.feature), c("partner", "new"))
    expect_equal(shares(df), list(`f 0` = c(0.5, 0.5), `f 1` = c(0.5, 0.5),
        `f 2` = c(0, 1), `g 0` = c(1, 0), `g 1` = c(0.5, 0.5)))
})

test_that("plotOrigin applies detection and drops empty samples", {
    # Above 1: P has a and b, f0 nothing, f1 b, c and d, f2 nothing.
    df <- plotOrigin(toy, "t", "id", "partner", detection = 1,
        align.to = "ev")$data
    expect_equal(shares(df), list(`f 1` = c(1, 0, 0, 2) / 3))
})

test_that("plotOrigin averages subjects by group", {
    p <- plotOrigin(toy, "t", "id", "partner", align.to = "ev",
        level = "group", group.by = "arm")
    df <- p$data[order(p$data$t, p$data$.feature), ]
    expect_equal(df$.value, c(0.125, 0.375, 0.125, 0.375, 0, 0, 0, 1))
    p <- plotOrigin(toy, "t", "id", "partner", align.to = "ev",
        layout = "area")
    ld <- layer_data(p, 1L)
    expect_equal(unique(as.vector(round(tapply(ld$ymax, ld$x, max), 10))),
        1)
})
