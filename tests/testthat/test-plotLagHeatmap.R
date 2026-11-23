library(SummarizedExperiment)
library(ggplot2)

# Feature a equals the target one sample later; b is the target itself.
set.seed(1)
v <- rnorm(20)
toy <- SummarizedExperiment(
    assays = list(counts = rbind(a = c(v[-1], 0), b = v)),
    colData = DataFrame(t = 1:20, y = v))

test_that("plotLagHeatmap checks its arguments", {
    expect_error(plotLagHeatmap(toy, "y", "t", max.lag = -1), "'max.lag'")
    expect_error(plotLagHeatmap(toy, "none", "t"))
})

test_that("plotLagHeatmap finds the lag of a shifted feature", {
    df <- plotLagHeatmap(toy, "y", "t", max.lag = 2)$data
    r <- function(f, k) df$.r[df$.feature == f & df$.lag == k]
    expect_equal(r("b", 0), 1)
    # a at sample i is the target at sample i + 1: a leads by one sample.
    expect_equal(r("a", 1), 1)
    expect_lt(abs(r("a", -1)), 1)
    expect_setequal(df$.lag, -2:2)
})

test_that("plotLagHeatmap pairs samples within subjects only", {
    toy$id <- rep(c("s1", "s2"), each = 10)
    # The last sample of s1 has no later sample in s1; pairing it with the
    # first sample of s2 would add an outlier.
    assay(toy)["a", 10] <- -100
    df <- plotLagHeatmap(toy, "y", "t", group = "id", max.lag = 1)$data
    expect_equal(df$.r[df$.feature == "a" & df$.lag == 1], 1)
})

test_that("plotLagHeatmap correlates within subjects", {
    # Within each subject a falls as y rises, but s2 is higher in both.
    se <- SummarizedExperiment(
        assays = list(counts = rbind(a = c(-v[1:10], 10 - v[11:20]))),
        colData = DataFrame(t = rep(1:10, 2), y = c(v[1:10], 10 + v[11:20]),
            id = rep(c("s1", "s2"), each = 10)))
    df <- plotLagHeatmap(se, "y", "t", group = "id", max.lag = 0)$data
    expect_equal(df$.r, -1)
})

test_that("plotLagHeatmap reads a feature target from assay.type", {
    se <- toy
    assay(se, "alt") <- rbind(a = rev(v), b = rev(v))
    df <- plotLagHeatmap(se, "b", "t", assay.type = "alt", max.lag = 0)$data
    expect_equal(df$.r, c(1, 1))
})

test_that("plotLagHeatmap takes the target from another experiment", {
    toy$id <- rep(c("s1", "s2"), each = 10)
    colnames(toy) <- paste0("p", 1:20)
    cd <- colData(toy)
    tgt <- SummarizedExperiment(matrix(cd$y, 1L,
        dimnames = list("y", colnames(toy))))
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        list(features = toy, target = tgt[, 20:1]), colData = cd)
    p <- plotLagHeatmap(mae, target = "y", time.col = "t", group = "id",
        experiment = "features", target.experiment = "target",
        max.lag = 1)
    q <- plotLagHeatmap(toy, "y", "t", group = "id", max.lag = 1)
    expect_equal(p$data$.r, q$data$.r)
    expect_error(plotLagHeatmap(mae, "y", "t", experiment = "none",
        target.experiment = "target"), "'experiment'")
    expect_error(plotLagHeatmap(mae, "y", "none", experiment = "features",
        target.experiment = "target"), "'time.col'")
})

test_that("plotLagHeatmap draws negative blue and positive red", {
    se <- SummarizedExperiment(assays = list(counts = rbind(b = v, c = -v)),
        colData = colData(toy))
    p <- plotLagHeatmap(se, "y", "t", max.lag = 0)
    fill <- setNames(layer_data(p)$fill, as.character(p$data$.feature))
    expect_identical(unname(fill[c("b", "c")]), c("#B2182B", "#2166AC"))
    expect_identical(p$theme$legend.position, "bottom")
})

test_that("plotLagHeatmap gives NA without a warning for a constant feature", {
    se <- SummarizedExperiment(assays = list(counts = rbind(b = v,
        flat = 1)), colData = colData(toy))
    expect_no_warning(p <- plotLagHeatmap(se, "y", "t", max.lag = 0))
    expect_true(is.na(p$data$.r[p$data$.feature == "flat"]))
    expect_identical(layer_data(p)$fill[p$data$.feature == "flat"],
        "grey50")
})

test_that("plotLagHeatmap reads counts by default", {
    se <- toy
    assays(se) <- list(other = assay(se))
    expect_error(plotLagHeatmap(se, "y", "t"), "'assay.type'")
})

test_that("plotLagHeatmap counts lags in visits for a visit factor", {
    # Visit 3 is missing: at lag 1, visit 2 has no partner and visit 4 is
    # not paired with visit 2.
    se <- toy[, -3]
    se$v <- factor(se$t, levels = 1:20)
    df <- plotLagHeatmap(se, "y", "v", max.lag = 1)$data
    expect_equal(df$.r[df$.feature == "a" & df$.lag == 1], 1)
    expect_match(plotLagHeatmap(se, "y", "v", max.lag = 1)$labels$x,
        "lag in visits")
    # By samples, visit 2 pairs with visit 4 and breaks the shift.
    df <- plotLagHeatmap(se, "y", "t", max.lag = 1)$data
    expect_lt(df$.r[df$.feature == "a" & df$.lag == 1], 1)
})
