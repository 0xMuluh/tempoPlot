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
