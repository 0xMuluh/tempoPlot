library(SummarizedExperiment)

data("hitchip1006", package = "miaTime")

# Three subjects of hitchip1006 with at least three samples each. In the
# MultiAssayExperiment, colData lives in the MAE: experiment "fwd" keeps
# only its own column "extra", and "rev" has the samples in reverse order.
ids <- names(which(table(hitchip1006$subject) >= 3L))[1:3]
tse <- hitchip1006[, hitchip1006$subject %in% ids]
tse_rev <- tse[, rev(seq_len(ncol(tse)))]
fwd <- tse
colData(fwd) <- DataFrame(extra = seq_len(ncol(tse)),
    time = -1, row.names = colnames(tse))
rev_se <- tse_rev
colData(rev_se) <- colData(rev_se)[, 0L]
mae <- MultiAssayExperiment::MultiAssayExperiment(
    list(fwd = fwd, rev = rev_se), colData = colData(tse))

test_that(".get_experiment joins colData by primary name", {
    se <- .get_experiment(mae, "rev")
    expect_identical(colnames(se), colnames(tse_rev))
    expect_equal(colData(se), colData(tse_rev))
    # Own columns are kept; colData(x) wins on a shared name.
    se <- .get_experiment(mae, 1L)
    expect_identical(se$extra, seq_len(ncol(tse)))
    expect_identical(se$time, tse$time)
})

test_that("MultiAssayExperiment methods plot the chosen experiment", {
    expect_equal(
        plotTrend(mae, "diversity", "time", "subject",
            experiment = "rev")$data,
        plotTrend(tse_rev, "diversity", "time", "subject")$data)
    expect_equal(
        plotTiles(mae, "time", "subject", assay.type = "counts", top = 3,
            experiment = 2L)$data,
        plotTiles(tse_rev, "time", "subject", assay.type = "counts",
            top = 3)$data)
    # experiment defaults to the first one.
    feature <- rownames(tse)[1L]
    expect_equal(
        plotLagHeatmap(mae, feature, "time", "subject",
            assay.type = "counts", top = 3, max.lag = 1)$data,
        plotLagHeatmap(tse, feature, "time", "subject",
            assay.type = "counts", top = 3, max.lag = 1)$data)
})

test_that("plotLagHeatmap takes the target from target.experiment", {
    # The target experiment holds a subset of the samples, reversed.
    keep <- rev(colnames(tse)[-(1:2)])
    tgt <- SummarizedExperiment(
        matrix(tse$diversity[match(keep, colnames(tse))], 1L,
            dimnames = list("div", keep)))
    m <- MultiAssayExperiment::MultiAssayExperiment(
        list(fwd = fwd, tgt = tgt), colData = colData(tse))
    se <- tse[, -(1:2)]
    se$div <- se$diversity
    expect_equal(
        plotLagHeatmap(m, "div", "time", "subject", assay.type = "counts",
            top = 3, max.lag = 1, target.experiment = "tgt")$data,
        plotLagHeatmap(se, "div", "time", "subject", assay.type = "counts",
            top = 3, max.lag = 1)$data)
    expect_error(plotLagHeatmap(m, "div", "time", "subject",
        target.experiment = "none"), "'target.experiment'")
})

test_that("MultiAssayExperiment methods check experiment", {
    for (e in list("none", 3L, 0L, -1L, c(1L, 2L), NA_character_)) {
        expect_error(plotTrend(mae, "diversity", "time", "subject",
            experiment = e), "'experiment' must be the name or index")
    }
    m <- MultiAssayExperiment::MultiAssayExperiment(
        list(mat = assay(tse)), colData = colData(tse))
    expect_error(plotTrend(m, "diversity", "time", "subject",
        experiment = "mat"), "'experiment' must be a SummarizedExperiment")
})
