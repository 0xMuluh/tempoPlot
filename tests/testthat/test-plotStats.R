library(SummarizedExperiment)
library(ggplot2)

data("hitchip1006", package = "miaTime")
tse <- miaTime::addBimodality(hitchip1006, assay.type = "counts")
bim <- setNames(rowData(tse)$bimodality, rownames(tse))

# Per-subject mean diversity with its range over visits.
subj <- tse[, tse$subject %in% names(which(table(tse$subject) >= 3))]
subj$mean <- ave(subj$diversity, subj$subject, FUN = mean)
subj$min <- ave(subj$diversity, subj$subject, FUN = min)
subj$max <- ave(subj$diversity, subj$subject, FUN = max)

test_that("plotStats checks its arguments", {
    expect_error(plotStats(tse, stat = "missing"), "'stat'")
    expect_error(plotStats(tse, stat = NULL), "'stat'")
    expect_error(plotStats(tse, stat = "bimodality", lower = "bimodality"),
        "together")
    expect_error(plotStats(tse, stat = "Phylum"), "numeric")
    expect_error(plotStats(tse, stat = "bimodality", top = 0), "'top'")
    expect_error(plotStats(tse, stat = "bimodality", colour.by = "missing"),
        "'colour.by'")
    # Subject statistics come from colData, which has no bimodality.
    expect_error(plotStats(tse, stat = "bimodality", group = "subject"),
        "colData")
    expect_error(plotStats(subj, stat = "mean", group = "subject",
        colour.by = "bmi_group"), "constant")
})

test_that("plotStats ranks features by their statistic, largest at the top", {
    ld <- layer_data(plotStats(tse, stat = "bimodality", top = 20), 1L)
    top20 <- sort(bim[order(-abs(bim))][1:20])
    expect_identical(nrow(ld), 20L)
    expect_equal(ld$x[order(ld$y)], unname(top20))
    expect_equal(nrow(layer_data(plotStats(tse, stat = "bimodality"), 1L)),
        sum(!is.na(bim)))
})

test_that("plotStats keeps the largest absolute statistics", {
    data("crohn_survival", package = "miaTime")
    set.seed(1)
    cs <- suppressWarnings(miaTime::addSurvival(crohn_survival,
        assay.type = "counts", time.col = "event_time", event.col = "event"))
    coef <- rowData(cs)$survival_coef
    expect_true(any(coef < 0))
    ld <- layer_data(plotStats(cs, stat = "survival_coef", top = 5), 1L)
    expect_setequal(ld$x, coef[order(-abs(coef))][1:5])
})

test_that("plotStats colours points and drops missing statistics", {
    p <- plotStats(tse, stat = "bimodality", top = 20, colour.by = "Phylum")
    phyla <- rowData(tse)$Phylum[match(names(sort(bim[order(-abs(bim))][1:20])),
        rownames(tse))]
    expect_length(unique(layer_data(p, 1L)$colour), length(unique(phyla)))

    na <- tse
    rowData(na)$bimodality[1:5] <- NA
    expect_identical(nrow(layer_data(plotStats(na, stat = "bimodality"), 1L)),
        sum(!is.na(rowData(na)$bimodality)))
})

test_that("plotStats draws one row per subject with its interval", {
    p <- plotStats(subj, stat = "mean", lower = "min", upper = "max",
        group = "subject", colour.by = "sex", size = 3)
    ranges <- layer_data(p, 1L)
    points <- layer_data(p, 2L)
    expect_identical(nrow(points), length(unique(subj$subject)))
    expect_true(all(points$size == 3))
    expect_equal(sort(ranges$xmin), sort(tapply(subj$diversity,
        subj$subject, min)[unique(as.character(subj$subject))]),
        ignore_attr = TRUE)
    expect_equal(sort(ranges$xmax), sort(tapply(subj$diversity,
        subj$subject, max)[unique(as.character(subj$subject))]),
        ignore_attr = TRUE)
    expect_identical(get_labs(p)$y, "subject")
})
