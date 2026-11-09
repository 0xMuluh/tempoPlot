library(SummarizedExperiment)
library(SingleCellExperiment)
library(ggplot2)

data("hitchip1006", package = "miaTime")
# Subjects sampled five times, without the profiles that are stored twice.
x <- hitchip1006[, !duplicated(t(assay(hitchip1006, "counts")))]
x <- x[, x$subject %in% names(which(table(x$subject) == 5))]
x$subject <- droplevels(x$subject)
pca <- prcomp(t(log1p(assay(x, "counts"))))$x[, 1:3]
reducedDim(x, "PCA") <- pca
n_subjects <- nlevels(x$subject)
ord <- order(x$subject, x$time)

test_that("plotPath checks its arguments", {
    expect_error(plotPath(x, "time", "subject"), "either")
    expect_error(plotPath(x, "time", "subject", dimred = "PCA",
        x.var = "diversity", y.var = "time"), "either")
    expect_error(plotPath(x, "time", "subject", x.var = "diversity"),
        "either")
    expect_error(plotPath(x, "time", "subject", dimred = "MDS"),
        "reducedDimNames")
    expect_error(plotPath(x, "time", "subject", dimred = "PCA",
        ncomponents = 1:3), "two")
    expect_error(plotPath(x, "time", "subject", dimred = "PCA",
        ncomponents = c(1, 9)), "two")
    expect_error(plotPath(x, "time", "subject", x.var = "missing",
        y.var = "diversity"), "'x.var' was not found")
    expect_error(plotPath(x, "time", "subject", dimred = "PCA",
        arrows = "yes"), "'arrows'")
})

test_that("plotPath joins each subject's samples in time order", {
    p <- plotPath(x, "time", "subject", dimred = "PCA")
    ld <- layer_data(p, 1L)
    expect_length(unique(ld$group), n_subjects)
    expect_equal(ld$x, unname(pca[ord, 1L]))
    expect_equal(ld$y, unname(pca[ord, 2L]))
    expect_identical(get_labs(p)$x, "PCA:1")
    expect_false(is.null(p$layers[[1L]]$geom_params$arrow))

    p <- plotPath(x, "time", "subject", dimred = "PCA",
        ncomponents = c("PC3", "PC2"), arrows = FALSE)
    ld <- layer_data(p, 1L)
    expect_equal(ld$x, unname(pca[ord, 3L]))
    expect_equal(ld$y, unname(pca[ord, 2L]))
    expect_null(p$layers[[1L]]$geom_params$arrow)
})

test_that("plotPath plots two variables", {
    genus <- "Bacteroides vulgatus et rel."
    p <- plotPath(x, "time", "subject", x.var = genus, y.var = "diversity",
        assay.type = "counts", linewidth = 2)
    ld <- layer_data(p, 1L)
    expect_equal(ld$x, unname(assay(x, "counts")[genus, ord]))
    expect_equal(ld$y, x$diversity[ord])
    expect_true(all(ld$linewidth == 2))
    expect_identical(get_labs(p)$y, "diversity")
})

test_that("plotPath starts paths at zero with origin = 'baseline'", {
    ld <- layer_data(plotPath(x, "time", "subject", dimred = "PCA",
        origin = "baseline"), 1L)
    first <- !duplicated(ld$group)
    expect_equal(ld$x[first], rep(0, n_subjects))
    expect_equal(ld$y[first], rep(0, n_subjects))
    # Displacements from the first sample are kept.
    expect_equal(diff(ld$x)[!first[-1L]],
        diff(unname(pca[ord, 1L]))[!first[-1L]])
})

test_that("plotPath colours by a column and breaks paths at gaps", {
    p <- plotPath(x, "time", "subject", dimred = "PCA",
        colour.by = "subject")
    expect_length(unique(layer_data(p, 1L)$colour), n_subjects)

    gaps <- unlist(tapply(x$time, x$subject, function(t) diff(sort(t))))
    p <- plotPath(x, "time", "subject", dimred = "PCA", max.gap = 0.3)
    expect_length(unique(layer_data(p, 1L)$group),
        n_subjects + sum(gaps > 0.3))
})
