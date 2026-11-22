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
        ncomponents = 3), "2 components")
    expect_error(plotPath(x, "time", "subject", dimred = "PCA",
        ncomponents = c(1, 9)), "2 components")
    expect_error(plotPath(x, "time", "subject", x.var = "missing",
        y.var = "diversity"), "'x.var' was not found")
    expect_error(plotPath(x, "time", "subject", dimred = "PCA",
        arrows = "yes"), "'arrows'")
    expect_error(plotPath(x, "time", "subject", x.var = rownames(x)[1L],
        y.var = "diversity", assay.type = "clr"),
        "'assay.type' must be one of assayNames")
})

test_that("plotPath joins each subject's samples in time order", {
    p <- plotPath(x, "time", "subject", dimred = "PCA")
    ld <- layer_data(p, 1L)
    expect_length(unique(ld$group), n_subjects)
    expect_equal(ld$x, unname(pca[ord, 1L]))
    expect_equal(ld$y, unname(pca[ord, 2L]))
    expect_identical(get_labs(p)$x, "PCA:1")
    expect_false(is.null(p$layers[[1L]]$geom_params$arrow))
    expect_identical(p$coordinates$ratio, 1)

    p <- plotPath(x, "time", "subject", dimred = "PCA",
        ncomponents = c("PC3", "PC2"), arrows = FALSE)
    ld <- layer_data(p, 1L)
    expect_equal(ld$x, unname(pca[ord, 3L]))
    expect_equal(ld$y, unname(pca[ord, 2L]))
    expect_null(p$layers[[1L]]$geom_params$arrow)
})

test_that("plotPath plots two variables", {
    genus <- "Bacteroides vulgatus et rel."
    # Features come from counts by default, not from the first assay.
    assays(x) <- c(list(other = assay(x, "counts") * 0), as.list(assays(x)))
    p <- plotPath(x, "time", "subject", x.var = genus, y.var = "diversity",
        linewidth = 2)
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
    p <- plotPath(x, "time", "subject", dimred = "PCA", colour.by = "sex")
    expect_setequal(layer_data(p, 1L)$colour, c("#0072B2", "#E69F00"))

    gaps <- unlist(tapply(x$time, x$subject, function(t) diff(sort(t))))
    p <- plotPath(x, "time", "subject", dimred = "PCA", max.gap = 0.3)
    expect_length(unique(layer_data(p, 1L)$group),
        n_subjects + sum(gaps > 0.3))
})

test_that("plotPath breaks paths at visits without a sample", {
    x$visit <- factor(ave(x$time, x$subject, FUN = rank), levels = 1:5)
    gap <- x[, x$visit != "3"]
    p <- plotPath(gap, "visit", "subject", dimred = "PCA")
    expect_length(unique(layer_data(p, 1L)$group), 2L * n_subjects)
    # Lag pairs stay within paths: visits 1-2 and 4-5.
    p <- plotPath(gap, "visit", "subject", x.var = "diversity", lag = 1)
    expect_identical(nrow(p$data), 2L * n_subjects)
    p <- plotPath(gap, "visit", "subject", dimred = "PCA", max.gap = 2)
    expect_length(unique(layer_data(p, 1L)$group), n_subjects)
})

test_that("plotPath draws one path without group", {
    one <- x[, x$subject == x$subject[1L]]
    p <- plotPath(one, "time", dimred = "PCA", origin = "baseline")
    ld <- layer_data(p, 1L)
    expect_length(unique(ld$group), 1L)
    expect_equal(ld$x[1L], 0)
    p <- plotPath(one, "time", x.var = "diversity", lag = 1)
    expect_identical(nrow(p$data), ncol(one) - 1L)
})

test_that("plotPath checks centroids, field, lag and animate", {
    expect_error(plotPath(x, "time", "subject", dimred = "PCA",
        centroids = 1), "'centroids'")
    expect_error(plotPath(x, "time", "subject", dimred = "PCA",
        field = 0), "'field'")
    expect_error(plotPath(x, "time", "subject", dimred = "PCA",
        lag = 1, centroids = TRUE), "cannot be combined")
})

test_that("plotPath draws mean paths per group and time", {
    x$visit <- factor(ave(x$time, x$subject, FUN = rank))
    p <- plotPath(x, "visit", "subject", dimred = "PCA",
        group.by = "sex", centroids = TRUE)
    cen <- p$layers[[3L]]$data
    pc <- SingleCellExperiment::reducedDim(x, "PCA")
    f1 <- x$sex == "female" & x$visit == "1"
    expect_equal(cen$.x[cen$.by == "female" & cen$.t == 1],
        mean(pc[f1, 1]))
    expect_true(all(layer_data(p, 1L)$alpha == 0.3))
    p <- plotPath(x, "visit", "subject", dimred = "PCA",
        centroids = TRUE)
    expect_setequal(p$layers[[3L]]$data$.t, seq_len(nlevels(x$visit)))
})

test_that("plotPath draws a movement field", {
    p <- plotPath(x, "time", "subject", dimred = "PCA", field = 1)
    move <- p$layers[[3L]]$data
    # One cell: the mean of all steps within subjects.
    pc <- SingleCellExperiment::reducedDim(x, "PCA")[, 1]
    ord <- order(x$subject, x$time)
    s <- x$subject[ord]
    d <- diff(pc[ord])[s[-1L] == s[-length(s)]]
    expect_identical(nrow(move), 1L)
    expect_equal(move$.dx, mean(d))
})

test_that("plotPath draws lag plots", {
    p <- plotPath(x, "time", "subject", x.var = "diversity", lag = 1)
    ord <- order(x$subject, x$time)
    s <- x$subject[ord]
    v <- x$diversity[ord]
    keep <- s[-1L] == s[-length(s)]
    expect_equal(sort(p$data$.later), sort(v[-1L][keep]))
    expect_identical(class(p$layers[[1L]]$geom)[1L], "GeomAbline")
    expect_error(plotPath(x, "time", "subject", x.var = "diversity",
        lag = 50), "more than 'lag'")
})

test_that("plotPath animates over time", {
    skip_if_not_installed("gganimate")
    p <- plotPath(x, "time", "subject", dimred = "PCA", animate = TRUE)
    expect_s3_class(p, "gganim")
})
