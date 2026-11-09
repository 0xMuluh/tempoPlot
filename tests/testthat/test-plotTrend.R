library(SummarizedExperiment)
library(SingleCellExperiment)
library(ggplot2)

data("minimalgut", package = "miaTime")
x <- minimalgut
n_subjects <- length(unique(x$StudyIdentifier))
n_times <- length(unique(x$Time_hr_num))

test_that("plotTrend checks its arguments", {
    expect_error(
        plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier",
            group.by = "missing"),
        "'group.by'"
    )
    expect_error(
        plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier",
            facet.by = "missing"),
        "'facet.by'"
    )
    expect_error(
        plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier",
            layout = "invalid")
    )
    expect_error(
        plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier",
            smooth = "invalid")
    )
    expect_error(
        plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier",
            show.individuals = "yes"),
        "'show.individuals'"
    )
    expect_error(
        plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier",
            baseline = "yes"),
        "should be one of"
    )
    expect_error(
        plotTrend(x, "Acetate", "Time_hr_num"),
        "duplicated"
    )
    x_chr <- x
    x_chr$Time_hr_num <- as.character(x_chr$Time_hr_num)
    expect_error(
        plotTrend(x_chr, "Acetate", "Time_hr_num", "StudyIdentifier"),
        "factor"
    )
    expect_error(
        plotTrend(x, "missing", "Time_hr_num", "StudyIdentifier"),
        "not found"
    )
})

test_that("plotTrend draws a line per subject and a point per sample", {
    p <- plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier")
    expect_s3_class(p, "ggplot")
    expect_length(p$layers, 2L)
    expect_length(unique(layer_data(p, 1L)$group), n_subjects)
    expect_identical(nrow(layer_data(p, 2L)), ncol(x))

    # ... goes to geom_line
    p_lw <- plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier",
        linewidth = 2.5)
    expect_true(all(layer_data(p_lw, 1L)$linewidth == 2.5))
})

test_that("plotTrend supports layout = 'box'", {
    p <- plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier",
        layout = "box")
    expect_s3_class(p, "ggplot")
    expect_length(p$layers, 3L)
    expect_identical(length(unique(layer_data(p, 1L)$x)), n_times)
    expect_length(unique(layer_data(p, 2L)$group), n_subjects)
    expect_identical(nrow(layer_data(p, 3L)), ncol(x))

    # ... goes to geom_boxplot
    p_box <- plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier",
        layout = "box", notch = TRUE)
    expect_true(all(layer_data(p_box, 1L)$notchupper >
        layer_data(p_box, 1L)$middle))
})

test_that("show.individuals = FALSE hides subject lines and points", {
    # In layout = 'line', smooth is left
    p_line <- plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier",
        layout = "line", show.individuals = FALSE, smooth = "mean")
    expect_length(p_line$layers, 1L)
    expect_identical(nrow(layer_data(p_line, 1L)), n_times)

    # In layout = 'box', boxes are left
    p_box <- plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier",
        layout = "box", show.individuals = FALSE)
    expect_length(p_box$layers, 1L)
    expect_identical(length(unique(layer_data(p_box, 1L)$x)), n_times)
})

test_that("smooth = 'mean' draws the mean per time point on top", {
    p <- plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier",
        smooth = "mean")
    expect_length(p$layers, 3L)
    ld_mean <- layer_data(p, 3L)
    expect_identical(nrow(ld_mean), n_times)

    # Check that calculated values match actual means
    expected_means <- as.vector(tapply(x$Acetate, x$Time_hr_num, mean))
    expect_equal(ld_mean$y, expected_means)
})

test_that("individuals fade under a summary, which is drawn thicker", {
    p <- plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier")
    expect_true(all(is.na(layer_data(p, 1L)$alpha)))
    p <- plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier",
        smooth = "mean")
    expect_true(all(layer_data(p, 1L)$alpha == 0.3))
    expect_true(all(layer_data(p, 2L)$alpha == 0.3))
    expect_gt(layer_data(p, 3L)$linewidth[1L],
        layer_data(p, 1L)$linewidth[1L])
    p <- plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier",
        layout = "box")
    expect_true(all(layer_data(p, 2L)$alpha == 0.3))

    # Individual points replace the outlier points of the boxes.
    data("hitchip1006", package = "miaTime")
    h <- hitchip1006[, hitchip1006$subject %in%
        names(which(table(hitchip1006$subject) >= 3))]
    h$visit <- factor(ave(h$time, h$subject, FUN = rank))
    outliers <- function(show) {
        p <- plotTrend(h, "diversity", "visit", "subject", layout = "box",
            show.individuals = show)
        sum(lengths(layer_data(p, 1L)$outliers))
    }
    expect_identical(outliers(TRUE), 0L)
    expect_gt(outliers(FALSE), 0L)
})

test_that("smooth = 'loess' overlays a loess curve", {
    p <- plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier",
        smooth = "loess")
    expect_length(p$layers, 3L)
    ld_loess <- layer_data(p, 3L)
    expect_gt(nrow(ld_loess), 0L)
    expect_s3_class(p$layers[[3L]]$stat, "StatSmooth")
})

test_that("baseline expresses values relative to the first sample", {
    # Value at the earliest time of each subject line.
    first_y <- function(p) {
        ld <- layer_data(p, 1L)
        vapply(split(ld, ld$group), function(d) d$y[which.min(d$x)],
            numeric(1))
    }
    expected <- c(difference = 0, fold = 1, log2fold = 0)
    labels <- c(difference = "Change in Acetate",
        fold = "Fold change in Acetate",
        log2fold = "log2 fold change in Acetate")
    for (b in names(expected)) {
        p <- plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier",
            baseline = b)
        expect_equal(unname(first_y(p)), rep(expected[[b]], n_subjects))
        expect_identical(get_labs(p)$y, labels[[b]])
    }

    # Last value of one bioreactor against its first, from the data.
    a <- x[, x$StudyIdentifier == "Bioreactor A"]
    a <- a[, order(a$Time_hr_num)]
    ratio <- a$Acetate[ncol(a)] / a$Acetate[1L]
    last_y <- function(b) {
        ld <- layer_data(plotTrend(a, "Acetate", "Time_hr_num",
            "StudyIdentifier", baseline = b), 1L)
        ld$y[which.max(ld$x)]
    }
    expect_equal(last_y("difference"), a$Acetate[ncol(a)] - a$Acetate[1L])
    expect_equal(last_y("fold"), ratio)
    expect_equal(last_y("log2fold"), log2(ratio))

    # Samples stored in reverse time order, starting after time 0.
    x_rev <- x[, rev(which(x$Time_hr_num >= 48))]
    p_rev <- plotTrend(x_rev, "Acetate", "Time_hr_num", "StudyIdentifier",
        baseline = "difference")
    expect_equal(unname(first_y(p_rev)), rep(0, n_subjects))

    # Ratios need positive first values.
    x_zero <- x
    x_zero$Acetate[which.min(x_zero$Time_hr_num)] <- 0
    expect_error(plotTrend(x_zero, "Acetate", "Time_hr_num",
        "StudyIdentifier", baseline = "fold"), "positive")
    expect_no_error(plotTrend(x_zero, "Acetate", "Time_hr_num",
        "StudyIdentifier", baseline = "difference"))
})

test_that("group.by colours subjects and smooths", {
    p <- plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier",
        group.by = "StudyIdentifier", smooth = "mean")
    expect_identical(
        length(unique(layer_data(p, 1L)$colour)),
        n_subjects
    )
    expect_identical(
        length(unique(layer_data(p, 2L)$colour)),
        n_subjects
    )
    expect_identical(
        length(unique(layer_data(p, 3L)$colour)),
        n_subjects
    )
})

test_that("facet.by facets the plot", {
    p <- plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier",
        facet.by = "Pertubation")
    expect_s3_class(p$facet, "FacetWrap")
})

test_that("max.gap breaks trajectories", {
    # Gaps in minimalgut Time_hr_num: consecutive points are 4, 16, etc.
    p_gap <- plotTrend(x, "Acetate", "Time_hr_num", "StudyIdentifier",
        max.gap = 10)
    expect_gt(
        length(unique(layer_data(p_gap, 1L)$group)),
        n_subjects
    )
})

test_that("plotTrend resolves assay features and reducedDim", {
    taxon <- rownames(x)[[1L]]
    p_feat <- plotTrend(x, taxon, "Time_hr_num", "StudyIdentifier",
        assay.type = "counts")
    expect_identical(nrow(layer_data(p_feat, 2L)), ncol(x))

    x_sce <- as(x, "SingleCellExperiment")
    reducedDim(x_sce, "PCA") <- matrix(
        seq_len(ncol(x_sce) * 2),
        ncol = 2
    )
    colnames(reducedDim(x_sce, "PCA")) <- c("PC1", "PC2")
    p_pca <- plotTrend(x_sce, "PCA:1", "Time_hr_num", "StudyIdentifier")
    expect_identical(nrow(layer_data(p_pca, 2L)), ncol(x_sce))
})

test_that("plotTrend handles group = NULL", {
    # Subset to distinct times
    x_one <- x[, x$StudyIdentifier == "Bioreactor A"]
    p_none <- plotTrend(x_one, "Acetate", "Time_hr_num")
    expect_s3_class(p_none, "ggplot")
    expect_length(unique(layer_data(p_none, 1L)$group), 1L)

    # Without group, all samples form one series in either layout.
    p_box_none <- plotTrend(x_one, "Acetate", "Time_hr_num", layout = "box")
    expect_length(unique(layer_data(p_box_none, 2L)$group), 1L)
})

test_that("plotTrend handles Date and factor time types", {
    x_date <- x
    colData(x_date)$date <- as.Date("2020-01-01") +
        as.integer(x_date$Time_hr_num)
    p_date <- plotTrend(x_date, "Acetate", "date", "StudyIdentifier")
    expect_s3_class(layer_scales(p_date)$x, "ScaleContinuousDate")

    x_fac <- x
    t <- colData(x_fac)$Time_hr_num
    colData(x_fac)$visit <- factor(
        paste0("H", t),
        levels = paste0("H", sort(unique(t)))
    )
    p_fac <- plotTrend(x_fac, "Acetate", "visit", "StudyIdentifier")
    expect_identical(
        sort(unique(as.integer(layer_data(p_fac, 2L)$x))),
        sort(unique(as.integer(x_fac$visit)))
    )
})
