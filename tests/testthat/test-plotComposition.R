library(SummarizedExperiment)
library(ggplot2)

data("minimalgut", package = "miaTime")
mg <- minimalgut
prop <- prop.table(assay(mg, "counts"), 2L)

data("hitchip1006", package = "miaTime")
# Subjects sampled five times, without the profiles that are stored twice.
h <- hitchip1006[, !duplicated(t(assay(hitchip1006, "counts")))]
h <- h[, h$subject %in% names(which(table(h$subject) == 5))]
h$visit <- factor(ave(h$time, h$subject, FUN = rank))

test_that("plotComposition checks its arguments", {
    expect_error(plotComposition(mg, "Time_hr_num", "StudyIdentifier",
        top = 0), "'top'")
    expect_error(plotComposition(mg, "Time_hr_num", "StudyIdentifier",
        layout = "muller"))
})

test_that("plotComposition stacks proportions that sum to one", {
    p <- plotComposition(mg, "Time_hr_num", "StudyIdentifier", top = 5)
    ld <- layer_data(p, 1L)
    tops <- tapply(ld$ymax, interaction(ld$PANEL, ld$x), max)
    expect_equal(as.vector(tops[!is.na(tops)]), rep(1, ncol(mg)))
    expect_length(unique(ld$PANEL), length(unique(mg$StudyIdentifier)))
})

test_that("plotComposition keeps the most abundant features and Other", {
    p <- plotComposition(mg, "Time_hr_num", "StudyIdentifier", top = 5)
    expected <- names(sort(rowMeans(prop), decreasing = TRUE))[1:5]
    expect_identical(levels(p$data$.feature), c(expected, "Other"))
    # The first six rows are the first sample.
    d <- p$data[1:6, ]
    expect_identical(as.character(d$.feature), c(expected, "Other"))
    expect_equal(d$.value[1:5], unname(prop[expected, 1L]))
    expect_equal(d$.value[6L], 1 - sum(prop[expected, 1L]))

    p <- plotComposition(mg, "Time_hr_num", "StudyIdentifier", top = 50)
    expect_identical(nlevels(p$data$.feature), nrow(mg))
})

test_that("plotComposition reads counts by default and selects features", {
    # The first assay is not counts. s2 is ten times deeper than s1.
    counts <- matrix(c(1, 2, 7, 60, 30, 10), 3L,
        dimnames = list(c("a", "b", "c"), c("s1", "s2")))
    toy <- SummarizedExperiment(
        assays = list(other = counts * 0 + 1, counts = counts),
        colData = DataFrame(t = c(1, 2)))
    expect_error(plotComposition(toy, "t", assay.type = "none"),
        "'assay.type'")
    # Top by mean plotted value: c by proportion (0.4 against 0.35 for a),
    # a by count (30.5 against 8.5).
    df <- plotComposition(toy, "t", top = 1)$data
    expect_identical(levels(df$.feature), c("c", "Other"))
    df <- plotComposition(toy, "t", top = 1, as.relative = FALSE)$data
    expect_identical(levels(df$.feature), c("a", "Other"))
    expect_equal(df$.value[df$.feature == "a"], c(1, 60))
    df <- plotComposition(toy, "t", features = c("b", "a"))$data
    expect_identical(levels(df$.feature), c("b", "a", "Other"))
    expect_error(plotComposition(toy, "t", features = "z"), "'features'")
})

test_that("plotComposition colours features in order and Other grey", {
    p <- plotComposition(mg, "Time_hr_num", "StudyIdentifier", top = 3)
    ld <- layer_data(p, 1L)
    fill <- as.vector(tapply(ld$fill, ld$group, unique))
    expect_identical(fill, c(.okabe_ito[1:3], "grey50"))
    expect_identical(p$theme$legend.position, "bottom")
    # More than six features: viridis, and Other stays grey.
    p <- plotComposition(mg, "Time_hr_num", "StudyIdentifier", top = 8)
    ld <- layer_data(p, 1L)
    fill <- as.vector(tapply(ld$fill, ld$group, unique))
    expect_identical(fill, c(grDevices::hcl.colors(8, "viridis"),
        "grey50"))
})

test_that("plotComposition sums into a feature already named Other", {
    # "Other" is the most abundant row, as after mia::subsetBy(collapse).
    toy <- SummarizedExperiment(
        assays = list(counts = matrix(c(5, 4, 3, 9, 6, 2, 1, 6), 4L,
            dimnames = list(c("a", "b", "c", "Other"), c("s1", "s2")))),
        colData = DataFrame(t = c(1, 2)))
    df <- plotComposition(toy, "t", top = 1)$data
    expect_identical(levels(df$.feature), c("a", "Other"))
    expect_equal(df$.value[df$.feature == "Other"], c(16, 9) / c(21, 15))
    df <- plotComposition(toy, "t", top = 5)$data
    expect_identical(levels(df$.feature), c("a", "b", "c", "Other"))
})

test_that("plotComposition averages groups at each time", {
    p <- plotComposition(h, "visit", "subject", top = 3, level = "group",
        group.by = "sex")
    hp <- prop.table(assay(h, "counts"), 2L)
    f <- levels(p$data$.feature)[1L]
    first <- h$visit == "1" & h$sex == "female"
    d <- p$data[p$data$visit == "1" & p$data$sex == "female" &
        p$data$.feature == f, ]
    expect_equal(d$.value, mean(hp[f, first]))
    expect_length(unique(layer_data(p, 1L)$PANEL), length(unique(h$sex)))
    # Without group.by, all subjects form one group in one panel.
    p <- plotComposition(h, "visit", "subject", top = 3, level = "group")
    d <- p$data[p$data$visit == "1" & p$data$.feature == f, ]
    expect_equal(d$.value, mean(hp[f, h$visit == "1"]))
    expect_s3_class(p$facet, "FacetNull")
})

test_that("plotComposition draws areas and keeps unused visits", {
    p <- plotComposition(mg, "Time_hr_num", "StudyIdentifier",
        layout = "area", alpha = 0.5)
    expect_s3_class(p$layers[[1L]]$geom, "GeomArea")
    expect_true(all(layer_data(p, 1L)$alpha == 0.5))

    h$visit <- factor(h$visit, levels = c(levels(h$visit), "6"))
    p <- plotComposition(h, "visit", "subject")
    expect_true("6" %in% layer_scales(p)$x$get_limits())
})

test_that("plotComposition stacks absolute values and streams", {
    p <- plotComposition(mg, "Time_hr_num", "StudyIdentifier", top = 3,
        as.relative = FALSE)
    d <- p$data
    s1 <- d[d$Time_hr_num == d$Time_hr_num[1] &
        d$StudyIdentifier == d$StudyIdentifier[1], ]
    expect_equal(sum(s1$.value),
        sum(assay(mg, "counts")[, colnames(mg)[mg$Time_hr_num ==
            s1$Time_hr_num[1] & mg$StudyIdentifier ==
            s1$StudyIdentifier[1]]]))
    expect_identical(p$labels$y, "counts")
    p <- plotComposition(mg, "Time_hr_num", "StudyIdentifier", top = 3,
        as.relative = FALSE, layout = "stream")
    d <- p$data
    expect_identical(class(p$layers[[1L]]$geom)[1L], "GeomRibbon")
    key <- list(d$Time_hr_num, d$StudyIdentifier)
    lo <- tapply(d$.ymin, key, min)
    hi <- tapply(d$.ymax, key, max)
    expect_equal(lo[!is.na(lo)], -hi[!is.na(hi)])
    # The first feature is on top, as in bars, and the centred y axis is
    # unlabelled.
    first <- d[d$.feature == levels(d$.feature)[1L], ]
    expect_equal(first$.ymax, as.vector(hi[cbind(
        as.character(first$Time_hr_num), first$StudyIdentifier)]))
    expect_s3_class(p$theme$axis.text.y, "element_blank")
    expect_error(plotComposition(mg, "Time_hr_num", "StudyIdentifier",
        as.relative = "yes"), "'as.relative'")
})
