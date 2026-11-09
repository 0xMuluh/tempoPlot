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
        rank = "missing"), "'rank'")
    expect_error(plotComposition(mg, "Time_hr_num", "StudyIdentifier",
        level = "group"), "group.by")
    expect_error(plotComposition(mg, "Time_hr_num", "StudyIdentifier",
        layout = "stream"))
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

test_that("plotComposition sums features of a rank", {
    p <- plotComposition(h, "visit", "subject", rank = "Phylum", top = 3)
    phyla <- rowsum(assay(h, "counts"), rowData(h)$Phylum)
    phyla <- prop.table(phyla, 2L)
    expected <- names(sort(rowMeans(phyla), decreasing = TRUE))[1:3]
    expect_identical(levels(p$data$.feature), c(expected, "Other"))
    expect_identical(get_labs(p)$fill, "Phylum")
})

test_that("plotComposition averages groups at each time", {
    p <- plotComposition(h, "visit", "subject", rank = "Phylum", top = 3,
        level = "group", group.by = "sex")
    phyla <- prop.table(rowsum(assay(h, "counts"), rowData(h)$Phylum), 2L)
    first <- h$visit == "1" & h$sex == "female"
    d <- p$data[p$data$visit == "1" & p$data$sex == "female" &
        p$data$.feature == "Firmicutes", ]
    expect_equal(d$.value, mean(phyla["Firmicutes", first]))
    expect_length(unique(layer_data(p, 1L)$PANEL), length(unique(h$sex)))
})

test_that("plotComposition draws areas and keeps unused visits", {
    p <- plotComposition(mg, "Time_hr_num", "StudyIdentifier",
        layout = "area", alpha = 0.5)
    expect_s3_class(p$layers[[1L]]$geom, "GeomArea")
    expect_true(all(layer_data(p, 1L)$alpha == 0.5))

    h$visit <- factor(h$visit, levels = c(levels(h$visit), "6"))
    p <- plotComposition(h, "visit", "subject", rank = "Phylum")
    expect_true("6" %in% layer_scales(p)$x$get_limits())
})
