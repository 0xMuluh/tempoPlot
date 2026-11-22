library(SummarizedExperiment)
library(ggplot2)

data("hitchip1006", package = "miaTime")

# Subjects with 3 or more samples. All start at time 0, so age at each
# sample gives start times that differ between subjects.
x <- hitchip1006[, hitchip1006$subject %in%
    names(which(table(hitchip1006$subject) >= 3))]
x$age.time <- x$age + x$time
n_subjects <- length(unique(x$subject))

test_that("plotSwimmer checks its arguments", {
    expect_error(plotSwimmer(x, "time", group = NULL), "'group'")
    expect_error(
        plotSwimmer(x, "time", "subject", colour.by = "missing"),
        "'colour.by'"
    )
    expect_error(plotSwimmer(x, "time", "subject", order.by = "last"))
    x_chr <- x
    x_chr$time <- as.character(x_chr$time)
    expect_error(plotSwimmer(x_chr, "time", "subject"), "factor")
})

test_that("plotSwimmer draws a point per sample and a lane per subject", {
    p <- plotSwimmer(x, "time", "subject")
    expect_s3_class(p, "ggplot")
    expect_identical(nrow(layer_data(p, 2L)), ncol(x))
    expect_length(unique(layer_data(p, 1L)$group), n_subjects)

    colData(x)[["sampling time"]] <- x$time
    p <- plotSwimmer(x, "sampling time", "subject", size = 4)
    expect_true(all(layer_data(p, 2L)$size == 4))
})

test_that("plotSwimmer orders subjects by first time, first at the top", {
    ld <- layer_data(plotSwimmer(x, "age.time", "subject"), 2L)
    # tapply() returns rows bottom to top.
    start <- rev(tapply(ld$x, ld$y, min))
    expect_gt(length(unique(start)), 1L)
    expect_true(all(diff(start) >= 0))
})

test_that("plotSwimmer orders subjects by number of samples", {
    ld <- layer_data(plotSwimmer(x, "time", "subject", order.by = "n"), 2L)
    n <- rev(as.vector(table(ld$y)))
    expect_gt(length(unique(n)), 1L)
    expect_true(all(diff(n) <= 0))
})

test_that("plotSwimmer breaks lanes at gaps larger than max.gap", {
    gaps <- unlist(tapply(x$time, x$subject, function(t) diff(sort(t))))
    p <- plotSwimmer(x, "time", "subject", max.gap = 5)
    expect_length(
        unique(layer_data(p, 1L)$group),
        n_subjects + sum(gaps > 5)
    )
})

test_that("plotSwimmer colours points by colour.by", {
    p <- plotSwimmer(x, "time", "subject", colour.by = "bmi_group")
    expect_length(
        unique(layer_data(p, 2L)$colour),
        length(unique(x$bmi_group))
    )
    expect_true(all(layer_data(p, 2L)$colour %in% c(.okabe_ito, "grey50")))
    p <- plotSwimmer(x, "time", "subject")
    expect_length(unique(layer_data(p, 2L)$colour), 1L)
    expect_identical(p$theme$legend.position, "bottom")
    p <- plotSwimmer(x, "time", "subject", colour.by = "age")
    expect_s3_class(p$scales$get_scales("colour"), "ScaleContinuous")
})

test_that("plotSwimmer keeps subjects with a single sample", {
    expect_warning(
        p <- plotSwimmer(hitchip1006, "time", "subject"),
        "single sample"
    )
    ld <- layer_data(p, 2L)
    expect_identical(nrow(ld), ncol(hitchip1006))
    expect_length(unique(ld$y), length(unique(hitchip1006$subject)))
})

test_that("plotSwimmer draws visits in level order with unused levels", {
    x$visit <- factor(x$time, levels = c(-1, sort(unique(x$time))))
    ld <- layer_data(plotSwimmer(x, "visit", "subject"), 2L)
    expect_identical(
        sort(unique(as.integer(ld$x))),
        sort(unique(as.integer(x$visit)))
    )
})

test_that("plotSwimmer hides the labels of more than 60 subjects", {
    few <- x[, x$subject %in% unique(x$subject)[1:60]]
    p <- plotSwimmer(few, "time", "subject")
    expect_false(inherits(p$theme$axis.text.y, "element_blank"))
    many <- hitchip1006[, hitchip1006$subject %in%
        names(which(table(hitchip1006$subject) >= 2))]
    p <- plotSwimmer(many, "time", "subject")
    expect_match(class(p$theme$axis.text.y)[1L], "element_blank")
})

test_that("plotSwimmer draws Date times on a continuous axis", {
    x$date <- as.Date("2000-01-01") + round(x$time * 30)
    p <- plotSwimmer(x, "date", "subject")
    expect_s3_class(layer_scales(p)$x, "ScaleContinuousDate")
    expect_equal(sort(layer_data(p, 2L)$x), sort(as.numeric(x$date)))
})

test_that("plotSwimmer checks missing and events", {
    expect_error(plotSwimmer(x, "time", "subject", missing = "skip"))
    expect_error(
        plotSwimmer(x, "time", "subject", missing = "show"),
        "requires 'max.gap'"
    )
    expect_error(
        plotSwimmer(x, "time", "subject", events = data.frame(end = 1)),
        "'start'"
    )
    expect_error(
        plotSwimmer(x, "time", "subject",
            events = data.frame(start = "a")),
        "type of 'time.col'"
    )
    expect_error(
        plotSwimmer(x, "time", "subject",
            events = data.frame(start = 2, end = 1)),
        "before"
    )
    x$visit <- factor(x$time)
    expect_error(
        plotSwimmer(x, "visit", "subject",
            events = data.frame(start = "not a visit")),
        "levels"
    )
})

# Subjects 1-3 have visits 1-3; subject 2 misses visit 2 and subject 3
# misses visit 3 (dropout).
v <- x[, seq_len(9L)]
v$id <- rep(c("s1", "s2", "s3"), each = 3L)
v$visit <- factor(rep(c("v1", "v2", "v3"), 3L), levels = c("v1", "v2", "v3"))
v <- v[, -c(5L, 9L)]

test_that("plotSwimmer breaks lanes at missed visits", {
    ld <- layer_data(plotSwimmer(v, "visit", "id"), 1L)
    # s1 one lane, s2 two lanes, s3 one lane
    expect_length(unique(ld$group), 4L)
})

test_that("plotSwimmer shows missing samples", {
    p <- plotSwimmer(v, "visit", "id", missing = "show")
    geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1))
    gap <- layer_data(p, which(geoms == "GeomSegment"))
    expect_identical(nrow(gap), 1L)
    expect_equal(as.numeric(c(gap$x, gap$xend)), c(1, 3))
    missed <- layer_data(p, which(geoms == "GeomPoint")[2L])
    expect_identical(nrow(missed), 2L)
    expect_setequal(as.numeric(missed$x), c(2, 3))

    gaps <- unlist(tapply(x$time, x$subject, function(t) diff(sort(t))))
    p <- plotSwimmer(x, "time", "subject", max.gap = 5, missing = "show")
    geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1))
    expect_identical(
        nrow(layer_data(p, which(geoms == "GeomSegment"))),
        sum(gaps > 5)
    )
})

test_that("plotSwimmer drops subjects with missing samples", {
    expect_message(
        p <- plotSwimmer(v, "visit", "id", missing = "drop"),
        "Removed 2 subjects"
    )
    expect_identical(levels(p$data$id), "s1")
    expect_error(
        plotSwimmer(v[, v$id != "s1"], "visit", "id", missing = "drop"),
        "All subjects"
    )
    # Factor subjects whose labels look like the codes of other subjects:
    # subject "3" (code 2) has a gap, so subject "2" must stay.
    v$id <- factor(c("2", "2", "2", "3", "3", "9", "9"))
    expect_message(
        p <- plotSwimmer(v, "visit", "id", missing = "drop"),
        "Removed 2 subjects"
    )
    expect_identical(levels(p$data$id), "2")
})

test_that("plotSwimmer draws events on subject rows and across subjects", {
    events <- data.frame(
        subject = c(NA, NA, "217", "218"),
        start = c(1, 2, 0, 3),
        end = c(NA, 4, 1, NA),
        label = c("a", "b", "c", "d")
    )
    # Events without a subject span all subjects.
    p <- plotSwimmer(x, "time", "subject", events = events[1:2, -1L])
    geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1))
    expect_identical(geoms[[1L]], "GeomRect")
    expect_true(is.infinite(layer_data(p, 1L)$ymin))
    expect_true("GeomVline" %in% geoms)
    expect_identical(p$labels$fill, "event")

    # Events of a subject are drawn on its row.
    p <- plotSwimmer(x, "time", "subject", events = events[3:4, ])
    rows <- levels(p$data$subject)
    band <- layer_data(p, 1L)
    expect_equal(band$ymin, match("217", rows) - 0.4)
    marker <- layer_data(p, length(p$layers))
    expect_equal(marker$y, match("218", rows) + 0.3)
    expect_identical(c(band$fill, marker$fill), .okabe_ito[1:2])

    # Events of absent subjects are not drawn.
    p <- plotSwimmer(x, "time", "subject",
        events = data.frame(subject = "none", start = 1))
    expect_length(p$layers, 2L)
})

test_that("plotSwimmer draws visit events at their levels", {
    p <- plotSwimmer(v, "visit", "id",
        events = data.frame(start = "v2", end = "v3"))
    band <- layer_data(p, 1L)
    expect_equal(as.numeric(c(band$xmin, band$xmax)), c(2, 3))
})

test_that("plotSwimmer draws one band per experiment", {
    cd <- DataFrame(id = c("s1", "s1", "s2", "s2"), t = c(0, 1, 0, 2),
        row.names = paste0("p", 1:4))
    e1 <- SummarizedExperiment(matrix(1, 1L, 4L,
        dimnames = list("f", paste0("p", 1:4))))
    e2 <- SummarizedExperiment(matrix(1, 1L, 2L,
        dimnames = list("g", c("p1", "p4"))))
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        list(a = e1, b = e2), colData = cd)
    p <- plotSwimmer(mae, "t", "id")
    df <- p$data
    expect_identical(nrow(df), 6L)
    expect_identical(as.character(df$.primary[df$.experiment == "b"]),
        c("p1", "p4"))
    # Both subjects start at 0; ties keep their order, s1 on top.
    expect_equal(sort(unique(df$.row[df$id == "s1"])), 2)
    expect_equal(unique(df$.y[df$.experiment == "a" & df$id == "s1"]),
        2.3)
    expect_error(plotSwimmer(mae, "none", "id"), "'time.col'")
})
