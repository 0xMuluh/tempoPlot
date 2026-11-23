skip_if_not_installed("patchwork")
library(SummarizedExperiment)
library(ggplot2)

# Subject a at times 1, 2, 3 and 5, subject b at 1 to 3; three features.
counts <- matrix(c(
    1, 2, 3, 4, 5, 6, 7,
    1, 1, 1, 1, 1, 1, 1,
    8, 7, 6, 5, 4, 3, 2), nrow = 3L, byrow = TRUE,
    dimnames = list(c("f1", "f2", "f3"), paste0("s", 1:7)))
se <- SummarizedExperiment(list(counts = counts), colData = DataFrame(
    subject = rep(c("a", "b"), c(4L, 3L)),
    time = c(1, 2, 3, 5, 1, 2, 3),
    load = c(10, 20, 30, 50, 1, 2, 3),
    state = c("x", "x", "y", "y", "z", "z", "z")))

x_range <- function(p) {
    ggplot_build(p)$layout$panel_params[[1L]]$x.range
}

test_that("plotTimeline checks the subject", {
    expect_error(plotTimeline(se, "time", "subject", subject = "c"),
        "'subject' must be one value")
    expect_error(plotTimeline(se, "time", "subject"),
        "'subject' must be one value")
    expect_error(plotTimeline(se, "time", "subject", subject = c("a", "b")),
        "'subject' must be one value")
    expect_error(plotTimeline(se, "time", "subject", subject = "a",
        state.col = "colour"), "'state.col'")
})

test_that("plotTimeline stacks one track per layer", {
    p <- plotTimeline(se, "time", "subject", subject = "a")
    expect_s3_class(p, "patchwork")
    expect_length(p, 1L)

    p <- plotTimeline(se, "time", "subject", subject = "a",
        y = c("load", "f1"), state.col = "state",
        events = data.frame(start = 2, end = 3, label = "course"))
    expect_length(p, 5L)
    # Composition, the two lines, states and events from top to bottom.
    expect_identical(p[[2L]]$labels$y, "load")
    expect_identical(p[[3L]]$labels$y, "f1")
    expect_identical(p[[4L]]$labels$fill, "state")
    expect_identical(layer_data(p[[5L]])$xmin, 2)
})

test_that("plotTimeline draws the subject only", {
    p <- plotTimeline(se, "time", "subject", subject = "b", y = "load",
        state.col = "state")
    expect_setequal(p[[1L]]$data$time, 1:3)
    expect_identical(layer_data(p[[2L]])$y, c(1, 2, 3))
    expect_identical(as.character(p[[3L]]$data$.state), rep("z", 3L))
    # Events of other subjects are left out, with their track.
    events <- data.frame(subject = c("a", "b"), start = c(2, 1),
        end = c(3, NA))
    p <- plotTimeline(se, "time", "subject", subject = "b", events = events)
    expect_length(p, 2L)
    expect_identical(layer_data(p[[2L]])$xintercept, 1)
    p <- plotTimeline(se, "time", "subject", subject = "b",
        events = events[1L, ])
    expect_length(p, 1L)
})

test_that("plotTimeline shares the x axis between tracks", {
    p <- plotTimeline(se, "time", "subject", subject = "a", y = "load",
        state.col = "state")
    # Half the spacing before 1 and after 5.
    ranges <- lapply(seq_along(p), function(i) x_range(p[[i]]))
    expect_identical(ranges, rep(list(c(0.5, 6)), 3L))
    breaks <- lapply(seq_along(p), function(i)
        ggplot_build(p[[i]])$layout$panel_params[[1L]]$x$breaks)
    expect_true(all(vapply(breaks, identical, logical(1), breaks[[1L]])))
    # Only the bottom track shows the time axis.
    hidden <- vapply(seq_along(p), function(i)
        inherits(p[[i]]$theme$axis.title.x, "element_blank"), logical(1))
    expect_identical(hidden, c(TRUE, TRUE, FALSE))

    # Events outside the samples widen the axis.
    p <- plotTimeline(se, "time", "subject", subject = "a",
        events = data.frame(start = -2, end = 0))
    expect_identical(x_range(p[[1L]]), c(-2, 6))
    expect_identical(x_range(p[[2L]]), c(-2, 6))
})

test_that("plotTimeline keeps all visits on every track", {
    v <- se[, se$subject == "b"]
    v$visit <- factor(c("V1", "V2", "V4"), levels = paste0("V", 1:4))
    p <- plotTimeline(v, "visit", state.col = "state",
        events = data.frame(start = "V3"))
    ranges <- lapply(seq_along(p), function(i) x_range(p[[i]]))
    expect_identical(ranges, rep(list(c(0.5, 4.5)), 3L))
    labels <- ggplot_build(p[[3L]])$layout$panel_params[[1L]]$x$get_labels()
    expect_identical(labels, paste0("V", 1:4))
    expect_equal(layer_data(p[[3L]])$xintercept, 3)
})
