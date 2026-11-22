library(SummarizedExperiment)
library(ggplot2)

# s1 at 0, 2, 5; s2 at 1, 3; s3 at 0, 4, 5.
toy <- SummarizedExperiment(
    assays = list(counts = matrix(c(
        0, 1, 1, 0, 0, 2, 0, 0), 1L,
        dimnames = list("f", paste0("c", 1:8)))),
    colData = DataFrame(
        id = c("s1", "s1", "s1", "s2", "s2", "s3", "s3", "s3"),
        t = c(0, 2, 5, 1, 3, 0, 4, 5),
        state = c("A", "A", "B", "A", "A", "B", "A", "B"),
        sick = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
        arm = c("x", "x", "x", "y", "y", "x", "x", "x")))
fits <- function(...) {
    df <- plotSurvival(toy, "t", "id", ...)$data
    df[df$.time > 0 | df$.censored > 0, ]
}

test_that("plotSurvival checks its arguments", {
    expect_error(plotSurvival(toy, "t", NULL), "'group'")
    expect_error(plotSurvival(toy, "t", "id", event = "arrival"),
        "'feature'")
    expect_error(plotSurvival(toy, "t", "id", event = "state.exit"),
        "'state.col'")
    expect_error(plotSurvival(toy, "t", "id", event = "outcome",
        event.col = "state"), "logical")
    expect_error(plotSurvival(toy, "t", "id", event = "censored"))
})

test_that("plotSurvival defines dropout from the last sampling time", {
    # Only s2 ends before time 5, after 2 units of follow-up.
    df <- fits()
    expect_equal(df$.time, c(2, 5))
    expect_equal(df$.surv, c(2 / 3, 2 / 3))
    expect_equal(df$.censored, c(0L, 2L))
})

test_that("plotSurvival times arrival and loss of a feature", {
    # Absent at first sample: s1 (arrives at 2), s2 (never, censored at 2).
    df <- fits(event = "arrival", feature = "f")
    expect_equal(df$.time, c(2))
    expect_equal(df$.surv, 0.5)
    # Present at first sample: s3, lost at 4.
    df <- fits(event = "loss", feature = "f")
    expect_equal(df$.time, 4)
    expect_equal(df$.surv, 0)
    # Above detection = 1, only s3's 2 is present: s1 never arrives.
    df <- fits(event = "arrival", feature = "f", detection = 1)
    expect_true(all(df$.surv == 1))
    expect_error(fits(event = "arrival", feature = "f", assay.type = "clr"),
        "'assay.type' must be one of assayNames")
    expect_error(fits(event = "arrival", feature = "f", detection = NA),
        "'detection'")
})

test_that("plotSurvival times state exits and outcomes by group", {
    # s1 leaves A at 5, s2 never, s3 leaves B at 4.
    p <- plotSurvival(toy, "t", "id", event = "state.exit",
        state.col = "state", group.by = "arm")
    x <- p$data[p$data$.by == "x", ]
    expect_equal(x$.time, c(0, 4, 5))
    expect_equal(x$.surv, c(1, 0.5, 0))
    expect_identical(p$labels$colour, "arm")
    expect_setequal(layer_data(p, 2L)$colour, .okabe_ito[1:2])
    expect_setequal(layer_data(p, 1L)$fill, .okabe_ito[1:2])
    expect_identical(p$theme$legend.position, "bottom")
    df <- fits(event = "outcome", event.col = "sick")
    expect_equal(df$.time[df$.surv < 1], 5)
})

test_that("plotSurvival draws the confidence band behind the curve", {
    # s1 leaves A at 5, s2 is censored at 2, s3 leaves B at 4.
    p <- plotSurvival(toy, "t", "id", event = "state.exit",
        state.col = "state")
    fit <- survival::survfit(survival::Surv(c(5, 2, 4), c(1, 0, 1)) ~ 1)
    expect_equal(p$data$.lower, c(1, fit$lower))
    expect_equal(p$data$.upper, c(1, fit$upper))
    expect_equal(p$data$.next, c(2, 4, 5, 5))
    expect_identical(class(p$layers[[1L]]$geom)[1L], "GeomRect")
})
