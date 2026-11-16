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
})

test_that("plotSurvival times state exits and outcomes by group", {
    # s1 leaves A at 5, s2 never, s3 leaves B at 4.
    p <- plotSurvival(toy, "t", "id", event = "state.exit",
        state.col = "state", group.by = "arm")
    x <- p$data[p$data$.by == "x", ]
    expect_equal(x$.time, c(0, 4, 5))
    expect_equal(x$.surv, c(1, 0.5, 0))
    expect_identical(p$labels$colour, "arm")
    df <- fits(event = "outcome", event.col = "sick")
    expect_equal(df$.time[df$.surv < 1], 5)
})
