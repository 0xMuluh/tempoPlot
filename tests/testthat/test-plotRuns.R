library(SummarizedExperiment)
library(ggplot2)

# s1 at 0, 1, 2, 4, 5 and s2 at 0, 2, 3. Feature f is present in s1 at
# 0, 1, 4, 5 and in s2 at 3; g never.
toy <- SummarizedExperiment(
    assays = list(counts = matrix(c(
        1, 1, 0, 2, 3, 0, 0, 1,
        0, 0, 0, 0, 0, 0, 0, 0), 2L, byrow = TRUE,
        dimnames = list(c("f", "g"), paste0("c", 1:8)))),
    colData = DataFrame(
        id = c("s1", "s1", "s1", "s1", "s1", "s2", "s2", "s2"),
        t = c(0, 1, 2, 4, 5, 0, 2, 3),
        state = c("A", "A", "B", NA, "B", "B", "A", "A"),
        arm = c("x", "x", "x", "x", "x", "y", "y", "y")))
# Runs of the values v of toy, sorted by subject and time.
runs <- function(v, max.gap = NULL) {
    df <- data.frame(id = toy$id, .t = toy$t, .i = seq_len(ncol(toy)))
    .find_runs(v, df, "id", max.gap, visit = FALSE)
}
# Kaplan-Meier steps after time 0 of one panel.
steps <- function(p, state) {
    df <- p$data[p$data$.state == state & p$data$.time > 0, ]
    df[, c(".time", ".surv", ".censored")]
}

test_that("plotRuns checks its arguments", {
    expect_error(plotRuns(toy, "t", NULL), "'group'")
    expect_error(plotRuns(toy, "t", "id", state.col = "none"), "'state.col'")
    expect_error(plotRuns(toy, "t", "id", features = "h"), "'features'")
    expect_error(plotRuns(toy, "t", "id", assay.type = "clr"),
        "'assay.type'")
    expect_error(plotRuns(toy, "t", "id", max.gap = 0), "'max.gap'")
    expect_error(plotRuns(toy, "t", "id", layout = "line"))
    toy$state <- NA
    expect_error(plotRuns(toy, "t", "id", state.col = "state"),
        "No sample")
})

test_that("runs last until the next run and are censored at the end", {
    f <- c("p", "p", "a", "p", "p", "a", "a", "p")
    # s1: p 0-1 until a at 2, a 2 until p at 4, p 4-5 censored;
    # s2: a 0-2 until p at 3, p 3 censored.
    expect_equal(runs(f), data.frame(.state = c("p", "a", "p", "a", "p"),
        .time = c(2, 2, 1, 3, 0), .status = c(1L, 1L, 0L, 1L, 0L),
        .i = c(1L, 3L, 4L, 6L, 8L)))
    # Gaps above 1.5 end runs, censored: s1 a at 2, s2 a at 0.
    r <- runs(f, max.gap = 1.5)
    expect_equal(r$.time, c(2, 0, 1, 0, 1, 0))
    expect_equal(r$.status, c(1L, 0L, 0L, 0L, 1L, 0L))
    # Missing values are skipped: s1 A 0-1 until B at 2, B 2-5 censored.
    expect_equal(runs(toy$state)$.time, c(2, 3, 2, 1))
    expect_equal(runs(toy$state)$.status, c(1L, 0L, 1L, 0L))
})

test_that("plotRuns pools features into present and absent panels", {
    p <- plotRuns(toy, "t", "id", features = "f")
    expect_identical(levels(p$data$.state), c("present", "absent"))
    # Present runs: 2 (seen), 1 and 0 (censored). The censored 0 leaves 2
    # at risk: 1 censored, then 2 ends.
    expect_equal(steps(p, "present")$.surv, c(1, 0))
    expect_equal(steps(p, "present")$.censored, c(1L, 0L))
    # Absent runs: 2 and 3 of f, seen, and 5 and 3 of g, censored.
    p <- plotRuns(toy, "t", "id")
    expect_equal(steps(p, "absent")$.time, c(2, 3, 5))
    expect_equal(steps(p, "absent")$.surv, c(0.75, 0.5, 0.5))
    # Above detection = 1, f is present only at times 4 and 5 of s1.
    p <- plotRuns(toy, "t", "id", features = "f", detection = 1)
    expect_equal(steps(p, "present")$.time, 1)
})

test_that("plotRuns draws state runs by group and on log axes", {
    # A: s1 2 (seen) and s2 1 (censored); B: s1 3 (censored), s2 2.
    p <- plotRuns(toy, "t", "id", state.col = "state")
    expect_equal(steps(p, "A")$.surv, c(1, 0))
    expect_equal(steps(p, "B")$.surv, c(0.5, 0.5))
    # Each run takes arm from its first sample.
    p <- plotRuns(toy, "t", "id", state.col = "state", group.by = "arm")
    expect_setequal(p$data$.by, c("x", "y"))
    expect_equal(p$labels$colour, "arm")
    # ccdf drops times and shares of zero: A keeps time 1, B times 2, 3.
    p <- plotRuns(toy, "t", "id", state.col = "state", layout = "ccdf")
    expect_equal(layer_data(p)$x, log10(c(1, 2, 3)))
    expect_equal(layer_data(p)$y, log10(c(1, 0.5, 0.5)))
})

test_that("plotRuns ends runs at skipped visits", {
    # Visits 1, 2, 4 and 5: g is absent in two runs of one step, censored.
    se <- toy["g", 1:4]
    se$visit <- factor(c(1, 2, 4, 5), levels = 1:5)
    p <- plotRuns(se, "visit", "id")
    expect_equal(p$data$.time, c(0, 1))
    expect_equal(p$data$.censored, c(0L, 2L))
})
