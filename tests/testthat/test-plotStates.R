library(SummarizedExperiment)
library(ggplot2)

# s1: A, A, B and a missing state; s2: B, B, A; s3: A, A, A, B.
toy <- SummarizedExperiment(
    assays = list(counts = matrix(0, 1L, 11L)),
    colData = DataFrame(
        id = rep(c("s1", "s2", "s3"), c(4, 3, 4)),
        state = c("A", "A", "B", NA, "B", "B", "A", "A", "A", "A", "B"),
        arm = rep(c("x", "y", "x"), c(4, 3, 4)),
        day = seq_len(11)))
share <- function(p, s, st) {
    p$data$.share[p$data$.subject == s & p$data$.state == st]
}

test_that("plotStates checks its arguments", {
    expect_error(plotStates(toy, "none", "id"), "'state.col'")
    expect_error(plotStates(toy, "state", NULL), "'group'")
    expect_error(plotStates(toy, "state", "id", group.by = "day"),
        "constant within subjects")
    expect_error(plotStates(toy, "state", "id", layout = "pie"))
    toy$state <- NA
    expect_error(plotStates(toy, "state", "id"), "No sample has a state")
})

test_that("plotStates shares the samples with a state", {
    p <- plotStates(toy, "state", "id")
    expect_equal(share(p, "s1", "A"), 2 / 3)
    expect_equal(share(p, "s2", "B"), 2 / 3)
    expect_equal(share(p, "s3", "A"), 3 / 4)
    expect_equal(p$labels$fill, "share")
})

test_that("plotStates orders subjects by dominant state or clustering", {
    # A first (s3 0.75, s1 0.67), then B (s2); the first at the top.
    p <- plotStates(toy, "state", "id")
    expect_identical(levels(p$data$.subject), c("s2", "s1", "s3"))
    # s1 and s3 have close shares and end up next to each other.
    p <- plotStates(toy, "state", "id", order.by = "cluster")
    expect_equal(abs(diff(match(c("s1", "s3"), levels(p$data$.subject)))),
        1L)
})

test_that("plotStates stacks bars and splits subjects by group", {
    p <- plotStates(toy, "state", "id", layout = "bar", group.by = "arm")
    ld <- layer_data(p)
    # One full bar per subject.
    expect_equal(tapply(ld$xmax, paste(ld$PANEL, ld$y), max), c(1, 1, 1),
        ignore_attr = TRUE)
    expect_equal(p$labels$fill, "state")
    expect_identical(p$data$.by[p$data$.subject == "s2"], c("y", "y"))
    expect_s3_class(p$facet, "FacetGrid")
})
