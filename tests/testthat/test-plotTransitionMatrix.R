library(SummarizedExperiment)
library(ggplot2)

# Two subjects: s1 goes A, A, B; s2 goes B, NA, A, A.
toy <- SummarizedExperiment(
    assays = list(counts = matrix(1, 1L, 7L)),
    colData = DataFrame(
        id = c("s1", "s1", "s1", "s2", "s2", "s2", "s2"),
        t = c(1, 2, 3, 1, 2, 3, 4),
        state = c("A", "A", "B", "B", NA, "A", "A"),
        arm = c("x", "x", "x", "y", "y", "y", "y")))
cell <- function(df, from, to, col = "n") {
    df[[col]][df$from == from & df$to == to]
}

test_that("plotTransitionMatrix checks its arguments", {
    expect_error(plotTransitionMatrix(toy, "state", "t", group = NULL),
        "'group'")
    expect_error(plotTransitionMatrix(toy, "none", "t", "id"), "state.col")
    expect_error(plotTransitionMatrix(toy, "t", "t", "id"), "factor")
    expect_error(plotTransitionMatrix(toy, "state", "t", "id",
        as.relative = "row"), "'as.relative'")
    toy$state <- NA_character_
    expect_error(plotTransitionMatrix(toy, "state", "t", "id"),
        "consecutive")
})

test_that("plotTransitionMatrix counts transitions within subjects", {
    df <- plotTransitionMatrix(toy, "state", "t", "id",
        as.relative = FALSE)$data
    # A->A twice (s1, s2), A->B once; B->NA and NA->A are skipped, and
    # s1's last B does not lead into s2.
    expect_identical(cell(df, "A", "A"), 2L)
    expect_identical(cell(df, "A", "B"), 1L)
    expect_identical(cell(df, "B", "A"), 0L)
    expect_identical(sum(df$n), 3L)
})

test_that("plotTransitionMatrix shows shares of each row", {
    df <- plotTransitionMatrix(toy, "state", "t", "id")$data
    expect_equal(cell(df, "A", "A", "value"), 2 / 3)
    expect_equal(cell(df, "A", "B", "value"), 1 / 3)
    expect_true(all(is.na(df$value[df$from == "B"])))
    expect_identical(cell(df, "A", "A", ".label"), "0.67")
})

test_that("plotTransitionMatrix splits panels by facet.by", {
    p <- plotTransitionMatrix(toy, "state", "t", "id", as.relative = FALSE,
        facet.by = "arm")
    df <- p$data
    expect_identical(cell(df[df$arm == "x", ], "A", "A"), 1L)
    expect_identical(cell(df[df$arm == "y", ], "A", "A"), 1L)
    expect_s3_class(p$facet, "FacetWrap")
})

test_that("plotTransitionMatrix keeps factor levels", {
    toy$state <- factor(toy$state, levels = c("B", "A", "C"))
    df <- plotTransitionMatrix(toy, "state", "t", "id")$data
    expect_identical(levels(df$from), c("B", "A", "C"))
    expect_identical(nrow(df), 9L)
})

test_that("plotTransitionMatrix skips transitions across gaps", {
    # s1 misses visit v3: its B at v4 does not follow A at v2.
    toy$v <- factor(c("v1", "v2", "v4", "v1", "v2", "v3", "v4"),
        levels = paste0("v", 1:4))
    df <- plotTransitionMatrix(toy, "state", "v", "id",
        as.relative = FALSE)$data
    expect_identical(cell(df, "A", "B"), 0L)
    expect_identical(sum(df$n), 2L)
    df <- plotTransitionMatrix(toy, "state", "v", "id", as.relative = FALSE,
        max.gap = 2)$data
    expect_identical(cell(df, "A", "B"), 1L)
    # Numeric times are skipped only beyond max.gap.
    toy$t <- c(1, 2, 9, 1, 2, 3, 4)
    df <- plotTransitionMatrix(toy, "state", "t", "id", as.relative = FALSE,
        max.gap = 3)$data
    expect_identical(cell(df, "A", "B"), 0L)
    expect_error(plotTransitionMatrix(toy, "state", "t", "id", max.gap = 0),
        "'max.gap'")
})

test_that("plotTransitionMatrix fills shares on the sequential scale", {
    p <- plotTransitionMatrix(toy, "state", "t", "id")
    ld <- layer_data(p, 1L)
    # Rows without transitions (from B) are NA and grey.
    expect_true(all(ld$fill[is.na(p$data$value)] == "grey50"))
    expect_identical(ld$fill[which.max(p$data$value)], "#FDE725")
    expect_identical(p$theme$legend.position, "bottom")
})
