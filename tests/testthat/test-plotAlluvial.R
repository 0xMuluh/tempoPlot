library(SummarizedExperiment)
library(ggplot2)

# s1 and s2: A -> B; s3: A -> A; s4: B and no second visit.
toy <- SummarizedExperiment(
    assays = list(counts = matrix(1, 1L, 7L)),
    colData = DataFrame(
        id = c("s1", "s1", "s2", "s2", "s3", "s3", "s4"),
        v = factor(c("v1", "v2", "v1", "v2", "v1", "v2", "v1")),
        state = c("A", "B", "A", "B", "A", "A", "B")))
flows <- function(df) {
    w <- tapply(as.character(df$.state), list(df$.flow, df$.time), c)
    setNames(df$.n[df$.time == "v1"], apply(w, 1L, paste, collapse = ">"))
}

test_that("plotAlluvial checks its arguments", {
    skip_if_not_installed("ggalluvial")
    expect_error(suppressWarnings(plotAlluvial(toy, "none", "v", "id")),
        "'state.col'")
    expect_error(plotAlluvial(toy, "state", "v", NULL), "'group'")
    expect_error(suppressWarnings(plotAlluvial(toy, "state", "v", "id",
        min.flow = 0)), "'min.flow'")
    expect_error(suppressWarnings(plotAlluvial(toy, "state", "v", "id",
        missing = "break")))
})

test_that("plotAlluvial counts subjects per path", {
    skip_if_not_installed("ggalluvial")
    p <- suppressWarnings(plotAlluvial(toy, "state", "v", "id"))
    f <- flows(p$data)
    expect_equal(f[["A>B"]], 2)
    expect_equal(f[["A>A"]], 1)
    expect_equal(f[["B>missing"]], 1)
    expect_identical(levels(p$data$.state), c("A", "B", "missing"))
    expect_s3_class(p, "ggplot")
})

test_that("plotAlluvial drops subjects and hides small flows", {
    skip_if_not_installed("ggalluvial")
    expect_message(
        p <- suppressWarnings(plotAlluvial(toy, "state", "v", "id",
            missing = "drop")),
        "Removed 1 subjects")
    expect_false("missing" %in% levels(p$data$.state))
    expect_message(
        p <- suppressWarnings(plotAlluvial(toy, "state", "v", "id",
            min.flow = 2)),
        "Hid 2 subjects")
    expect_identical(names(flows(p$data)), "A>B")
    expect_error(suppressWarnings(plotAlluvial(toy, "state", "v", "id",
        min.flow = 5)), "No path")
})
