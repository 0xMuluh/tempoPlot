library(SummarizedExperiment)
library(ggplot2)

# s1 sees a, then b, then a and c; s2 sees nothing, then b.
toy <- SummarizedExperiment(
    assays = list(counts = matrix(c(
        1, 0, 2, 0, 0,
        0, 3, 0, 0, 1,
        0, 0, 1, NA, 0), 3L, byrow = TRUE,
        dimnames = list(c("a", "b", "c"), paste0("c", 1:5)))),
    colData = DataFrame(id = c("s1", "s1", "s1", "s2", "s2"),
        t = c(1, 2, 3, 1, 2), arm = c("x", "x", "x", "y", "y")),
    rowData = DataFrame(genus = c("g1", "g1", "g2")))

test_that("plotAccumulation checks its arguments", {
    expect_error(plotAccumulation(toy, "t", "id", assay.type = "none"),
        "'assay.type'")
    expect_error(plotAccumulation(toy, "t", "id", detection = "a"),
        "'detection'")
    expect_error(plotAccumulation(toy, "t", "id", rank = "none"), "'rank'")
})

test_that("plotAccumulation counts distinct features over time", {
    df <- plotAccumulation(toy, "t", "id")$data
    expect_equal(df$.n[df$.subject == "s1"], c(1, 2, 3))
    expect_equal(df$.n[df$.subject == "s2"], c(0, 1))
    df <- plotAccumulation(toy, "t", "id", detection = 1)$data
    expect_equal(df$.n[df$.subject == "s1"], c(0, 1, 2))
    df <- plotAccumulation(toy, "t", "id", rank = "genus")$data
    expect_equal(df$.n[df$.subject == "s1"], c(1, 1, 2))
    expect_identical(plotAccumulation(toy, "t", "id",
        rank = "genus")$labels$y, "cumulative genus")
})

test_that("plotAccumulation adds group means", {
    p <- plotAccumulation(toy, "t", "id", group.by = "arm")
    expect_identical(class(p$layers[[2L]]$stat)[1L], "StatSummary")
    expect_identical(p$labels$colour, "arm")
})
