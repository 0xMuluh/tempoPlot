library(SummarizedExperiment)
library(ggplot2)

toy <- SummarizedExperiment(
    assays = list(counts = matrix(1, 1L, 5L)),
    colData = DataFrame(id = c("s1", "s1", "s2", "s3", "s3"),
        a = c(1, 1, 2, 3, 3), b = c(5, 5, 4, NA, NA),
        sex = c("f", "f", "m", "m", "m"), chr = "x"))

test_that("plotSubjectScatter checks its arguments", {
    expect_error(plotSubjectScatter(toy, "a", "b", NULL), "'group'")
    expect_error(plotSubjectScatter(toy, "none", "b", "id"), "'x.stat'")
    expect_error(plotSubjectScatter(toy, "a", "chr", "id"), "numeric")
    expect_error(plotSubjectScatter(toy, "a", "b", "id", quadrants = 1),
        "'quadrants'")
    toy$a[1] <- 9
    expect_error(plotSubjectScatter(toy, "a", "b", "id"), "constant")
})

test_that("plotSubjectScatter draws one point per subject", {
    p <- plotSubjectScatter(toy, "a", "b", "id", colour.by = "sex")
    # s3 has no b and is dropped.
    expect_identical(p$data$id, c("s1", "s2"))
    expect_identical(nrow(layer_data(p)), 2L)
    expect_identical(p$labels$colour, "sex")
})

test_that("plotSubjectScatter draws quadrants", {
    p <- plotSubjectScatter(toy, "a", "b", "id", quadrants = TRUE)
    expect_equal(layer_data(p, 1L)$xintercept, 1.5)
    expect_equal(layer_data(p, 2L)$yintercept, 4.5)
    p <- plotSubjectScatter(toy, "a", "b", "id", quadrants = c(0, 1))
    expect_equal(layer_data(p, 1L)$xintercept, 0)
    expect_length(plotSubjectScatter(toy, "a", "b", "id")$layers, 1L)
})
