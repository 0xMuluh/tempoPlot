library(SummarizedExperiment)
library(ggplot2)

toy <- SummarizedExperiment(
    assays = list(counts = matrix(1, 1L, 4L)),
    colData = DataFrame(id = c("s1", "s1", "s2", "s2")))
colnames(toy) <- paste0("c", 1:4)
o <- matrix(c(0, .9, .5, .4, .9, 0, .6, .3, .5, .6, 0, .8, .4, .3, .8, 0),
    4L, dimnames = list(colnames(toy), colnames(toy)))
S4Vectors::metadata(toy)$o <- o
S4Vectors::metadata(toy)$d <- 1 - o

test_that("plotDOC checks its arguments", {
    expect_error(plotDOC(toy, "none", "d"), "'dissimilarity'")
    expect_error(plotDOC(toy, "o", "d", pairs = "within"), "'group'")
})

test_that("plotDOC draws each pair once", {
    df <- plotDOC(toy, "o", "d")$data
    expect_identical(nrow(df), 6L)
    expect_equal(df$.dissimilarity, 1 - df$.overlap)
    w <- plotDOC(toy, "o", "d", group = "id", pairs = "within")$data
    expect_setequal(w$.overlap, c(0.9, 0.8))
    b <- plotDOC(toy, "o", "d", group = "id", pairs = "between")$data
    expect_setequal(b$.overlap, c(0.5, 0.4, 0.6, 0.3))
})
