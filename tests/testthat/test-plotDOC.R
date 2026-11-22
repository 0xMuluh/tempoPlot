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
    expect_error(plotDOC(toy, "none", "d"), "'overlap'")
    expect_error(plotDOC(toy, "o", "none"), "'dissimilarity'")
    expect_error(plotDOC(toy, "o", "d", pairs = "within"), "'group'")
})

test_that("plotDOC fades points above 1000 pairs and draws the LOWESS", {
    set.seed(1)
    many <- data.frame(.overlap = runif(1500), .dissimilarity = runif(1500))
    p <- .plot_doc(many, "o", "d")
    expect_equal(unique(layer_data(p, 1L)$alpha), 2 / 3)
    curve <- lowess(many$.overlap, many$.dissimilarity)
    expect_equal(layer_data(p, 2L)$y, curve$y[order(curve$x)])
    expect_identical(unique(layer_data(p, 2L)$colour), .okabe_ito[1L])
    expect_identical(p$theme$legend.position, "bottom")
    p <- .plot_doc(many[1:100, ], "o", "d", alpha = 0.5)
    expect_equal(unique(layer_data(p, 1L)$alpha), 0.5)
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
