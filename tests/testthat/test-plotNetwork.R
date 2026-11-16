library(SummarizedExperiment)
library(SingleCellExperiment)
library(ggplot2)

toy <- SingleCellExperiment(assays = list(counts = matrix(1, 3L, 2L,
    dimnames = list(c("a", "b", "c"), c("s1", "s2")))))
S4Vectors::metadata(toy)$e <- data.frame(
    from = c("a", "b", "a", "a"), to = c("b", "c", "b", "c"),
    weight = c(0.8, -0.5, 0.7, 0.9), t = c("t1", "t1", "t2", "t2"))

test_that("plotNetwork checks its arguments", {
    skip_if_not_installed("igraph")
    expect_error(plotNetwork(toy, "none"), "'edges'")
    expect_error(plotNetwork(toy, "e", split.by = "none"), "'split.by'")
    expect_error(plotNetwork(toy, "e", diff = TRUE), "two levels")
})

test_that("plotNetwork draws edges from metadata in panels", {
    skip_if_not_installed("igraph")
    p <- plotNetwork(toy, "e", split.by = "t")
    df <- p$layers[[1L]]$data
    expect_identical(nrow(df), 4L)
    expect_identical(as.character(df$.sign), c("positive", "negative",
        "positive", "positive"))
    # The same node sits at the same place in every panel.
    a <- df[df$from == "a", ]
    expect_identical(length(unique(a$.x)), 1L)
    expect_s3_class(p$facet, "FacetWrap")
})

test_that("plotNetwork marks gained, lost and kept edges", {
    skip_if_not_installed("igraph")
    df <- plotNetwork(toy, "e", split.by = "t", diff = TRUE)$layers[[1L]]$data
    ch <- setNames(as.character(df$.change), paste(df$from, df$to))
    expect_identical(ch[["a b"]], "kept")
    expect_identical(ch[["b c"]], "only t1")
    expect_identical(ch[["a c"]], "only t2")
    expect_identical(nrow(df), 3L)
})

test_that("plotNetwork reads rowPairs", {
    skip_if_not_installed("igraph")
    rowPair(toy, "net") <- S4Vectors::SelfHits(c(1L, 2L), c(2L, 3L),
        nnode = 3L, value = c(0.5, -0.2))
    df <- plotNetwork(toy, "net")$layers[[1L]]$data
    expect_identical(paste(df$from, df$to), c("a b", "b c"))
    expect_equal(df$weight, c(0.5, -0.2))
})
