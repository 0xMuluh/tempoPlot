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
    # Negative blue and positive orange, also when all edges are positive.
    expect_identical(layer_data(p)$colour,
        c("#E69F00", "#0072B2", "#E69F00", "#E69F00"))
    pos <- toy
    S4Vectors::metadata(pos)$e$weight <- 1
    expect_true(all(layer_data(plotNetwork(pos, "e"))$colour == "#E69F00"))
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

test_that("plotNetwork matches undirected edges in either orientation", {
    skip_if_not_installed("igraph")
    rev <- toy
    S4Vectors::metadata(rev)$e$from[3:4] <- c("b", "c")
    S4Vectors::metadata(rev)$e$to[3:4] <- c("a", "a")
    df <- plotNetwork(rev, "e", split.by = "t", diff = TRUE)$layers[[1L]]$data
    ch <- setNames(as.character(df$.change), paste(df$from, df$to))
    expect_identical(ch[["a b"]], "kept")
    expect_identical(nrow(df), 3L)
    # Directed, b -> a is a different edge from a -> b.
    df <- plotNetwork(rev, "e", split.by = "t", diff = TRUE,
        directed = TRUE)$layers[[1L]]$data
    expect_identical(nrow(df), 4L)
    expect_identical(as.character(df$.change[df$from == "b" &
        df$to == "a"]), "only t2")
})

test_that("plotNetwork draws arrows and self-loops of directed edges", {
    skip_if_not_installed("igraph")
    S4Vectors::metadata(toy)$m <- data.frame(from = c("a", "a", "b"),
        to = c("a", "b", "a"), weight = c(5, 2, 1))
    p <- plotNetwork(toy, "m", directed = TRUE)
    expect_false(is.null(p$layers[[1L]]$geom_params$arrow))
    expect_identical(paste(p$layers[[1L]]$data$from, p$layers[[1L]]$data$to),
        c("a b", "b a"))
    loop <- p$layers[[2L]]$data
    expect_identical(class(p$layers[[2L]]$geom)[1L], "GeomPoint")
    expect_identical(c(loop$from, loop$to), c("a", "a"))
    expect_null(plotNetwork(toy, "e")$layers[[1L]]$geom_params$arrow)
})

test_that("plotNetwork names nodes of unnamed colPairs by index", {
    skip_if_not_installed("igraph")
    colnames(toy) <- NULL
    colPair(toy, "net") <- S4Vectors::SelfHits(1L, 2L, nnode = 2L)
    df <- plotNetwork(toy, "net")$layers[[1L]]$data
    expect_identical(c(df$from, df$to), c("1", "2"))
})
