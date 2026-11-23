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

test_that("plotNetwork draws a symmetric matrix of undirected edges", {
    p <- plotNetwork(toy, "e", split.by = "t", layout = "matrix")
    df <- p$data
    # t1 has a-b and b-c, t2 a-b and a-c: each filled in both orientations.
    expect_setequal(paste(df$.split, df$from, df$to), c("t1 a b", "t1 b a",
        "t1 b c", "t1 c b", "t2 a b", "t2 b a", "t2 a c", "t2 c a"))
    expect_equal(df$weight[df$.split == "t1" & df$from == "c"], -0.5)
    expect_s3_class(p$facet, "FacetWrap")
    # Nodes sorted by name: a in the first column and the top row.
    expect_identical(p$scales$get_scales("x")$get_limits(), c("a", "b", "c"))
    expect_identical(p$scales$get_scales("y")$get_limits(), c("c", "b", "a"))
    # Negative blue, positive red.
    ld <- layer_data(p)
    rgb <- col2rgb(ld$fill[df$weight == -0.5])
    expect_true(all(rgb["blue", ] > rgb["red", ]))
    rgb <- col2rgb(ld$fill[which.max(df$weight)])
    expect_gt(rgb["red", 1L], rgb["blue", 1L])
    # A table listing both orientations collapses to one tile each.
    both <- toy
    S4Vectors::metadata(both)$e <- data.frame(from = c("a", "b"),
        to = c("b", "a"), weight = 0.3)
    expect_identical(nrow(plotNetwork(both, "e", layout = "matrix")$data), 2L)
})

test_that("plotNetwork matrix keeps directions and the diagonal", {
    # Illustrative interactions: self-limitation on the diagonal, a
    # inhibits b, b promotes a.
    S4Vectors::metadata(toy)$glv <- data.frame(from = c("a", "b", "a", "b"),
        to = c("a", "b", "b", "a"), weight = c(-1, -2, -0.5, 0.4),
        p = c(0.001, 0.2, 0.01, 0.3))
    p <- plotNetwork(toy, "glv", directed = TRUE, layout = "matrix",
        p.col = "p")
    df <- p$data
    expect_identical(nrow(df), 4L)
    expect_equal(df$weight[df$from == "b" & df$to == "a"], 0.4)
    expect_identical(sum(df$from == df$to), 2L)
    dots <- p$layers[[2L]]$data
    expect_identical(paste(dots$from, dots$to), c("a a", "a b"))
    expect_identical(p$labels$caption, "Dots: p < 0.05")
    expect_error(plotNetwork(toy, "glv", layout = "matrix", p.col = "from"),
        "'p.col'")
    expect_error(plotNetwork(toy, "glv", p.col = "p"), "layout")
    expect_error(plotNetwork(toy, "glv", layout = "matrix", p.col = "p",
        threshold = 2), "'threshold'")
    # Undirected, b-a and a-b are the same pair with two weights.
    expect_error(plotNetwork(toy, "glv", layout = "matrix"), "one row per")
})

test_that("plotNetwork matrix fills gained, lost and kept edges", {
    p <- plotNetwork(toy, "e", split.by = "t", diff = TRUE,
        layout = "matrix")
    ch <- setNames(as.character(p$data$.change),
        paste(p$data$from, p$data$to))
    expect_identical(ch[c("a b", "b a", "b c", "a c")],
        c(`a b` = "kept", `b a` = "kept", `b c` = "only t1",
            `a c` = "only t2"))
    expect_null(p$facet$params$facets)
    expect_identical(p$labels$fill, "t")
})
