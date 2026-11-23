library(SummarizedExperiment)
library(ggplot2)

# Four samples of two subjects; arm splits s1 from s2.
toy <- SummarizedExperiment(
    assays = list(counts = matrix(1, 1L, 4L)),
    colData = DataFrame(id = c("s1", "s1", "s2", "s2"),
        arm = c("x", "x", "y", "x")))
colnames(toy) <- paste0("c", 1:4)
d <- matrix(c(
    0, 1, 5, 6,
    1, 0, 7, 8,
    5, 7, 0, 2,
    6, 8, 2, 0), 4L, dimnames = list(colnames(toy), colnames(toy)))
S4Vectors::metadata(toy)$d <- as.dist(d)

test_that("plotIndividuality checks its arguments", {
    expect_error(plotIndividuality(toy, "none", "id"), "'dissimilarity'")
    expect_error(plotIndividuality(toy, "d", NULL), "'group'")
    expect_error(plotIndividuality(toy, "d", "id", layout = "hist"))
    S4Vectors::metadata(toy)$bad <- matrix(0, 2L, 2L)
    expect_error(plotIndividuality(toy, "bad", "id"), "sample by sample")
    one <- toy
    one$id <- paste0("s", 1:4)
    expect_error(plotIndividuality(one, "d", "id"), "more than one sample")
})

test_that("plotIndividuality splits pairs into within and between", {
    df <- plotIndividuality(toy, "d", "id")$data
    expect_identical(nrow(df), 6L)
    expect_setequal(df$dissimilarity[df$pair == "within"], c(1, 2))
    expect_setequal(df$dissimilarity[df$pair == "between"], c(5, 6, 7, 8))
    expect_setequal(layer_data(plotIndividuality(toy, "d", "id"))$colour,
        c("#0072B2", "#E69F00"))
})

test_that("plotIndividuality matches samples by name", {
    S4Vectors::metadata(toy)$d <- as.dist(d[4:1, 4:1])
    df <- plotIndividuality(toy, "d", "id")$data
    expect_setequal(df$dissimilarity[df$pair == "within"], c(1, 2))
})

test_that("plotIndividuality keeps pairs within facet.by", {
    p <- plotIndividuality(toy, "d", "id", layout = "box", facet.by = "arm")
    # Pairs in arm x: c1-c2 (within), c1-c4, c2-c4 (between).
    expect_setequal(p$data$dissimilarity, c(1, 6, 8))
    expect_true(all(p$data$arm == "x"))
    expect_s3_class(p$facet, "FacetWrap")
    expect_identical(class(p$layers[[1L]]$geom)[1L], "GeomBoxplot")
    p <- plotIndividuality(toy, "d", "id", layout = "violin")
    expect_identical(class(p$layers[[1L]]$geom)[1L], "GeomViolin")
})

test_that("plotIndividuality splits between pairs by group.by", {
    p <- plotIndividuality(toy, "d", "id", layout = "box", group.by = "arm")
    # c1-c4 and c2-c4 share arm x; c1-c3 and c2-c3 do not.
    df <- p$data
    expect_identical(levels(df$pair),
        c("within", "between, same arm", "between, different arm"))
    expect_setequal(df$dissimilarity[df$pair == "between, same arm"], c(6, 8))
    expect_setequal(df$dissimilarity[df$pair == "between, different arm"],
        c(5, 7))
})

test_that("plotIndividuality layout time summarises pairs per window", {
    toy$w <- c(1, 1, 1, 2)
    expect_error(plotIndividuality(toy, "d", "id", layout = "time"),
        "'time.col'")
    toy$chr <- c("a", "a", "a", "b")
    expect_error(plotIndividuality(toy, "d", "id", layout = "time",
        time.col = "chr"), "numeric, Date")
    p <- plotIndividuality(toy, "d", "id", layout = "time", time.col = "w")
    geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1))
    expect_identical(unname(geoms), c("GeomRibbon", "GeomLine"))
    # Window 1 holds c1-c2 (within, 1), c1-c3 (5) and c2-c3 (7); the
    # single sample of window 2 gives no pair.
    band <- layer_data(p, 1L)
    expect_equal(band$x, c(1, 1))
    expect_equal(band$ymin[order(band$group)], c(1, 5.5))
    expect_equal(band$ymax[order(band$group)], c(1, 6.5))
    line <- layer_data(p, 2L)
    expect_equal(line$y[order(line$group)], c(1, 6))
    expect_identical(line$colour[order(line$group)], .okabe_ito[1:2])
    expect_identical(p$labels$y, "d")
})
