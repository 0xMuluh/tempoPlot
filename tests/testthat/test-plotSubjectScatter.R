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

test_that("plotSubjectScatter uses the shared theme and colour scales", {
    p <- plotSubjectScatter(toy, "a", "b", "id", colour.by = "sex")
    expect_identical(layer_data(p)$colour, c("#0072B2", "#E69F00"))
    expect_identical(p$theme$legend.position, "bottom")
    p <- plotSubjectScatter(toy, "a", "b", "id", colour.by = "a")
    expect_identical(layer_data(p)$colour, c("#440154", "#FDE725"))
})

test_that("plotSubjectScatter draws quadrants", {
    p <- plotSubjectScatter(toy, "a", "b", "id", quadrants = TRUE)
    expect_equal(layer_data(p, 1L)$xintercept, 1.5)
    expect_equal(layer_data(p, 2L)$yintercept, 4.5)
    p <- plotSubjectScatter(toy, "a", "b", "id", quadrants = c(0, 1))
    expect_equal(layer_data(p, 1L)$xintercept, 0)
    expect_length(plotSubjectScatter(toy, "a", "b", "id")$layers, 1L)
})

test_that("plotSubjectScatter draws subject loadings of a decomposition", {
    S4Vectors::metadata(toy) <- list(
        tem = list(A_hat = cbind(PC1 = c(s1 = 0.1, s2 = 0.2, s3 = 0.3),
            PC2 = c(0.4, 0.5, 0.6), PC3 = c(0.7, 0.8, 0.9)),
            r_square = c(0.5, 0.3, 0.1)),
        par = list(Fac = list(matrix(1:4, 2L,
            dimnames = list(c("s2", "s9"), NULL)))),
        anon = list(A_hat = matrix(1, 3L, 2L)))
    p <- plotSubjectScatter(toy, group = "id", colour.by = "sex",
        source = "tem", ncomponents = c(1, 3))
    expect_identical(as.character(p$data$.subject), c("s1", "s2", "s3"))
    expect_equal(p$data$.y, c(0.7, 0.8, 0.9))
    expect_identical(p$labels$x, "PC1 (R\u00b2 0.50)")
    expect_identical(p$labels$y, "PC3 (R\u00b2 0.10)")
    q <- plotSubjectScatter(toy, group = "id", source = "tem",
        ncomponents = c("PC1", "PC3"))
    expect_equal(q$data$.y, p$data$.y)
    # Only subjects of x are drawn.
    p <- plotSubjectScatter(toy, group = "id", source = "par")
    expect_identical(as.character(p$data$.subject), "s2")
    expect_equal(c(p$data$.x, p$data$.y), c(1, 3))
    expect_error(plotSubjectScatter(toy, group = "id", source = "none"),
        "'source'")
    expect_error(plotSubjectScatter(toy, group = "id", source = "tem",
        ncomponents = c(1, 4)), "'ncomponents'")
    expect_error(plotSubjectScatter(toy, group = "id", source = "anon"),
        "row names")
})
