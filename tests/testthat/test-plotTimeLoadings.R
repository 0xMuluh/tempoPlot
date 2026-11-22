library(SummarizedExperiment)
library(ggplot2)

toy <- SummarizedExperiment(assays = list(counts = matrix(1, 1L, 2L)))
S4Vectors::metadata(toy) <- list(
    tem = list(Phi_hat = cbind(PC1 = c(0, 1, 2), PC2 = c(3, 2, 1)),
        time_Phi = c(0, 5, 10), r_square = c(0.6, 0.25)),
    par = list(Fac = list(matrix(1, 2L, 2L), matrix(1, 3L, 2L),
        matrix(1:8, 4L, 2L, dimnames = list(c(1, 2, 4, 8), NULL)))),
    mat = data.frame(time = c(1, 2), a = c(0.1, 0.2)),
    bad = "x")

test_that("plotTimeLoadings checks its arguments", {
    expect_error(plotTimeLoadings(toy, "none"), "'source'")
    expect_error(plotTimeLoadings(toy, "bad"), "not a TEMPTED")
    expect_error(plotTimeLoadings(toy, "tem", ncomponents = 3),
        "'ncomponents'")
})

test_that("plotTimeLoadings reads TEMPTED results", {
    p <- plotTimeLoadings(toy, "tem")
    expect_equal(p$data$.time, c(0, 5, 10, 0, 5, 10))
    expect_equal(p$data$.value, c(0, 1, 2, 3, 2, 1))
    expect_identical(levels(p$data$.component),
        c("PC1 (R\u00b2 0.60)", "PC2 (R\u00b2 0.25)"))
    p <- plotTimeLoadings(toy, "tem", ncomponents = 1)
    expect_identical(nlevels(p$data$.component), 1L)
    p <- plotTimeLoadings(toy, "tem", ncomponents = "PC2")
    expect_equal(p$data$.value, c(3, 2, 1))
})

test_that("plotTimeLoadings reads PARAFAC models and matrices", {
    p <- plotTimeLoadings(toy, "par")
    expect_equal(unique(p$data$.time), c(1, 2, 4, 8))
    expect_equal(p$data$.value, 1:8)
    p <- plotTimeLoadings(toy, "mat")
    expect_equal(p$data$.time, c(1, 2))
    expect_identical(levels(p$data$.component), "a")
})

test_that("plotTimeLoadings plots a real TEMPTED fit", {
    skip_if_not_installed("tempted")
    data("minimalgut", package = "miaTime")
    tse <- minimalgut
    counts <- t(assay(tse, "counts"))
    invisible(capture.output(res <- tempted::tempted_all(counts,
        timepoint = tse$Time_hr_num,
        subjectID = as.character(tse$StudyIdentifier), r = 2,
        threshold = 1, pseudo = 0.5, do_ratio = FALSE)))
    S4Vectors::metadata(tse)$tempted <- res
    p <- plotTimeLoadings(tse, "tempted")
    expect_equal(unique(p$data$.time), res$time_Phi)
    expect_equal(p$data$.value, as.vector(res$Phi_hat))
})
