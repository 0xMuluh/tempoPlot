library(SummarizedExperiment)
library(ggplot2)

# Three subjects at visits 1 and 2; s3's second sample is closer to s1.
toy <- SummarizedExperiment(
    assays = list(counts = matrix(1, 1L, 6L)),
    colData = DataFrame(id = rep(c("s1", "s2", "s3"), each = 2),
        v = factor(rep(c("1", "2"), 3))))
colnames(toy) <- paste0("c", 1:6)
pos <- c(0, 1, 10, 11, 20, 2)
d <- as.matrix(dist(pos))
dimnames(d) <- list(colnames(toy), colnames(toy))
S4Vectors::metadata(toy)$d <- d

test_that("plotFingerprint checks its arguments", {
    expect_error(plotFingerprint(toy, "none", "v", "id", "1", "2"),
        "'dissimilarity'")
    expect_error(plotFingerprint(toy, "d", "v", "id", "9", "2"),
        "'reference.time'")
})

test_that("plotFingerprint identifies subjects by their own sample", {
    p <- plotFingerprint(toy, "d", "v", "id", "1", "2")
    df <- p$data
    expect_identical(df$.subject, c("s1", "s2", "s3"))
    expect_equal(df$own, c(1, 1, 18))
    expect_equal(df$other, c(9, 9, 2))
    expect_identical(df$identified, c(TRUE, TRUE, FALSE))
    expect_match(p$labels$subtitle, "2 of 3")
})
