library(SummarizedExperiment)
library(ggplot2)

# s1 sampled at 0, 1, 3; s2 at 0 and 2.
toy <- SummarizedExperiment(
    assays = list(counts = matrix(1, 1L, 5L)),
    colData = DataFrame(id = c("s1", "s2", "s1", "s1", "s2"),
        t = c(3, 0, 0, 1, 2)))
colnames(toy) <- paste0("c", 1:5)
d <- as.matrix(dist(c(10, 20, 11, 15, 21)))
dimnames(d) <- list(colnames(toy), colnames(toy))
S4Vectors::metadata(toy)$d <- d

test_that("plotRecurrence checks its arguments", {
    expect_error(plotRecurrence(toy, "none", "t", "id"), "'dissimilarity'")
    expect_error(plotRecurrence(toy, "d", "t", NULL), "'group'")
    expect_error(plotRecurrence(toy, "d", "t", "id", subjects = "s9"),
        "'subjects'")
    expect_error(plotRecurrence(toy, "d", "t", "id", threshold = "a"),
        "'threshold'")
})

test_that("plotRecurrence draws all pairs of each subject", {
    df <- plotRecurrence(toy, "d", "t", "id")$data
    expect_identical(as.vector(table(df$id)), c(9L, 4L))
    s1 <- df[df$id == "s1", ]
    # Edges of s1 at 0, 1, 3: -0.5, 0.5, 2, 4.
    expect_setequal(s1$xmin, c(-0.5, 0.5, 2))
    expect_setequal(s1$xmax, c(0.5, 2, 4))
    # Pair (t = 0, t = 3) is samples c3 and c1: |11 - 10| = 1.
    expect_equal(s1$dissimilarity[s1$xmin == -0.5 & s1$ymin == 2], 1)
    expect_true(all(s1$dissimilarity[s1$xmin == s1$ymin] == 0))
    expect_identical(plotRecurrence(toy, "d", "t", "id")$theme$aspect.ratio,
        1)
})

test_that("plotRecurrence selects subjects and marks recurrences", {
    p <- plotRecurrence(toy, "d", "t", "id", subjects = "s2",
        threshold = 1)
    expect_identical(levels(p$data$id), "s2")
    expect_identical(nrow(p$data), 4L)
    ld <- layer_data(p)
    # Only the diagonal and the 20-21 pair are within 1.
    expect_identical(sum(ld$fill == "#0072B2"), 4L)
    ld <- layer_data(plotRecurrence(toy, "d", "t", "id", threshold = 0.5))
    expect_identical(sort(unique(ld$fill)), c("#0072B2", "grey80"))
})

test_that("plotRecurrence handles a single sample, visits and Dates", {
    toy$v <- factor(toy$t, levels = c(0, 1, 2, 3))
    p <- plotRecurrence(toy, "d", "v", "id")
    # Visits 1 and 3 of s2 leave visit 2 empty.
    expect_setequal(p$data$xmin[p$data$id == "s2"], c(0.5, 2.5))
    toy$date <- as.Date("2020-01-01") + toy$t
    p <- plotRecurrence(toy, "d", "date", "id")
    expect_s3_class(p$data$xmin, "Date")
    one <- toy[, 1:4]
    S4Vectors::metadata(one)$d <- d[1:4, 1:4]
    expect_warning(p <- plotRecurrence(one, "d", "t", "id",
        subjects = "s2"), "single sample")
    expect_identical(nrow(p$data), 1L)
})

test_that("plotRecurrence colours dissimilarities on viridis", {
    p <- plotRecurrence(toy, "d", "t", "id")
    ld <- layer_data(p)
    # The diagonal is the lowest dissimilarity.
    expect_true(all(ld$fill[p$data$dissimilarity == 0] == "#440154"))
    expect_identical(p$theme$legend.position, "bottom")
})
