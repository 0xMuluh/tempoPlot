# Extracted from test-plotFeatureScatter.R:33

# prequel ----------------------------------------------------------------------
library(SummarizedExperiment)
library(ggplot2)
toy <- SummarizedExperiment(
    assays = list(counts = matrix(c(
        1, 3, 0, 4,
        10, 30, 20, 40,
        100, 300, 200, 400), 3L, byrow = TRUE,
        dimnames = list(c("a", "b", "c"), paste0("c", 1:4)))),
    colData = DataFrame(arm = c("x", "x", "y", "y")),
    rowData = DataFrame(score = c(1, 2, 3), label = c("p", "q", "r")))

# test -------------------------------------------------------------------------
se <- SummarizedExperiment(
        assays = list(counts = outer(c(1, 10, 100), c(1, 3, 2, 4))))
p <- plotFeatureScatter(se, "mean", "variance", fit = "power")
