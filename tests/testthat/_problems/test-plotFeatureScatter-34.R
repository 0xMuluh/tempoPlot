# Extracted from test-plotFeatureScatter.R:34

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
p <- plotFeatureScatter(toy, "mean", "variance", fit = "power")
curve <- p$layers[[2L]]$data
expect_match(unique(curve$.label), "slope 2.00")
