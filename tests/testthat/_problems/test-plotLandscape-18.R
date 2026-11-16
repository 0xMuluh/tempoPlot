# Extracted from test-plotLandscape.R:18

# prequel ----------------------------------------------------------------------
library(SummarizedExperiment)
library(ggplot2)
set.seed(1)
toy <- SummarizedExperiment(
    assays = list(counts = matrix(1, 1L, 60L)),
    colData = DataFrame(v = c(rnorm(30, -2, 0.3), rnorm(30, 2, 0.3)),
        arm = rep(c("a", "b"), 30)))

# test -------------------------------------------------------------------------
p <- plotLandscape(toy, "v", bw = 0.3)
