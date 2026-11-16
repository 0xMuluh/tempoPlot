# Extracted from test-plotLagHeatmap.R:21

# prequel ----------------------------------------------------------------------
library(SummarizedExperiment)
library(ggplot2)
set.seed(1)
v <- rnorm(20)
toy <- SummarizedExperiment(
    assays = list(counts = rbind(a = c(v[-1], 0), b = v)),
    colData = DataFrame(t = 1:20, y = v))

# test -------------------------------------------------------------------------
df <- plotLagHeatmap(toy, "y", "t", max.lag = 2)$data
r <- function(f, k) df$.r[df$.feature == f & df$.lag == k]
expect_equal(r("b", 0), 1)
expect_equal(r("a", -1), 1)
