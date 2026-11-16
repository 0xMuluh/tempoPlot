# Extracted from test-plotRanks.R:31

# prequel ----------------------------------------------------------------------
library(SummarizedExperiment)
library(ggplot2)
toy <- SummarizedExperiment(
    assays = list(counts = matrix(c(
        9, 1, 5,
        5, 9, 1,
        1, 5, 9,
        0, 0, 0), 4L, byrow = TRUE,
        dimnames = list(c("a", "b", "c", "d"), paste0("c", 1:3)))),
    colData = DataFrame(t = c(1, 2, 3)),
    rowData = DataFrame(genus = c("g1", "g1", "g2", "g3")))

# test -------------------------------------------------------------------------
df <- plotRanks(toy, "t", top = 1, rank = "genus")$data
