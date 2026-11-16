# Extracted from test-plotFingerprint.R:15

# prequel ----------------------------------------------------------------------
library(SummarizedExperiment)
library(ggplot2)
toy <- SummarizedExperiment(
    assays = list(counts = matrix(1, 1L, 6L)),
    colData = DataFrame(id = rep(c("s1", "s2", "s3"), each = 2),
        v = rep(c("1", "2"), 3)))
colnames(toy) <- paste0("c", 1:6)
pos <- c(0, 1, 10, 11, 20, 2)
S4Vectors::metadata(toy)$d <- as.matrix(dist(pos))

# test -------------------------------------------------------------------------
expect_error(plotFingerprint(toy, "none", "v", "id", "1", "2"),
        "'dissimilarity'")
