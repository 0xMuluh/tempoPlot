# Extracted from test-plotDirection.R:33

# prequel ----------------------------------------------------------------------
library(SummarizedExperiment)
library(SingleCellExperiment)
library(ggplot2)
toy <- SingleCellExperiment(
    assays = list(counts = matrix(1, 1L, 5L)),
    colData = DataFrame(id = c("s1", "s1", "s1", "s2", "s2"),
        t = c(0, 1, 2, 0, 1)))
reducedDim(toy, "PC") <- cbind(PC1 = c(0, 1, 1, 0, -2),
    PC2 = c(0, 0, 3, 0, 0))

# test -------------------------------------------------------------------------
reducedDim(toy, "PC")[5, ] <- 0
df <- plotDirection(toy, "t", "id", "PC")$data
expect_identical(nrow(df), 2L)
toy$arm <- c("a", "a", "b", "b", "b")
p <- plotDirection(toy, "t", "id", "PC", group.by = "arm")
expect_setequal(p$data$arm, c("a", "b"))
