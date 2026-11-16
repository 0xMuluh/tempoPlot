# Extracted from test-plotRecurrence.R:48

# prequel ----------------------------------------------------------------------
library(SummarizedExperiment)
library(ggplot2)
toy <- SummarizedExperiment(
    assays = list(counts = matrix(1, 1L, 5L)),
    colData = DataFrame(id = c("s1", "s2", "s1", "s1", "s2"),
        t = c(3, 0, 0, 1, 2)))
colnames(toy) <- paste0("c", 1:5)
d <- as.matrix(dist(c(10, 20, 11, 15, 21)))
dimnames(d) <- list(colnames(toy), colnames(toy))
S4Vectors::metadata(toy)$d <- d

# test -------------------------------------------------------------------------
toy$v <- factor(toy$t, levels = c(0, 1, 2, 3))
p <- plotRecurrence(toy, "d", "v", "id")
expect_setequal(p$data$xmin[p$data$id == "s2"], c(0.5, 2))
