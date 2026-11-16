# Extracted from test-plotPath.R:85

# prequel ----------------------------------------------------------------------
library(SummarizedExperiment)
library(SingleCellExperiment)
library(ggplot2)
data("hitchip1006", package = "miaTime")
x <- hitchip1006[, !duplicated(t(assay(hitchip1006, "counts")))]
x <- x[, x$subject %in% names(which(table(x$subject) == 5))]
x$subject <- droplevels(x$subject)
pca <- prcomp(t(log1p(assay(x, "counts"))))$x[, 1:3]
reducedDim(x, "PCA") <- pca
n_subjects <- nlevels(x$subject)
ord <- order(x$subject, x$time)

# test -------------------------------------------------------------------------
expect_error(plotPath(tse, "time", "subject", dimred = "PCoA",
        centroids = 1), "'centroids'")
